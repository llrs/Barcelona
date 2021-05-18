# Prepare macros to explore raw data
library("dplyr")
library("purrr")
library("stringr")
library("org.Hs.eg.db")
library("integration")

seqtab.nochim <- readRDS("data/ASV.RDS")

ASV <- colnames(seqtab.nochim)
counts_ASV <- seqtab.nochim
colnames(counts_ASV) <- NULL

tab <- t(counts_ASV)
# Remove trailing numbers
colnames(tab) <- gsub("_S.*", "", colnames(tab))
colnames(tab) <- gsub("_p.*", "", colnames(tab))

microorganism <-  readRDS("output/taxonomy_ASV.RDS")$tax
microorganism <- cbind(microorganism, "ASV" = rownames(microorganism))
microorganism <- as.data.frame(microorganism, stringsAsFactors = FALSE)
microorganism$rowname <- as.character(seq_len(nrow(microorganism)))
rownames(microorganism) <- as.character(seq_len(nrow(microorganism)))

# From the QC step
meta <- readRDS("output/info_samples.RDS")
meta <- meta[match(colnames(tab), meta$Name), ]

# filter (keep in mind that it should be on the same order)
if (!all(colnames(tab) == meta$Name)) {
  stop("Reorder the samples to match the condition!!")
}
bcn <- tab[, meta$Study %in% c("BCN", "Controls")]
meta <- meta[meta$Study %in% c("BCN", "Controls"), ]
meta$Counts <- colSums(tab)

# Remove duplicate samples
replicates <- table(meta$Original)
replicate_samples <- meta[meta$Original %in% names(replicates[replicates > 1]), ]

## But keeping those more sequenced
keepDup <- replicate_samples %>%
  group_by(Original) %>%
  filter(Counts == max(Counts)) %>%
  arrange(desc(abs(Counts))) %>%
  ungroup()

nam <- c(names(replicates[replicates == 1]), keepDup$Name)
otus <- bcn[, colnames(bcn) %in% nam]

# Working with RNAseq
# From Juanjo: The original counts are ok, but I need to remove the reseq samples as they
# have different length and bias the PCA
conn <- gzfile("data/TNF.all.samples.original.counts.tsv.gz")
rna <- read.table(conn, sep = "\t", check.names = FALSE)

rna <- rna[ , !grepl(" reseq$", colnames(rna))] # Remove as said
colnames(rna)[grep("[Ww]", colnames(rna))] <- tolower(colnames(rna)[grep("[Ww]", colnames(rna))])

correct_bcn <- function(x) {
  if (length(x) > 1) {
    a <- str_pad(x[1], width = 3, pad = "0")
    b <- str_pad(x[2], width = 3, pad = "0")
    x <- paste(a, b, sep = "-w")
  }
  x
}

colnames2 <- colnames(rna) %>%
  str_split("-w") %>% # Ready for BCN
  map(correct_bcn) %>%
  unlist() %>%
  gsub("-T-TR-", "-T-DM-", .) # Ready for TRIM
colnames(rna) <- colnames2


# Filter the samples ####
rna2 <- rna[, colnames(rna) %in% meta$Original]
meta2 <- droplevels(meta[meta$Original %in% colnames(rna2), ])
OTUs2 <- otus[, colnames(otus) %in% meta2$Name]

colnames(OTUs2) <- meta2$Original[match(colnames(OTUs2), meta2$Name)]

# Reorder samples to match!
meta2 <- meta2[match(colnames(rna2), meta2$Original), ]
meta3 <- readRDS("output/refined_meta_wo_out.RDS")
meta3$IBD <- as.character(meta3$IBD)
meta3$IBD[meta3$IBD == "CONTROL"] <- "C"
meta3$sample_location[meta3$IBD == "CONTROL"] <- "colon"
# If using wo outliers we removed 2 samples
# stopifnot(all(meta3$Original == meta2$Original))

A <- readRDS("data/RGCCA_data_wo_out.RDS")

meta <- A$Meta %>%
  dplyr::select(Original:UC_endoscopic_remission) %>% # Remove processing columns
  mutate(sample_location = case_when(Exact_location == "ileum" ~ "Ileum",
                           !is.na(Exact_location) ~ "Colon"),
         IBD = case_when(is.na(IBD) ~ "CONTROL",
                         IBD == "UC" ~ "UC",
                         IBD == "CD" ~ "CD",
                         TRUE ~ "CONTROL"),
         inv = case_when(Exact_location == tolower(gsub(" COLON$", "", Aftected_area)) ~ "involved",
                         TRUE ~ "not_involved")) %>%
  select_if(function(x)any(!is.na(x))) # Select just those with information

# Common data ####
m2 <- meta %>%
  dplyr::select(Original, IBD, SEX, sample_location, Patient_ID, ANTITNF_responder, Age, diagTime, AgeDiag, `Phenotype CD`, Ulcers) %>%
  arrange(IBD, sample_location, SEX, Age)
rownames(m2) <- m2$Original

rna2 <- rna2[, match(m2$Original, colnames(rna2))]
OTUs2 <- OTUs2[, match(colnames(rna2), colnames(OTUs2))]


write.table(m2,
            file = "output/GETS_sample.tsv", quote = FALSE, sep = "\t", row.names = FALSE, na = " ")
df <- data.frame(type = c("SAMPLEINFO", "SAMPLEINFO", "SAMPLEINFO", "SAMPLEINFO", "SAMPLEINFO", "SAMPLEINFO", "SAMPLEINFO"),
                 column = c("IBD", "SEX", "Ulcers", "Age", "diagTime", "AgeDiag", "sample_location"),
                 value = c("", "male", "yes", "", "", "", "ileum"),
                 color = c("AUTO", "ORANGE", "RED", "GRADIENT_RED", "GRADIENT_RED", "GRADIENT_RED", "MAGENTA"))
write.table(df, file = "output/GETS_colors.tsv", quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)

# RNASeq ####


rna3 <- integration::filter_RNAseq(rna2)
out <- mapIds(org.Hs.eg.db, keys = trimVer(rownames(rna3)), column = "SYMBOL", keytype = "ENSEMBL")
out <- data.frame(ensembl = trimVer(rownames(rna3)), gene = out)


cm <- cbind(genes = rownames(rna3), integration::norm_RNAseq(rna3))
write.table(cm, file = "output/GETS_matrix_RNASeq.tsv", quote = FALSE, sep = "\t", row.names = FALSE)
system2("gzip", args = "-kf output/GETS_matrix_RNASeq.tsv") # Compress it to upload to website
write.table(out, file = "output/GETS_gene_RNASeq.tsv", quote = FALSE, sep = "\t", row.names = FALSE, na = "")
system2("gzip", args = "-kf output/GETS_gene_RNASeq.tsv") # Compress it to upload to website


# Microbiome ####

out <- microorganism
k <- apply(OTUs2, 1, var) != 0
rownames(OTUs2) <- as.character(seq_len(nrow(OTUs2)))
out <- out[k, ]
OTUs3 <- OTUs2[k, ]
cm <- cbind(genes = rownames(OTUs3), integration::norm_RNAseq(OTUs3))
write.table(cm, file = "output/GETS_matrix_ASV.tsv", quote = FALSE, sep = "\t", row.names = FALSE)
system2("gzip", args = "-kf output/GETS_matrix_ASV.tsv") # Compress it to upload to website
write.table(out, file = "output/GETS_ASV.tsv", quote = FALSE, sep = "\t", row.names = FALSE, na = "")
system2("gzip", args = "-kf output/GETS_ASV.tsv") # Compress it to upload to website
