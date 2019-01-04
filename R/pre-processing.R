library("dplyr")
library("metagenomeSeq")
library("integration")
library("stringr")
library("purrr")

tab <- read.delim("data/Partek_Michigan3_Kraken_Classified_genus.tsv", check.names = FALSE)
colnames(tab) <- gsub("_S.*", "", colnames(tab)) # Remove trailing numbers
counts <- tab[, -1]
genus <- tab[, 1, FALSE]

# From the QC step
meta <- readRDS("info_samples.RDS")
meta$Counts <- colSums(counts)

# filter (keep in mind that it should be on the same order)
if (!all(colnames(counts) == meta$Name)) {
  stop("Reorder the samples to match the condition!!")
}
bcn <- counts[, meta$Study %in% c("BCN", "Controls")]
meta <- meta[meta$Study %in% c("BCN", "Controls"), ]

# Remove duplicate samples
replicates <- table(meta$Original)
replicate_samples <- meta[meta$Original %in% names(replicates[replicates > 1]), ]

## By keeping those more sequenced
keepDup <- replicate_samples %>%
  group_by(Original) %>%
  filter(Counts == max(Counts)) %>%
  arrange(desc(abs(Counts))) %>%
  ungroup()

nam <- c(names(replicates[replicates == 1]), keepDup$Name)
otus <- bcn[, colnames(bcn) %in% nam]
meta <- meta[meta$Name %in% nam, ]

# normalize 16S
OTUs <- norm_otus(otus, genus)


# Working with RNAseq
conn <- gzfile("data/voom.RNAseq.data.all.cal.noduplications.tsv.gz")
rna <- read.table(conn, sep = "\t", check.names = FALSE)
close(conn)

colnames(rna) <- gsub(" reseq$", "", colnames(rna))
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

# Filter them
rna2 <- rna[, colnames(rna) %in% meta$Original]
meta2 <- meta[meta$Original %in% colnames(rna2), ]
OTUs2 <- OTUs[, colnames(OTUs) %in% meta2$Name]

A <- list("RNAseq" = t(rna2), "Micro" = t(OTUs2), "Meta" = meta2)
A[1:2] <- clean_unvariable(A[1:2]) # Just the numeric ones
saveRDS(A, "data/RGCCA_data.RDS")
