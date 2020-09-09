# load data ####
library("readxl")
library("dplyr")
library("metagenomeSeq")
library("integration")
library("stringr")
library("purrr")
library("ggplot2")
library("lubridate")
library("org.Hs.eg.db")
{
tab <- read.delim("data/Partek_Michigan3_Kraken_Classified_family.tsv", check.names = FALSE)
colnames(tab) <- gsub("_S.*", "", colnames(tab)) # Remove trailing numbers
counts <- tab[, -1]
family_all <- tab[, 1, FALSE]
rownames(counts) <- family_all[, 1]
write.csv(family_all, "data/family.csv")

# From the QC step
meta <- readRDS("data_out/info_samples.RDS")
meta$Counts <- colSums(counts)

# filter (keep in mind that it should be on the same order)
if (!all(colnames(counts) == meta$Name)) {
  stop("Reorder the samples to match the condition!!")
}
bcn <- counts[, meta$Study %in% c("BCN", "Controls")]
meta <- meta[meta$Study %in% c("BCN", "Controls"), ]
m <- readRDS("data_out/refined_meta_all.RDS")

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
meta3 <- readRDS("data_out/refined_meta.RDS")
meta3$IBD <- as.character(meta3$IBD)
meta3$IBD[is.na(meta3$IBD)] <- "C"
stopifnot(all(meta3$Original == meta2$Original))
meta2 <- meta3
OTUs2 <- OTUs2[match(colnames(rna2), colnames(OTUs2))]
rownames(meta2) <- 1:nrow(meta2)

OTUs2 <- as.matrix(OTUs2)
rna2 <- as.matrix(rna2)
}
# * ileum ####
{
keep_ileum <- meta2$Exact_location %in% "ileum"
rna_ileum <- rna2[, keep_ileum]
otus_ileum <- OTUs2[, keep_ileum]

abundance <- 0.005 # 0.5%
ab <- prop.table(otus_ileum, 2)
b_ileum <- rowSums(ab > abundance)

otus_ileum <- norm_RNAseq(otus_ileum)
rna_ileum <- norm_RNAseq(rna_ileum)
rna_ileum <- filter_RNAseq(rna_ileum)
}
{
# * colon ####
keep_colon <- !keep_ileum
keep_uc <- meta2$IBD %in% "UC"
keep_c <- meta2$IBD %in% "C"
rna_colon <- rna2[, keep_colon]
otus_colon <- OTUs2[, keep_colon]

abundance <- 0.005 # 0.5%
ab <- prop.table(otus_colon, 2)
b_colon <- rowSums(ab > abundance)

otus_colon <- norm_RNAseq(otus_colon)
rna_colon <- norm_RNAseq(rna_colon)
rna_colon <- filter_RNAseq(rna_colon)
}
{
# ** Colon CD C ####
rna_colon_CD <- rna2[, keep_colon  & !keep_uc]
otus_colon_CD <- OTUs2[, keep_colon  & !keep_uc]

abundance <- 0.005 # 0.5%
ab <- prop.table(otus_colon_CD, 2)
b_colon_CD <- rowSums(ab > abundance)

otus_colon_CD <- norm_RNAseq(otus_colon_CD)
rna_colon_CD <- norm_RNAseq(rna_colon_CD)
rna_colon_CD <- filter_RNAseq(rna_colon_CD)
}
{
# ** Colon UC C ####
rna_colon_UC <- rna2[, keep_colon & (keep_uc | keep_c)]
otus_colon_UC <- OTUs2[, keep_colon & (keep_uc | keep_c)]

abundance <- 0.005 # 0.5%
ab <- prop.table(otus_colon_UC, 2)
b_colon_UC <- rowSums(ab > abundance)

otus_colon_UC <- norm_RNAseq(otus_colon_UC)
rna_colon_UC <- norm_RNAseq(rna_colon_UC)
rna_colon_UC <- filter_RNAseq(rna_colon_UC)
}
{
# * All ####
abundance <- 0.005 # 0.5%
a <- prop.table(OTUs2, 2)
b_all <- rowSums(a > abundance)

OTUs_all <- norm_RNAseq(OTUs2)
rna_all <- norm_RNAseq(rna2)
rna_all <- filter_RNAseq(rna_all)
}
{
# Filter by model ####
model <- readRDS("data_out/model2b_sgcca.RDS")
w_rna <- model$a[[1]][, 1]
names_rna <- names(w_rna)[w_rna != 0]
w_dna <- model$a[[2]][, 1]

}
# Select options ####
{
sOTUs2 <- OTUs_all
srna2 <- rna_all
b <- b_all
header <- "20200908_all_family_"

fOTUS2 <- sOTUs2[b != 0, ]
frna2 <- srna2[rownames(srna2) %in% names_rna, ]

# Fix names ####
# Gene names instead of ENSEMBL
s <- mapIds(org.Hs.eg.db, keys = trimVer(rownames(frna2)),
            keytype = "ENSEMBL", column = "SYMBOL")
frna2 <- frna2[!is.na(s), ]
rownames(frna2) <- s[!is.na(s)]
# Tax family instead of numbers
rownames(fOTUS2) <- as.character(family_all[b != 0, 1])

# Correlations ####
df <- expand.grid(family = rownames(fOTUS2), genes = rownames(frna2),
                  stringsAsFactors = FALSE)
df$r <- 0
df$p.value <- 1
df <- arrange(df, family)

pdf(paste0("Figures/", header, "correlations.pdf"))
for (i in seq_len(nrow(df))) {
  family <- df$family[i]
  gene <- df$genes[i]

  x <- frna2[gene, ]
  y <- fOTUS2[family, ]

  names(x) <- NULL
  names(y) <- NULL

  # Filter based on the number of pairwise values existing
  d <- rbind(x, y)
  pairwise <- apply(d, 2, function(z){all(!is.na(z))})
  k <- sum(pairwise, na.rm = TRUE)
  if (k/length(pairwise) < 0.15 & k < 4) {
    next
  }

  try({
    co <- cor.test(x, y, use = "spearman",
                   use = "pairwise.complete.obs")
    df$p.value[i] <- co$p.value
    df$r[i] <- co$estimate
    if (co$p.value < 0.05) {
      plot(x, y, xlab = gene, ylab = family,
           main = paste0("Correlation: ",
                         round(co$estimate, 4), "\n",
                         "p-value: ", round(co$p.value, 4)))
    }
  },
  silent = TRUE
  )
}
dev.off()
df <- df[!is.na(df$r), ]
saveRDS(df, paste0("data_out/", header, "correlations.RDS"))

# As decided on July remove some after manually inspecting the above files.
skip_micro <- c("Tsukamurellaceae", "Cyclobacteriaceae", "Beutenbergiaceae",
                "Conexibacteriaceae", "Dermacoccaceae", "Nocardiaceae",
                "Thermaceae", "Thermomicrobiaceae", "Beijerinckiaceae",
                "Campylobacteraceae", "Halanaerobiaceae")
skip_gene <- c("GIMD1")

# * Plot distributions ####
df %>%
  group_by(family) %>%
  summarise(n = sum(p.value < 0.05)) %>%
  filter(n != 0) %>%
  arrange(desc(n)) %>%
  ggplot() +
  geom_histogram(aes(n)) +
  theme_minimal() +
  labs(x = element_blank(), title = "Significant correlations by microorganism", y = "n")

df %>%
  group_by(genes) %>%
  summarise(n = sum(p.value < 0.05)) %>%
  filter(n != 0) %>%
  arrange(desc(n)) %>%
  ggplot() +
  geom_histogram(aes(n)) +
  theme_minimal() +
  labs(x = element_blank(), title = "Significant correlations by gene", y = "n")
df %>%
  filter(p.value < 0.05) %>%
  ggplot() +
  geom_histogram(aes(abs(r)), binwidth = 0.005) +
  scale_y_log10() +
  theme_minimal() +
  labs(title = "Distribution of significant correlations",
       x = "Correlation (absolute value)", y = "n")

# ** Higher correlations ####
q <- quantile(abs(df$r[df$p.value < 0.05 & !is.na(df$p.value)]), 0.99)
sum(abs(df$r) > q & !is.na(df$p.value))

# Redo just those correlations above the threshold to be able to check that they are fit
subDF <- df[abs(df$r) > q & !is.na(df$p.value), ]
subDF <- subDF[!subDF$family %in% skip_micro & !subDF$genes %in% skip_gene, ]
subDF <- subDF[order(subDF$family, subDF$p.value, decreasing = FALSE), ]
# write.csv(meta, "data_out/20200629_refined_meta.csv", na = "", row.names = FALSE)
# write.csv(subDF, header, na = "", row.names = FALSE)
write.csv(subDF, paste0("data_out/", header, "high_correlations_family.csv"),
          na = "", row.names = FALSE)

subMeta <- meta2[match(colnames(frna2), meta2$Original), ]
loc <- ifelse(subMeta$Exact_location == "ileum", "ileum", "colon")
loc <- factor(loc, levels = c("colon", "ileum"))
ibd <- as.character(subMeta$IBD)
ibd[is.na(ibd)] <- "C"
ibd <- factor(ibd, levels = c("C", "CD", "UC"))
act <- subMeta$Activity
act[is.na(act)] <- "INACTIVE"
pch <- ifelse(act == "ACTIVE", 15, 0) + as.numeric(loc) -1

pdf(paste0("Figures/", header, "high_correlations_family.pdf"))
for (i in seq_len(nrow(subDF))) {
  family <- subDF$family[i]
  gene <- subDF$genes[i]

  x <- frna2[gene, ]
  # x <- rm_outliers(x, qrna)
  y <- fOTUS2[family, ]
  # y <- rm_outliers(y, qdna)
  names(x) <- NULL
  names(y) <- NULL
  rm_samples <- !is.na(x) & !is.na(y)
  loc2 <- loc[rm_samples]
  ibd2 <- ibd[rm_samples]
  g <- x[rm_samples]
  o <- y[rm_samples]

  if (var(g, use = "everything") < .Machine$double.eps) {
    next
  }
  if (var(o, use = "everything") < .Machine$double.eps) {
    next
  }

  try({
    co <- cor.test(g, o, use = "spearman",
                   use = "pairwise.complete.obs")
    if (co$p.value >= 0.05) {
      next
    }
    plot(g, o, xlab = gene, ylab = family, pch = pch,
         col = ibd2,
         main = paste0("Correlation: ",
                       round(co$estimate, 4), "\n",
                       "p-value: ", round(co$p.value, 4)))
    legend("top",
           legend = levels(ibd),
           fill = as.factor(levels(ibd)))
    legend("right",
           legend = c("Inactive colon", "Inactive ileum",
                      "Active colon", "Active ileum"),
           pch = c(0, 1, 15, 16))
  },
  silent = FALSE
  )
}
dev.off()
}

subDF %>%
  count(family, sort = TRUE) %>%
  head()
subDF %>%
  count(genes, sort = TRUE) %>%
  head()
