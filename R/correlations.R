{# load data ####
library("readxl")
library("openxlsx")
library("dplyr")
library("metagenomeSeq")
library("integration")
library("stringr")
library("purrr")
library("ggplot2")
library("lubridate")
library("org.Hs.eg.db")
library("dplyr")
}
{
# ASVs
seqtab.nochim <- readRDS("data/ASV.RDS")

ASV <- colnames(seqtab.nochim)
counts_ASV <- seqtab.nochim
colnames(counts_ASV) <- NULL

tab <- t(counts_ASV)
# Remove trailing numbers
colnames(tab) <- gsub("_S.*", "", colnames(tab))
colnames(tab) <- gsub("_p.*", "", colnames(tab))

# Taxonomy of the samples
microorganism <-  readRDS("output/taxonomy_ASV.RDS")$tax
microorganism <- cbind(microorganism, "ASV" = rownames(microorganism))
microorganism <- as.data.frame(microorganism, stringsAsFactors = FALSE)
microorganism$rowname <- as.character(seq_len(nrow(microorganism)))
rownames(microorganism) <- as.character(seq_len(nrow(microorganism)))
mit <- rownames(microorganism)[microorganism[, "Family"] == "Mitochondria"]

# Collapse the data so that there is one Genus per sample, not multiple
tab2 <- tab
colnames(tab2) <- seq_along(colnames(tab))
tab2 <- cbind(tab2, microorganism)

groups <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
rm_micro <- tab2 %>%
  group_by(!!!groups) %>%
  summarize(across(where(is.numeric), sum)) %>%
  as.data.frame()
# Change the first number according to the number of columns on group_by +1
colnames(rm_micro)[(length(groups) + 1):ncol(rm_micro)] <- colnames(tab)
counts <- rm_micro[(length(groups) + 1):ncol(rm_micro)]
microorganism <- rm_micro[, 1:length(groups)] # And here the same number

# From the QC step
meta <- readRDS("output/info_samples.RDS")
colnames(counts) <- gsub("\\.[0-9]", "", colnames(counts))
meta <- meta[match(colnames(counts), meta$Name), ]

# filter (keep in mind that it should be on the same order)
if (!all(colnames(counts) == meta$Name)) {
  stop("Reorder the samples to match the condition!!")
}
bcn <- counts[, meta$Study %in% c("BCN", "Controls")]
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
meta3 <- readRDS("output/refined_meta.RDS")
meta3$IBD <- as.character(meta3$IBD)
meta3$IBD[meta3$IBD == "CONTROL"] <- "C"
# If using wo outliers we removed 2 samples
stopifnot(all(meta3$Original == meta2$Original))
meta2 <- meta3
OTUs2 <- OTUs2[, match(colnames(rna2), colnames(OTUs2))]
rownames(meta2) <- 1:nrow(meta2)

OTUs2 <- as.matrix(OTUs2)
rownames(OTUs2) <- as.character(seq_len(nrow(OTUs2)))
rna2 <- as.matrix(rna2)
rownames(rna2) <- trimVer(rownames(rna2))
}
{
# * ileum: CD + C ####
keep_ileum <- meta2$Exact_location %in% "ileum"
rna_ileum <- rna2[, keep_ileum]
otus_ileum <- OTUs2[, keep_ileum]

abundance <- 0.005 # 0.5%
ab <- prop.table(otus_ileum, 2)
b_ileum <- rowSums(ab > abundance)

otus_ileum_norm <- norm_RNAseq(otus_ileum)
rna_ileum_norm <- norm_RNAseq(rna_ileum)
rna_ileum_norm <- filter_RNAseq(rna_ileum_norm)
}
{
# * colon ####
keep_colon <- !keep_ileum
keep_uc <- meta2$IBD %in% "UC"
keep_c <- meta2$IBD %in% "C" | is.na(meta2$IBD)
rna_colon <- rna2[, keep_colon]
otus_colon <- OTUs2[, keep_colon]

abundance <- 0.005 # 0.5%
ab <- prop.table(otus_colon, 2)
b_colon <- rowSums(ab > abundance)

otus_colon_norm <- norm_RNAseq(otus_colon)
rna_colon_norm <- norm_RNAseq(rna_colon)
rna_colon_norm <- filter_RNAseq(rna_colon_norm)
}
{
# ** Colon CD + C ####
rna_colon_CD <- rna2[, keep_colon  & !keep_uc]
otus_colon_CD <- OTUs2[, keep_colon  & !keep_uc]

abundance <- 0.005 # 0.5%
ab <- prop.table(otus_colon_CD, 2)
b_colon_CD <- rowSums(ab > abundance)

otus_colon_CD_norm <- norm_RNAseq(otus_colon_CD)
rna_colon_CD_norm <- norm_RNAseq(rna_colon_CD)
rna_colon_CD_norm <- filter_RNAseq(rna_colon_CD_norm)
}
{
# ** Colon UC + C ####
rna_colon_UC <- rna2[, keep_colon & (keep_uc | keep_c)]
otus_colon_UC <- OTUs2[, keep_colon & (keep_uc | keep_c)]

abundance <- 0.005 # 0.5%
ab <- prop.table(otus_colon_UC, 2)
b_colon_UC <- rowSums(ab > abundance)

otus_colon_UC_norm <- norm_RNAseq(otus_colon_UC)
rna_colon_UC_norm <- norm_RNAseq(rna_colon_UC)
rna_colon_UC_norm <- filter_RNAseq(rna_colon_UC_norm)
}
{
# * All ####
abundance <- 0.005 # 0.5%
a <- prop.table(OTUs2, 2)
b_all <- rowSums(a > abundance)
otus_all <- OTUs2
OTUs_all_norm <- norm_RNAseq(OTUs2)

rna_all <- rna2
rna_all_norm <- norm_RNAseq(rna2)
rna_all_norm <- filter_RNAseq(rna_all_norm)
}
{
# Filter by model ####
model <- readRDS("output/model2b2_sgcca_b.RDS")
w_rna <- model$a[[1]][, 1]
names_rna <- trimVer(names(w_rna)[w_rna != 0])
w_dna <- model$a[[2]][, 1]

}

# Functions ####
outliers <- function(x, k = 3){
  q <- quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE)
  iqr <- q[3] - q[1]
  x > q[1] + k*iqr | x < q[3] - k*iqr
}

remove_na <- function(x){
  x[!is.na(x)]
}

correlations_all <- function(otus_norm, otus, rna_norm, rna, b, header,
                             meta, names_rna, families) {

  stopifnot(length(unique(meta2$IBD[meta2$Original %in% colnames(otus_norm)]))  >= 2)

  fOTUS2 <- otus_norm[b != 0, ]
  frna2 <- rna_norm[rownames(rna_norm) %in% names_rna, ]

  # Internal options to filter samples
  proportion_lost <- 0.15
  samples_correlated <- 20

  # Fix names
  # Gene names instead of ENSEMBL
  s <- mapIds(org.Hs.eg.db, keys = rownames(frna2),
              keytype = "ENSEMBL", column = "SYMBOL")
  frna2 <- frna2[!is.na(s), ]
  rownames(frna2) <- s[!is.na(s)]
  # Tax family instead of numbers
  rownames(fOTUS2) <- as.character(families[b != 0, 1])
  rownames(otus) <- as.character(families[ , 1])

  # Correlations
  df <- expand.grid(family = rownames(fOTUS2), genes = rownames(frna2),
                    stringsAsFactors = FALSE)
  df$r <- 0
  df$p.value <- 1
  df$samples <- 0
  df <- arrange(df, family)


  # * Plot correlations
  pdf(paste0("Figures/", header, "correlations.pdf"))
  for (i in seq_len(nrow(df))) {
    family <- df$family[i]
    if (is.na(family) || family == "NA" || !family %in% rownames(fOTUS2)) {
      next
    }
    gene <- df$genes[i]
    if (is.na(gene) || gene == "NA") {
      next
    }
    message("Plots family: ", family, " Gene: ", gene)
    x <- frna2[gene, ]
    y <- fOTUS2[family, ]

    names(x) <- NULL
    names(y) <- NULL

    # Filter based on that the original matrix had 0
    x_remove <- rna[names(gene), ] != 0
    y_remove <- otus[family, ] != 0
    xx <- remove_na(x[x_remove & y_remove])
    yy <- remove_na(y[x_remove & y_remove])

    stopifnot(length(xx) == length(yy))
    if (length(xx)/length(x) < proportion_lost &
        length(xx) < samples_correlated | var(xx) == 0 | var(yy) == 0) {
      next
    }
    df$samples[i] <- length(xx)
    try({
      co <- cor.test(xx, yy, use = "spearman", use = "pairwise.complete.obs")
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
  saveRDS(df, paste0("output/", header, "correlations.RDS"))

  # As decided on July remove some after manually inspecting the above files.
  # skip_micro <- c("Tsukamurellaceae", "Cyclobacteriaceae", "Beutenbergiaceae",
  #                 "Conexibacteriaceae", "Dermacoccaceae", "Nocardiaceae",
  #                 "Thermaceae", "Thermomicrobiaceae", "Beijerinckiaceae",
  #                 "Campylobacteraceae", "Halanaerobiaceae")
  skip_micro <- c("Mitochondria")
  # skip_gene <- c("GIMD1")

  # ** Higher correlations
  # On 10/09/2020 decided to remove the 0/lowest offset of the correlation,
  # but plot it anyway.
  # Also to include the controls in all the correlations

  # Redo just those correlations above the threshold to be able to check that they are fit
  subDF <- df[df$p.value < 0.05, ]
  q <- quantile(abs(subDF$r), 0.95)
  subDF <- subDF[abs(subDF$r) >= q, ]
  subDF <- subDF[!subDF$family %in% skip_micro, ] # & !subDF$genes %in% skip_gene
  subDF <- subDF[order(subDF$family, subDF$p.value, decreasing = FALSE), ]
  write.csv(subDF, paste0("output/", header, "high_correlations_family.csv"),
            na = "", row.names = FALSE)
  write.xlsx(subDF, file = paste0("output/", header, "high_correlations_family.xlsx"),
             colNames = TRUE)

  subMeta <- meta2[match(colnames(frna2), meta2$Original), ]
  loc <- ifelse(subMeta$Exact_location == "ileum", "ileum", "colon")
  loc <- factor(loc, levels = c("colon", "ileum"))
  ibd <- as.character(subMeta$IBD)
  ibd[is.na(ibd)] <- "C"
  ibd <- factor(ibd, levels = c("C", "CD", "UC"))
  act <- subMeta$Activity
  act[is.na(act)] <- "INACTIVE"
  pch <- ifelse(act == "ACTIVE", 15, 0) + as.numeric(loc) - 1

  pdf(paste0("Figures/", header, "high_correlations_family.pdf"))
  for (i in seq_len(nrow(subDF))) {
    family <- subDF$family[i]
    if (is.na(family) || family == "NA") {
      next
    }
    gene <- subDF$genes[i]
    if (is.na(gene) || gene == "NA") {
      next
    }
    message("Plots high family: ", family, " Gene: ", gene)
    x <- frna2[gene, ]
    y <- fOTUS2[family, ]
    names(x) <- NULL
    names(y) <- NULL

    # Filter based on that the original matrix had 0
    x_remove <- rna[names(gene), ] != 0
    y_remove <- otus[family, ] != 0
    keep1 <- x_remove & y_remove
    xx <- remove_na(x[keep1])
    yy <- remove_na(y[keep1])

    stopifnot(length(xx) == length(yy))
    if (length(xx)/length(x) < proportion_lost &
        length(xx) < samples_correlated | var(xx) == 0 |
        var(yy) == 0) {
      next
    }

    try({
      co <- cor.test(xx, yy, use = "spearman",
                     use = "pairwise.complete.obs")
      plot(x, y, xlab = gene, ylab = family, pch = pch,
           col = ibd,
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
  return(df)
}

# Run for each ####
rownames(otus_all) <- microorganism[, "Genus"]
rownames(OTUs_all_norm) <- microorganism[, "Genus"]
date <- format(lubridate::today(), "%Y%m%d")
df_all <- correlations_all(otus_norm = OTUs_all_norm,
                 otus = otus_all,
                 rna_norm = rna_all_norm,
                 rna = rna,
                 b = b_all,
                 header = paste0(date, "_all_genus_"),
                 meta = meta2,
                 names_rna = names_rna,
                 families = microorganism[, "Genus", drop = FALSE])
rownames(otus_colon) <- microorganism[, "Genus"]
df_colon <- correlations_all(otus_norm = otus_colon_norm,
                 otus = otus_colon,
                 rna_norm = rna_colon_norm,
                 rna = rna_colon,
                 b = b_colon,
                 header = paste0(date, "_colon_genus_"),
                 meta = meta2,
                 names_rna = names_rna,
                 families = microorganism[, "Genus", drop = FALSE])
df_colon_UC <- correlations_all(otus_norm = otus_colon_UC_norm,
                 otus = otus_colon_UC,
                 rna_norm = rna_colon_UC_norm,
                 rna = rna_colon_UC,
                 b = b_colon_UC,
                 header = paste0(date, "_colon_UC_genus_"),
                 meta = meta2,
                 names_rna = names_rna,
                 families = microorganism[, "Genus", drop = FALSE])
df_colon_CD <- correlations_all(otus_norm = otus_colon_CD_norm,
                 otus = otus_colon_CD,
                 rna_norm = rna_colon_CD_norm,
                 rna = rna_colon_CD,
                 b = b_colon_CD,
                 header = paste0(date, "_colon_CD_genus_"),
                 meta = meta2,
                 names_rna = names_rna,
                 families = microorganism[, "Genus", drop = FALSE])
df_ileum <- correlations_all(otus_norm = otus_ileum_norm,
                 otus = otus_ileum,
                 rna_norm = rna_ileum_norm,
                 rna = rna_ileum,
                 b = b_ileum,
                 header = paste0(date, "_ileum_genus_"),
                 meta = meta2,
                 names_rna = names_rna,
                 families = microorganism[, "Genus", drop = FALSE])

dev.off()
subDF <- df[abs(df$r) > q & !is.na(df$p.value), ]
write.csv(subDF, paste0("output/", header, "high_correlations_family.csv"),
          na = "", row.names = FALSE)

subDF %>%
  count(family, sort = TRUE) %>%
  head()
subDF %>%
  count(genes, sort = TRUE) %>%
  head()

