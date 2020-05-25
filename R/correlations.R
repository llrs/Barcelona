library("ggplot2")
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

tab <- read.delim("data/Partek_Michigan3_Kraken_Classified_genus.tsv", check.names = FALSE)
colnames(tab) <- gsub("_S.*", "", colnames(tab)) # Remove trailing numbers
counts <- tab[, -1]
genus <- tab[, 1, FALSE]
write.csv(genus, "data/genus.csv")

# From the QC step
meta <- readRDS("data_out/info_samples.RDS")
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


# Filter the samples
rna2 <- rna[, colnames(rna) %in% meta$Original]
meta2 <- droplevels(meta[meta$Original %in% colnames(rna2), ])
OTUs2 <- otus[, colnames(otus) %in% meta2$Name]

colnames(OTUs2) <- meta2$Original[match(colnames(OTUs2), meta2$Name)]

# Reorder samples to match!
meta2 <- meta2[match(colnames(rna2), meta2$Original), ]
OTUs2 <- OTUs2[match(colnames(rna2), colnames(OTUs2))]

OTUs2 <- as.matrix(OTUs2)
rna2 <- as.matrix(rna2)

abundance <- 0.005 # 0.5%
a <- prop.table(OTUs2, 2)
b <- rowSums(a > abundance)

OTUs2 <- norm_RNAseq(OTUs2)
rna2 <- norm_RNAseq(rna2)
rna2 <- filter_RNAseq(rna2)

# Filter by model ####
model <- models2 <- readRDS("data_out/models2.RDS")
model <- model[[3]]
w_rna <- model$a[[1]][, 1]
names_rna <- names(w_rna)[w_rna != 0]
w_dna <- model$a[[2]][, 1]


fOTUS2 <- OTUs2[w_dna != 0 & b != 0, ]
frna2 <- rna2[rownames(rna2) %in% names_rna, ]


# Fix names ####
# Gene names instead of ENSEMBL
s <- mapIds(org.Hs.eg.db, keys = trimVer(rownames(frna2)), keytype = "ENSEMBL", column = "SYMBOL")
frna2 <- frna2[!is.na(s), ]
rownames(frna2) <- s[!is.na(s)]
# Tax genus instead of numbers
rownames(fOTUS2) <- as.character(genus[w_dna != 0  & b != 0, 1])


# Correlations ####
df <- expand.grid(genus = rownames(fOTUS2), genes = rownames(frna2),
                  stringsAsFactors = FALSE)
df$r <- 0
df$p.value <- 1
# Remove big outliers
qrna <- quantile(frna2, c(0.05, 0.95))
qdna <- quantile(fOTUS2, c(0.05, 0.95))

rm_outliers <- function(x, quantiles) {
  x[x < quantiles[1]] <- NA
  x[x > quantiles[2]] <- NA
  x[x == 0] <- NA
  x
}

less_precision <- function(x) {
  var(x, na.rm = TRUE) < .Machine$double.eps
}

outliers <- function(x) {
  # Filter based on NA values on NAs present and amount of values
  if (sum(!is.na(x)) < 3 | sum(!is.na(x))/length(x) < 0.15 | less_precision(x)) {
    TRUE
  } else {
    FALSE
  }
}

pdf("Figures/correlations_genus.pdf")
for (i in seq_len(nrow(df))) {
  genus <- df$genus[i]
  gene <- df$genes[i]


  x <- frna2[gene, ]
  x <- rm_outliers(x, qrna)
  if (outliers(x)) {
    next
  }

  y <- fOTUS2[genus, ]
  y <- rm_outliers(y, qdna)
  if (outliers(y)) {
    next
  }
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

      plot(x, y, xlab = gene, ylab = genus,
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
saveRDS(df, "data_out/correlations.RDS")

#  Not workking because n < length p.values
# df$fdr <- p.adjust(df$p.value, method = "fdr",
                   # n = max(c(nrow(frna2), nrow(fOTUS2))))

# Plot distributions ####
df %>%
  group_by(genus) %>%
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

# Select a threshold ####
q <- quantile(abs(df$r[df$p.value < 0.05 & !is.na(df$p.value)]), 0.99)
sum(abs(df$r) > q & !is.na(df$p.value))

# Redo just those correlations above the threshold to be able to check that they are fit
subDF <- df[abs(df$r) > q & !is.na(df$p.value), ]
subDF <- subDF[order(subDF$p.value, decreasing = FALSE), c("genes", "genus")]
meta <- readRDS("data_out/refined_meta.RDS")
write.csv(meta, "data_out/refined_meta.csv", na = "", row.names = FALSE)

subMeta <- meta[match(colnames(frna2), meta$Original), ]
loc <- ifelse(meta$Exact_location == "ileum", "ileum", "colon")
loc <- factor(loc, levels = c("colon", "ileum"))
ibd <- as.character(subMeta$IBD)
ibd[is.na(ibd)] <- "C"
ibd <- factor(ibd, levels = c("C", "CD", "UC"))

pdf("Figures/high_correlations_genus.pdf")
for (i in seq_len(nrow(subDF))) {
  genus <- subDF$genus[i]
  gene <- subDF$genes[i]


  x <- frna2[gene, ]
  x <- rm_outliers(x, qrna)
  y <- fOTUS2[genus, ]
  y <- rm_outliers(y, qdna)
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
    plot(g, o, xlab = gene, ylab = genus, pch = 14+as.numeric(loc2),
         col = ibd2,
         main = paste0("Correlation: ",
                       round(co$estimate, 4), "\n",
                       "p-value: ", round(co$p.value, 4)))
    legend("top",
           legend = levels(ibd),
           fill = as.factor(levels(ibd)))
    legend("right",
           legend = levels(loc),
           pch = 14 + as.numeric(as.factor(levels(loc))))
  },
  silent = FALSE
  )
}
dev.off()
subDF <- df[abs(df$r) > q & !is.na(df$p.value), ]
write.csv(subDF, "data_out/high_correlations_genus.csv",
          na = "", row.names = FALSE)

subDF %>%
  count(genus, sort = TRUE)
subDF %>%
  count(genes, sort = TRUE)
