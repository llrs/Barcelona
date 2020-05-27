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

# From the QC step
meta <- readRDS("data_out/refined_meta_all.RDS")
ccounts <- colSums(counts)
meta <- meta[meta$Name %in% colnames(counts), ]

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

# Filter the samples
OTUs2 <- otus[, colnames(otus) %in% meta$Name]

colnames(OTUs2) <- meta$Original[match(colnames(OTUs2), meta$Name)]

# Reorder samples to match!
OTUs2 <- OTUs2[, match(meta$Original, colnames(OTUs2))]
stopifnot(all(colnames(OTUs2) == meta$Original))
OTUs2 <- as.matrix(OTUs2)

rownames(meta) <- colnames(OTUs2)
fd <- AnnotatedDataFrame(genus)
# Create the object
MR <- newMRexperiment(
  OTUs2,
  phenoData = AnnotatedDataFrame(meta),
  featureData = fd
)

filterData(MR, present = 10, depth = 1000)
p <- cumNormStatFast(MR)
MR <- cumNorm(MR, p = p)
assayData(MR)$relative <- assayData(MR)$counts / rowSums(assayData(MR)$counts)
assayData(MR)$prevalence <- as.matrix(assayData(MR)$counts != 0)

rareFeatures <- which(rowSums(MRcounts(MR) > 0) < 10)
MRtrim <- MR[-rareFeatures, ]
MRp <- cumNormStat(MRtrim, pFlag = TRUE, main = "Trimmed lung data")
MRtrim <- cumNorm(MRtrim, p = MRp)

normFactor <- normFactors(MRtrim)
assayData(MRtrim)$relative <- assayData(MRtrim)$counts / rowSums(assayData(MRtrim)$counts)
assayData(MRtrim)$prevalence <- as.matrix(assayData(MRtrim)$counts != 0)

normFactor <- log2(normFactor / median(normFactor) + 1)

# Set up variables for the model
ts <- pData(MRtrim)$Time
ts[is.na(ts)] <- "0"
mod <- data.frame(Time = paste0("t", ts))
pData(MRtrim)$ileum <- ifelse(meta$Exact_location == "ileum", "Ileum", "colon")
pData(MRtrim)$ileum[is.na(pData(MRtrim)$ileum )] <- "colon" # Checked manually on the database
ibd <- as.character(pData(MRtrim)$IBD)
ibd[is.na(ibd)] <- "C"
ibd[ibd == "CONTROL"] <- "C"
ibd <- factor(ibd, levels = c("C", "CD", "UC"))
mod <- model.matrix(~ts+pData(MRtrim)$ileum+ibd)
colnames(mod) <- c("(Intercept)", "ts14", "ts46", "Ileum", "CD","UC")

# Comparisions ####
fit <- fitZig(
  obj = MRtrim, mod = mod, useCSSoffset = FALSE,
  block = pData(MRtrim)$ID,
  control = zigControl(maxit = 1, verbose = FALSE)
)

zigFit <- slot(fit, "fit")

contrast.matrix <- makeContrasts(ftw = "ts46-ts14",
                                 ftf = "ts14-ts46",
                                 fw = "ts14-ts46",
                                 ileum_vs_colon = "(Intercept-Ileum)-Ileum",
                                 C_vs_CD = "UC-CD",
                                 C_vs_UC = "CD-UC",
                                 levels = mod)

apply(mod %*% contrast.matrix, 2, table)
fit2 <- contrasts.fit(zigFit, contrast.matrix)
fit2 <- eBayes(fit2)

dt <- decideTests(fit2)
summary(dt)
microb <- dt[, 4, drop = TRUE] != 0

topTable(fit2, coef = "ileum_vs_colon", number = 43)

# Prevalence ####

rownames(OTUs2) <- genus[, 1]
meta$ileum <- ifelse(meta$Exact_location == "ileum", "Ileum", "colon")
meta$Time[is.na(meta$Time)] <- "C"

gen <- comb_prevalence(OTUs2, meta, c("Time"))
gen_se <- comb_prevalence(OTUs2, meta, c("Time", "SEX"))
tre <- comb_prevalence(OTUs2, meta, c("treatment"))
se <- comb_prevalence(OTUs2, meta, c("SEX"))
std <- comb_prevalence(OTUs2, meta, c("Study"))
ibdm <- comb_prevalence(OTUs2, meta, c("IBD", "Time", "ileum"))
ibd <- comb_prevalence(OTUs2, meta, c("IBD", "ileum"))
loc <- comb_prevalence(OTUs2, meta, c("ileum"))
