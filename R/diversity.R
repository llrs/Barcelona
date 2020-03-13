# Calculate the alpha diversity
library("tidyr")
library("ggplot2")
library("dplyr")
library("purrr")
library("vegan")
library("stringr")
library("phyloseq")
library("ggedit")


tab <- read.delim("data/Partek_Michigan3_Kraken_Classified_genus.tsv",
                  check.names = FALSE)
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
conn <- gzfile("data/voom.RNAseq.data.all.cal.noduplications.tsv.gz")
# conn <- gzfile("data/TNF.all.samples.original.counts.tsv.gz") # TODO See if this is a good choice
rna <- read.table(conn, sep = "\t", check.names = FALSE)

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


# Filter the samples
rna2 <- rna[, colnames(rna) %in% meta$Original]
meta2 <- droplevels(meta[meta$Original %in% colnames(rna2), ])
OTUs2 <- otus[, colnames(otus) %in% meta2$Name]

colnames(OTUs2) <- meta2$Original[match(colnames(OTUs2), meta2$Name)]

# Reorder samples to match!
meta2 <- meta2[match(colnames(rna2), meta2$Original), ]
OTUs2 <- OTUs2[match(colnames(rna2), colnames(OTUs2))]

A <- readRDS("data/RGCCA_data.RDS")
otus <- OTUs2
genus <- read.csv("data/genus.csv", row.names = 1)

meta <- A$Meta
levels(meta$IBD) <- c(levels(meta$IBD), "C")
meta$IBD[is.na(meta$IBD)] <- "C"
rownames(meta) <- colnames(otus)
phyloseq <- phyloseq(otu_table(otus, taxa_are_rows = TRUE),
              sample_data(meta),
              tax_table(genus))


alpha_meas <- c("Simpson")
theme_set(theme_minimal())
p <- plot_richness(o, "SEX", "IBD", measures = alpha_meas)
remove_geom(p, 'point', 1) + geom_jitter()
q <- plot_richness(o, "Time", "IBD", measures = alpha_meas)
remove_geom(q, 'point', 1) + geom_jitter()
r <- plot_richness(o, "IBD", measures = alpha_meas)
remove_geom(r, 'point', 1) + geom_jitter()
s <- plot_richness(o, "Exact_location", "IBD", measures = alpha_meas)
remove_geom(s, 'point', 1) + geom_jitter()
u <- plot_richness(o, "ANTITNF_responder", "IBD", measures = alpha_meas)
remove_geom(u, 'point', 1) + geom_jitter()


d <- distance(o, "jaccard")


# Alpha diversity ####
alpha <- prop.table(otus, 2)*100
a <- as.data.frame(alpha)
a$otus <- genus[, 1]
a <- pivot_longer(a, colnames(alpha))

b <- a %>%
  group_by(name) %>%
  mutate(cumsum = cumsum(value),
         otus = forcats::fct_lump(otus, n = 20, w = value),
         n = 1:n()) %>%
  ungroup() %>%
  arrange(name)


  ggplot(b) +
  geom_col(aes(name, value, fill = otus), col = "transparent") +
  guides(fill = FALSE) +
  theme_minimal()

# Beta diversity ####
beta <- vegdist(otus, method = "jaccard")
cmd <- cmdscale(d = beta)
plot(cmd, col = as.factor(A$Meta$SEX))

# metagenomeSeq ####
MR <- phyloseq_to_metagenomeSeq(o) # For testing and comparing data
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
settings <- zigControl(maxit = 10, verbose = TRUE)
mod <- model.matrix(~ 0 + IBD, data = pData(MRtrim))
mod <- cbind(mod, ID = pData(MRtrim)$ID, Time = as.numeric(pData(MRtrim)$Time))
fit <- fitZig(
  obj = MRtrim, mod = mod, useCSSoffset = FALSE,
  control = settings
)

eb <- fit$eb
