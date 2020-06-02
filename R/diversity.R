# Calculate the alpha diversity
library("tidyr")
library("ggplot2")
library("dplyr")
library("purrr")
library("vegan")
library("stringr")
library("phyloseq")
library("metagenomeSeq")
library("limma")
library("ggedit")
library("readxl")
library("integration")
library("lubridate")
library("org.Hs.eg.db")
library("forcats")

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
genus <- read.csv("data/genus.csv", row.names = 1)

meta$ileum <- ifelse(meta$Exact_location == "ileum", "Ileum", "colon")
meta$ileum[is.na(meta$ileum )] <- "colon"
meta$SEX[is.na(meta$SEX) | meta$SEX == ""] <- "female"
meta$Time[is.na(meta$Time)] <- "C"

rownames(meta) <- colnames(OTUs2)
otus <- as.matrix(OTUs2)
phyloseq <- phyloseq(otu_table(OTUs2, taxa_are_rows = TRUE),
              sample_data(meta),
              tax_table(genus))


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
  theme_minimal() +
  labs(x = element_blank(), y = "%") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank())


alpha_meas <- c("Simpson", "Shannon")
theme_set(theme_minimal())
richness <- estimate_richness(phyloseq)
richness$Sample <- meta$Original
richness <- merge(richness, meta, by.x = "Sample", by.y = "Original")
ggplot(richness) +
  geom_col(aes(Shannon, fct_reorder(Sample, IBD), fill = IBD, col = IBD)) +
  labs(y = "Sample") +
  theme(panel.grid.major.y = element_blank())
ggplot(richness) +
  geom_jitter(aes(x = SEX,  y = Shannon, col = IBD, shape = IBD))

p <- plot_richness(phyloseq, "SEX", "IBD", measures = alpha_meas)
remove_geom(p, 'point', 1) + geom_jitter()
ggsave("Figures/alpha_simpson_sex_ibd.png")
q <- plot_richness(phyloseq, "Time", "IBD", measures = alpha_meas)
remove_geom(q, 'point', 1) + geom_jitter()
ggsave("Figures/alpha_simpson_time_ibd.png")
r <- plot_richness(phyloseq, "IBD", measures = alpha_meas)
remove_geom(r, 'point', 1) + geom_jitter()
s <- plot_richness(phyloseq, "Exact_location", "IBD", measures = alpha_meas)
remove_geom(s, 'point', 1) + geom_jitter()
ggsave("Figures/alpha_simpson_location_ibd.png")
s <- plot_richness(phyloseq, "ileum", "IBD", measures = alpha_meas)
remove_geom(s, 'point', 1) + geom_jitter()
ggsave("Figures/alpha_simpson_location_ibd.png")
u <- plot_richness(phyloseq, "ANTITNF_responder", "IBD", measures = alpha_meas)
remove_geom(u, 'point', 1) + geom_jitter()
ggsave("Figures/alpha_simpson_responders_ibd.png")

# Beta diversity ####
# Remove empty lines
e <- apply(otus, 1, function(x){sum(x == 0)})
beta <- vegdist(t(otus[e != 194, ]), method = "jaccard")
cmd <- cmdscale(d = beta)
png("Figures/beta_jaccard_sex.png")
plot(cmd, col = as.factor(meta$SEX), main = "Beta diversity",
     xlab = "Dimension 1", ylab = "Dimension 2")
legend("topleft", legend = c("male", "female"), fill = c("red", "black"))
dev.off()
png("Figures/beta_jaccard_location.png")
plot(cmd, col = as.factor(meta$ileum), main = "Beta diversity",
     xlab = "Dimension 1", ylab = "Dimension 2")
legend("topleft", legend = c("ileum", "colon"), fill = c("red", "black"))
dev.off()
png("Figures/beta_jaccard_location.png")
plot(cmd, col = as.factor(meta$IBD), main = "Beta diversity",
     xlab = "Dimension 1", ylab = "Dimension 2")
legend("topleft", legend = c("CD", "C", "UC"),
       fill = c("red", "black", "green"))
dev.off()

# metagenomeSeq ####
MR <- phyloseq_to_metagenomeSeq(phyloseq) # For testing and comparing data
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
settings <- zigControl(maxit = 10, verbose = FALSE)
mod <- model.matrix(~ 0 + IBD + SEX, data = pData(MRtrim))
Time <- as.numeric(pData(MRtrim)$Time)
Time[is.na(Time)] <- 0
Ileum <- ifelse(pData(MRtrim)$Exact_location == "ileum", 1, 0)
mod <- cbind(mod, Time = Time, Ileum = Ileum)
fit <- fitZig(
  obj = MRtrim, mod = mod, useCSSoffset = FALSE,
  control = settings,
  block =  pData(MRtrim)$ID
)
zigFit <- slot(fit, "fit")
finalMod <- slot(fit, "fit")$design
contrast.matrix <- makeContrasts(CD = IBDCD - IBDC,
                                 UC = IBDUC - IBDC,
                                 Ileum_vs_Colon = Ileum,
                                 Male_vs_femal = SEXmale,
                                 levels = finalMod)
fit2 <- contrasts.fit(zigFit, contrast.matrix)
fit2 <- eBayes(fit2)
dt <- decideTests(fit2)
summary(dt)
dt <- as.matrix(dt)
dt[dt < 0 ] <- 1
library("UpSetR")
upset(as.data.frame(dt), keep.order = FALSE, order.by = "freq", nsets = 50)

CD <- topTreat(fit2, coef = "CD", number = Inf)
UC <- topTreat(fit2, coef = "UC", number = Inf)
Loc <- topTreat(fit2, coef = "Ileum_vs_Colon", number = Inf)
sex <- topTreat(fit2, coef = "Male_vs_femal", number = Inf)
