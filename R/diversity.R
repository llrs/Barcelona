# Calculate the alpha diversity
# Preparation steps ####
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

{
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
# The missing values of Exact location
meta$ileum[is.na(meta$ileum )] <- "colon"
meta$SEX[is.na(meta$SEX) | meta$SEX == ""] <- "female"
meta$Time[is.na(meta$Time)] <- "C"

rownames(meta) <- colnames(OTUs2)
otus <- as.matrix(OTUs2)
phyloseq <- phyloseq(otu_table(OTUs2, taxa_are_rows = TRUE),
              sample_data(meta),
              tax_table(genus))
}
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
{method <- "bray"
beta <- vegdist(t(otus[e != 194, ]), method = method)
}
{cmd <- cmdscale(d = beta)
stopifnot(all(rownames(cmd) == meta$Original))
#  * All samples ####
png(paste0("Figures/beta_", method, "_sex.png"))
plot(cmd, col = as.factor(meta$SEX), main = "Beta diversity",
     xlab = "Dimension 1", ylab = "Dimension 2")
legend("topleft", legend = c("male", "female"), fill = c("red", "black"))
dev.off()
png(paste0("Figures/beta_", method, "_location.png"))
plot(cmd, col = as.factor(meta$ileum), main = "Beta diversity",
     xlab = "Dimension 1", ylab = "Dimension 2")
legend("topleft", legend = c("ileum", "colon"), fill = c("red", "black"))
dev.off()
png(paste0("Figures/beta_", method, "_location.png"))
ibds <- meta$IBD
ibds[is.na(ibds)] <- "CONTROL"
plot(cmd, col = as.factor(ibds), main = "Beta diversity",
     xlab = "Dimension 1", ylab = "Dimension 2")
legend("topleft", legend = c("CD", "C", "UC"),
       fill = c("black", "red", "green"))
dev.off()}

# * Comparing beta diversity ####
# Using refined meta because we want only those that were used on the integration
{m <- readRDS("data_out/refined_meta.RDS")
m$IBD <- as.character(m$IBD)
m$IBD[is.na(m$IBD)] <- "CONTROL"
b <- as.matrix(beta)
f <- function(m, l1, l2) {
  out <- m[l1, l2]
  stopifnot(ncol(out) == length(l2))
  stopifnot(nrow(out) == length(l1))
  diag(out) <- NA
  dim(out) <- NULL
  out[!is.na(out)]
}
controls <- m$Original[m$IBD %in% "CONTROL"]
cases <- m$Original[!m$IBD %in% "CONTROL"]
UC <- m$Original[m$IBD %in% "UC"]
CD <- m$Original[m$IBD %in% "CD"]
male <- m$Original[m$SEX %in% "male"]
female <- m$Original[m$SEX %in% "female"]
ileum <- m$Original[m$Exact_location == "ileum"]
colon <- m$Original[m$Exact_location != "ileum"]

controls_beta <- f(b, controls, controls)
UC_beta <- f(b, UC, UC)
CD_beta <- f(b, CD, CD)
controls_ileum_beta <- f(b, intersect(controls, ileum),
                         intersect(controls, ileum))
UC_ileum_beta <- f(b, intersect(UC, ileum), intersect(UC, ileum))
CD_ileum_beta <- f(b, intersect(CD, ileum), intersect(CD, ileum))
controls_colon_beta <- f(b, intersect(controls, colon),
                         intersect(controls, colon))
UC_colon_beta <- f(b, intersect(UC, colon), intersect(UC, colon))
CD_colon_beta <- f(b, intersect(CD, colon), intersect(CD, colon))

ungroup_list <- function(x) {
  data.frame(group = rep(names(x), lengths(x)),
             beta = unlist(x, FALSE, FALSE), stringsAsFactors = FALSE)
}
}
{
pdf(paste0("Figures/beta_diversity_", method, ".pdf"))
l1 <- list("C" = controls_beta, UC = UC_beta, CD = CD_beta,
           "C vs UC" = f(b, controls, UC),
           "C vs CD" = f(b, controls, CD),
           "UC vs CD" = f(b, UC, CD))
# boxplot(x = l1,
#         main = "Beta diversity", ylab = "Jaccard similarity",
#         ylim = c(0, 1), frame.plot = FALSE)
df1 <- ungroup_list(l1)
comparisons <- list(c("C", "UC"), c("C", "CD"), c("UC", "CD"),
                    c("C vs CD", "C vs UC"))
p1 <- ggplot(df1,
       aes(x = fct_relevel(group, names(l1)), y = beta)) +
  geom_violin(fill = NA) +
  geom_boxplot(width = 0.1, fill = NA) +
  ggpubr::stat_compare_means(comparisons = comparisons) +
  labs(x = element_blank(), y = paste0(method, " dissimilarity"),
       title = "Comparing all samples") +
  theme_minimal()
print(p1)
# Significant differences between UC and CD and between UC and CD.
l2 <- list(
  "C" = controls_colon_beta,
  "UC" = UC_colon_beta,
  "CD" = CD_colon_beta,
  "C vs CD" = f(b, intersect(controls, colon),
                intersect(CD, colon)),
  "C vs UC" = f(b, intersect(controls, colon),
                intersect(UC, colon)),
  "UC vs CD" = f(b, intersect(UC, colon),
                 intersect(CD, colon)))
# On colon nothing is significant
# boxplot(x = l2,
#   main = "Beta diversity on colon", ylab = "Jaccard similarity",
#   ylim = c(0, 1),frame.plot = FALSE)
df2 <- ungroup_list(l2)
comparisons <- list(c("C", "UC"), c("C", "CD"), c("UC", "CD"),
                    c("C vs CD", "C vs UC"))
p2 <- ggplot(df2,
       aes(x = fct_relevel(group, names(l2)), y = beta)) +
  geom_violin(fill = NA) +
  geom_boxplot(width = 0.1, fill = NA) +
  ggpubr::stat_compare_means(comparisons = comparisons) +
  labs(x = element_blank(), y = paste0(method, " dissimilarity"),
       title = "Comparing colon samples") +
  theme_minimal()
print(p2)
l3 <- list(
  "C" = controls_ileum_beta,
  "CD" = CD_ileum_beta,
  "C vs CD" = f(b, intersect(controls, ileum),
                intersect(CD, ileum)))
# boxplot(x = l3,
#   main = "Beta diversity on ileum", ylab = "Jaccard similarity",
#   # ylim = c(0, 1),
#   frame.plot = FALSE)
# Not significant differences between CD and C on the ileum alone.
df3 <- ungroup_list(l3)
comparisons <- list(c("C", "CD"))
p3 <- ggplot(df3,
       aes(x = fct_relevel(group, names(l3)), y = beta)) +
  geom_violin(fill = NA) +
  geom_boxplot(width = 0.1, fill = NA) +
  ggpubr::stat_compare_means(comparisons = comparisons) +
  labs(x = element_blank(), y = paste0(method, " dissimilarity"),
       title = "Comparing ileum samples") +
  theme_minimal()
print(p3)
dev.off()
}
f2()

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
