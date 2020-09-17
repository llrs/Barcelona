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
# Read ####
tab <- read.delim("data/Partek_Michigan3_Kraken_Classified_family.tsv", check.names = FALSE)
colnames(tab) <- gsub("_S.*", "", colnames(tab)) # Remove trailing numbers
counts <- tab[, -1]
microorganism <- tab[, 1, FALSE]

# From the QC step
meta <- readRDS("data_out/refined_meta.RDS")
otus <- bcn[, colnames(bcn) %in% meta$Name]
colnames(otus) <- meta$Original[match(colnames(otus), meta$Name)]

# Reorder samples to match!
otus <- otus[, match(meta$Original, colnames(otus))]
rownames(meta) <- meta$Original
stopifnot(all(colnames(otus) == meta$Original))
microorganism <- read.csv("data/family.csv", row.names = 1)

meta$ileum <- ifelse(meta$Exact_location == "ileum", "Ileum", "Colon")
# The missing values of Exact location
meta$ileum[is.na(meta$ileum )] <- "Colon"
meta$SEX[is.na(meta$SEX) | meta$SEX == ""] <- "female"
meta$Time[is.na(meta$Time)] <- "C"

otus <- as.matrix(otus)
rownames(otus) <- seq_len(nrow(otus))
phyloseq <- phyloseq(otu_table(otus, taxa_are_rows = TRUE),
              sample_data(meta),
              tax_table(as.matrix(microorganism)))
}
{
# Alpha diversity ####
alpha <- prop.table(otus, 2)*100
a <- as.data.frame(alpha)
a$otus <- genus[, 1]
a <- pivot_longer(a, colnames(alpha))

b <- a %>%
  group_by(name) %>%
  mutate(cumsum = cumsum(value)) %>%
  ungroup() %>%
  arrange(name) %>%
  filter(value != 0) %>%
  mutate(otus = forcats::fct_lump(otus, n = 7, w = value))
b <- merge(b, meta, by.x = "name", by.y = "Original") %>%
  arrange(IBD, ileum) %>%
  mutate(name = factor(name, levels = unique(name)))

ggplot(b) +
  geom_col(aes(name, value, fill = otus, col = otus)) +
  guides(col = FALSE) +
  theme_minimal() +
  labs(x = element_blank(), y = "%", title = "Abundance on samples",
       fill = "Families") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank()) +
  scale_fill_brewer(type= "qual") +
  scale_color_brewer(type = "qual")

# Decided to do Shannon on the 10/09/2020
# Should be the Shannon and the Simpson effective but I couldn't find how to calculate them
date <- "20200917"
alpha_meas <- c("Simpson", "Shannon")
theme_set(theme_minimal())
richness <- estimate_richness(phyloseq)
richness <- cbind(richness, meta)
richness$Time[is.na(richness$Time)] <- "C"
richness$Time <- as.factor(richness$Time)
richness$Time <- fct_relevel(richness$Time, "C", after = 0)
richness$IBD <- as.factor(richness$IBD)
richness$IBD <- fct_relevel(richness$IBD, "CONTROL", after = 0)
richness <- richness %>%
  arrange(IBD, ileum, SEX) %>%
  mutate(Original = factor(Original, levels = unique(Original)))
ggplot(richness) +
  geom_col(aes(Shannon, Original, fill = IBD, col = IBD)) +
  labs(y = "Sample") +
  theme(panel.grid.major.y = element_blank())
ggplot(richness) +
  geom_jitter(aes(x = SEX, y = Shannon, col = IBD, shape = IBD),
    height = 0, width = 0.3) +
  labs(x = element_blank())
ggplot(richness) +
  geom_jitter(aes(x = SEX, y = Simpson, col = IBD, shape = IBD),
    height = 0, width = 0.3) +
  labs(x = element_blank())

# p <- plot_richness(phyloseq, "SEX", "IBD", measures = alpha_meas)
# remove_geom(p, 'point', 1) + geom_jitter()

r2 <- pivot_longer(richness, cols = Chao1:Fisher, names_to = "Alpha diversity")
richness_rel <- filter(r2, `Alpha diversity` %in% c("Shannon", "Simpson")) %>%
  mutate(effective = case_when(
    `Alpha diversity` == "Shannon" ~ exp(value),
    `Alpha diversity` == "Simpson" ~ 1/value,
  ))

ggplot(richness_rel,aes(SEX, effective, col = IBD, shape = IBD)) +
  geom_boxplot(alpha = 0, outlier.size = 0) +
  geom_point(position = position_jitterdodge(jitter.height = 0)) +
  facet_grid(`Alpha diversity`~ ileum, scale = "free") +
  labs(y = "Alpha diversity", x = element_blank())
ggsave(paste0("Figures/", date, "_alpha_sex_location.png"))

ggplot(richness_rel, aes(Time, effective, col = IBD, shape = IBD)) +
  geom_boxplot(alpha = 0, outlier.size = 0) +
  geom_point(position = position_jitterdodge( jitter.height = 0)) +
  facet_grid(`Alpha diversity`~ ileum, scale = "free") +
  labs(y = "Alpha diversity", x = element_blank())
ggsave(paste0("Figures/", date, "_alpha_time_location.png"))

ggplot(richness_rel, aes(IBD, effective, col = Time, shape = Time)) +
  geom_boxplot(alpha = 0, outlier.size = 0) +
  geom_point(position = position_jitterdodge( jitter.height = 0)) +
  facet_grid(`Alpha diversity` ~ ileum, scale = "free", drop = TRUE) +
  labs(y = "Alpha diversity", x = element_blank())
ggsave(paste0("Figures/", date, "_alpha_location_ibd.png"))

ggplot(richness_rel, aes(ANTITNF_responder, effective, col = IBD, shape = IBD)) +
  geom_boxplot(alpha = 0, outlier.size = 0) +
  geom_point(position = position_jitterdodge( jitter.height = 0)) +
  facet_wrap(~`Alpha diversity`, scales = "free") +
  labs(y = "Alpha diversity", x = element_blank(), title = "Responders")
ggsave(paste0("Figures/", date, "_alpha_responders_ibd.png"))

ggplot(richness_rel, aes(ANTITNF_responder, effective, col = IBD, shape = IBD)) +
  geom_boxplot(alpha = 0, outlier.size = 0) +
  geom_point(position = position_jitterdodge( jitter.height = 0)) +
  facet_grid(`Alpha diversity` ~ ileum, scales = "free") +
  labs(y = "Alpha diversity", x = element_blank(), title = "Responders")
ggsave(paste0("Figures/", date, "_alpha_responders_ibd2.png"))
}
# Beta diversity ####
# Remove empty lines
e <- apply(otus, 1, function(x){sum(x == 0)})
{
method <- "bray"
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
