library("phyloseq")
library("ggplot2")
library("dplyr")
library("tidyr")
library("vegan")
# Species ####
## cut 3000 ####
pcut3000 <- read.delim("data/20210526_Partek_BCN_uncut_Kraken_Classified_species.txt",
                       check.names = FALSE, row.names = 1,
                       stringsAsFactors = FALSE)
pcut3000 <- as.matrix(pcut3000)
clean_names <- function(x) {
  x <- gsub("^[0-9]+\\.", "", x)
  gsub("_S[0-9]+$", "", x)
}


colnames(pcut3000) <- clean_names(colnames(pcut3000))
fs <- strcapture("(([0-9]{3}-w[0-9]{3}|C[0-9]-T-DM-.+)_?(p[0-9])?)", colnames(pcut3000),
                 proto = data.frame(file = character(),
                                    sample = character(),
                                    replicate = character(), stringsAsFactors = FALSE))
fs$replicate[!nzchar(fs$replicate)] <- NA
fs$OTUs <- colSums(pcut3000)
fs2 <- group_by(fs, sample) %>%
  filter(OTUs == max(OTUs))
otus <- pcut3000

meta <- readRDS("output/refined_meta_all.RDS")
meta <- meta[match(colnames(otus), meta$Original), ]
meta$ileum <- ifelse(meta$Exact_location == "ileum", "Ileum", "Colon")
meta$ileum[is.na(meta$ileum )] <- "Colon"
meta$SEX[is.na(meta$SEX) | meta$SEX == ""] <- "female"
meta$Time[is.na(meta$Time)] <- "C"
meta$IBD <- as.character(meta$IBD)
meta$IBD[is.na(meta$IBD)] <- "CONTROL"
meta$Activity <- ifelse(is.na(meta$Activity), "INACTIVE", meta$Activity)

meta <- mutate(meta,
               active = case_when(
                 IBD == "CONTROL" ~ "INACTIVE",
                 IBD == "CD" & `CDEIS_partial` <= 4 ~ "INACTIVE",
                 IBD == "UC" & `partial MAYO UC` <= 1 ~ "INACTIVE",
                 TRUE ~ "ACTIVE"))

rownames(meta) <- colnames(otus)
phyloseq <- phyloseq(otu_table(otus, taxa_are_rows = TRUE),
                     sample_data(meta)
                     # tax_table(as.matrix(family))
)


pdf("Figures/rarefaction_kraken_uncut.pdf")
v2 <- vegan::rarecurve(t(otus), sample = min(colSums(otus)), step = 200, label = FALSE,
                       main = "Rarefaction kraken min")
v2 <- vegan::rarecurve(t(otus), sample = median(colSums(otus)), step = 200, label = FALSE,
                       main = "Rarefaction kraken median")
dev.off()

theme_set(theme_bw())
beta <- estimate_richness(phyloseq)
res <- cbind(beta, meta)
res <- mutate(res, IBD = forcats::fct_relevel(IBD, c("CONTROL", "CD", "UC")))
res <- pivot_longer(res, cols = Chao1:Fisher, names_to = "Alpha diversity")
richness_rel <- filter(res, `Alpha diversity` %in% c("Shannon", "Simpson")) %>%
  mutate(effective = case_when(
    `Alpha diversity` == "Shannon" ~ exp(value),
    `Alpha diversity` == "Simpson" ~ 1/value,
  ))
ggplot(richness_rel, aes(ileum, effective, col = ANTITNF_responder)) +
  geom_boxplot(alpha = 0, outlier.size = 0) +
  geom_point(position = position_jitterdodge(jitter.height = 0, jitter.width = 1/4)) +
  facet_wrap(~ `Alpha diversity` , scales = "free_y", drop = TRUE) +
  labs(y = "Alpha diversity", x = element_blank(), title = "Species diversity",
       subtitle = "all samples") +
  theme_minimal()
ggsave("Figures/diversity_uncut_species_kraken.png")

richness_rel %>%
  filter(`Alpha diversity` == "Shannon") %>%
  ggplot(aes(Activity, effective, col = IBD)) +
  geom_boxplot(alpha = 0, outlier.size = 0) +
  geom_point(position = position_jitterdodge(jitter.height = 0, jitter.width = 1/4)) +
  facet_wrap( ~ ileum, drop = TRUE) +
  labs(y = "Shannon alpha diversity", x = element_blank(), title = "Diversity kraken at species level",
       subtitle = "all samples") +
  theme_minimal() +
  scale_y_continuous(limits = c(0, NA))
ggsave("Figures/uncut_kraken_species_diversity.png")


seqs <- read.delim("output/reads.txt", sep = " ",header = FALSE)
fs <- strcapture("(([0-9]{3}-w[0-9]{3}|C[0-9]-T-DM-.+)_?(p[0-9])?_.*)", seqs$V1,
                 proto = data.frame(file = character(),
                                    sample = character(),
                                    replicate = character(), stringsAsFactors = FALSE))
fs$replicate[!nzchar(fs$replicate)] <- NA
fs2 <- merge(seqs, fs, by.x = "V1", by.y = "file")
r <- richness_rel[richness_rel$`Alpha diversity` == "Shannon", ]
y <- merge(fs2[c(TRUE, FALSE), ], r, by.x = "sample", by.y = "Original")

pdf("Figures/diversity_and_reads.pdf")
r1 <- cor.test(y$V2, y$effective)
cor_plots_labels <- function(x) {
  paste("cor =", round(x$estimate, 3), "p-value =", round(x$p.value, 4))
}
plot(y$V2, y$effective, xlab = "Reads", ylab = "Shannon alpha diversity",
     main = paste0("Diversity and reads\n", cor_plots_labels(r1)))
abline(v = 3000, col = "red")
abline(v = 2*sd(seqs$V2)+mean(seqs$V2), col = "red")
r2 <- cor.test(log10(y$V2), y$effective)
plot(log10(y$V2), y$effective, xlab = "Reads (log10)", ylab = "Shannon alpha diversity",
     main = paste0("Diversity and reads\n", cor_plots_labels(r2)))
abline(v = log10(3000), col = "red")
abline(v = log10(2*sd(seqs$V2)+mean(seqs$V2)), col = "red")
dev.off()

