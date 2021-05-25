library("phyloseq")
library("ggplot2")
library("dplyr")
library("tidyr")
# Species ####
## cut 3000 ####
pcut3000 <- read.delim("data/20210520_Partek_BCN_cut3000_Kraken_Classified_species.txt",
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
otus <- pcut3000[, fs2$file]
colnames(otus) <- fs2$sample
saveRDS(fs2, "output/samples_cut3000.RDS")

meta <- readRDS("output/refined_meta_all.RDS")
meta <- meta[match(colnames(otus), meta$Original), ]
rownames(meta) <- colnames(otus)
meta$Loc <- ifelse(meta$Exact_location == "ileum", "ileum", "colon")
meta$Activity <- ifelse(is.na(meta$Activity), "INACTIVE", meta$Activity)
phyloseq <- phyloseq(otu_table(otus, taxa_are_rows = TRUE),
                     sample_data(meta)
                     # tax_table(as.matrix(family))
)
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
ggplot(richness_rel, aes(Loc, effective, col = ANTITNF_responder)) +
  geom_boxplot(alpha = 0, outlier.size = 0) +
  geom_point(position = position_jitterdodge(jitter.height = 0, jitter.width = 1/4)) +
  facet_wrap(~ `Alpha diversity` , scales = "free_y", drop = TRUE) +
  labs(y = "Alpha diversity", x = element_blank(), title = "Phylum diversity",
       subtitle = "cutoff 3000 reads") +
  theme_minimal()
ggsave("Figures/diversity_phylum_cut_3000_kraken.png")

richness_rel %>%
  filter(`Alpha diversity` == "Shannon") %>%
  ggplot(aes(Activity, effective, col = IBD)) +
  geom_boxplot(alpha = 0, outlier.size = 0) +
  geom_point(position = position_jitterdodge(jitter.height = 0, jitter.width = 1/4)) +
  facet_wrap( ~ Loc, drop = TRUE) +
  labs(y = "Shannon alpha diversity", x = element_blank(), title = "Diversity kraken",
       subtitle = "cutoff 3000 reads") +
  theme_minimal() +
  scale_y_continuous(limits = c(0, NA))
ggsave("Figures/cut3000_kraken_phylum_diversity.png")
