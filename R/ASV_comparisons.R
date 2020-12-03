library("phyloseq")
library("ggplot2")
library("tidyverse")
ASV <- readRDS("data/ASV.RDS")


counts_ASV <- ASV
colnames(counts_ASV) <- NULL

ASV_counts <- t(counts_ASV)

# Load and prepare metadata for the samples:
meta <- readRDS("data_out/refined_meta.RDS")
# The missing values of Exact location
meta$ileum <- ifelse(meta$Exact_location == "ileum", "Ileum", "Colon")
meta$ileum[is.na(meta$ileum )] <- "Colon"
meta$SEX[is.na(meta$SEX) | meta$SEX == ""] <- "female"
meta$Time[is.na(meta$Time)] <- "C"
meta$IBD <- as.character(meta$IBD)
meta$IBD[is.na(meta$IBD)] <- "CONTROL"

# Remove unnecessary samples to fit just the relevant
colnames(ASV_counts) <- gsub(pattern = "_S.*_L001_R.*", replacement = "",
                             x = colnames(ASV_counts))
stopifnot(any(meta$Name %in% colnames(ASV_counts))) # Check there aren't any missing samples
ASV_counts <- ASV_counts[ , colnames(ASV_counts) %in% meta$Name]
ASV_counts <- ASV_counts[, match(colnames(ASV_counts),  meta$Name)]
saveRDS(ASV_counts, file = "data_out/refined_ASV.RDS")
colnames(ASV_counts) <- meta$Original

# Diversity ####
rownames(ASV_counts) <- paste0("sp", seq_len(nrow(ASV_counts)))
rownames(meta) <- meta$Original
phyloseq <- phyloseq(otu_table(ASV_counts, taxa_are_rows = TRUE),
                     sample_data(meta))

alpha_meas <- c("Simpson", "Shannon")
theme_set(theme_minimal())
richness <- estimate_richness(phyloseq)
richness <- cbind(richness, meta)
richness$Time <- as.factor(richness$Time)
richness$Time <- forcats::fct_relevel(richness$Time, "C", after = 0)
richness$IBD <- as.factor(richness$IBD)
richness$IBD <- forcats::fct_relevel(richness$IBD, "CONTROL", after = 0)
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

date <- "20201105"
r2 <- pivot_longer(richness, cols = Chao1:Fisher, names_to = "Alpha diversity")
richness_rel <- filter(r2, `Alpha diversity` %in% c("Shannon", "Simpson")) %>%
  mutate(effective = case_when(
    `Alpha diversity` == "Shannon" ~ exp(value),
    `Alpha diversity` == "Simpson" ~ 1/value,
  ),
  resp_w14_46 = case_when(
    Time %in% c("14", "46") ~ "14/46",
    TRUE ~ as.character(Time)),
  resp_w14_46 = fct_relevel(resp_w14_46, c("C", "0", "14/46"))
  )
richness_rel %>%
  ggplot(aes(x = resp_w14_46, y = effective, col = ANTITNF_responder)) +
  geom_boxplot(alpha = 0, outlier.size = 0) +
  geom_point(position = position_jitterdodge(jitter.height = 0)) +
  facet_grid(`Alpha diversity`~ ileum, scale = "free") +
  labs(y = "Alpha diversity", x = element_blank())
ggsave(filename = paste0("Figures/", date, "_ASV_alpha_w14_46_location.png"))

richness_rel %>%
  ggplot(aes(x = resp_w14_46, y = effective, col = ANTITNF_responder)) +
  geom_boxplot(alpha = 0, outlier.size = 0) +
  geom_point(position = position_jitterdodge(jitter.height = 0)) +
  facet_wrap(~ `Alpha diversity`, scale = "free_y") +
  labs(y = "Alpha diversity", x = element_blank())
ggsave(filename = paste0("Figures/", date, "_ASV_alpha_w14_46.png"))

ggplot(richness_rel,aes(SEX, effective, col = IBD, shape = IBD)) +
  geom_boxplot(alpha = 0, outlier.size = 0) +
  geom_point(position = position_jitterdodge(jitter.height = 0)) +
  facet_grid(`Alpha diversity`~ ileum, scale = "free") +
  labs(y = "Alpha diversity", x = element_blank())
ggsave(paste0("Figures/", date, "_ASV_alpha_sex_location.png"))

ggplot(richness_rel, aes(Time, effective, col = IBD, shape = IBD)) +
  geom_boxplot(alpha = 0, outlier.size = 0) +
  geom_point(position = position_jitterdodge( jitter.height = 0)) +
  facet_grid(`Alpha diversity`~ ileum, scale = "free") +
  labs(y = "Alpha diversity", x = element_blank())
ggsave(paste0("Figures/", date, "_ASV_alpha_time_location.png"))

ggplot(richness_rel, aes(IBD, effective, col = Time, shape = Time)) +
  geom_boxplot(alpha = 0, outlier.size = 0) +
  geom_point(position = position_jitterdodge( jitter.height = 0)) +
  facet_grid(`Alpha diversity` ~ ileum, scale = "free", drop = TRUE) +
  labs(y = "Alpha diversity", x = element_blank())
ggsave(paste0("Figures/", date, "_ASV_alpha_location_ibd.png"))

ggplot(richness_rel, aes(ANTITNF_responder, effective, col = IBD, shape = IBD)) +
  geom_boxplot(alpha = 0, outlier.size = 0) +
  geom_point(position = position_jitterdodge( jitter.height = 0)) +
  facet_wrap(~`Alpha diversity`, scales = "free") +
  labs(y = "Alpha diversity", x = element_blank(), title = "Responders")
ggsave(paste0("Figures/", date, "_ASV_alpha_responders_ibd.png"))

ggplot(richness_rel, aes(ANTITNF_responder, effective, col = IBD, shape = IBD)) +
  geom_boxplot(alpha = 0, outlier.size = 0) +
  geom_point(position = position_jitterdodge( jitter.height = 0)) +
  facet_grid(`Alpha diversity` ~ ileum, scales = "free") +
  labs(y = "Alpha diversity", x = element_blank(), title = "Responders")
ggsave(paste0("Figures/", date, "_ASV_alpha_responders_ibd2.png"))
richness_rel %>%
  filter(`Alpha diversity` == "Shannon") %>%
  ggplot(aes(IBD, effective, col = ileum, shape = ileum)) +
  geom_boxplot(alpha = 0, outlier.size = 0) +
  geom_point(position = position_jitterdodge( jitter.height = 0)) +
  # facet_grid(`Alpha diversity` ~ ileum, scales = "free") +
  labs(y = "Shannon Effective", x = element_blank(), title = "Alpha diversity")
ggsave(paste0("Figures/", date, "_ASV_alpha_ibd_location.png"))

# Comparisons #####
