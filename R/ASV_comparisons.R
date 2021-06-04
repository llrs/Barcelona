library("phyloseq")
library("ggplot2")
library("tidyverse")
ASV <- readRDS("data/ASV.RDS")

pdf("Figures/rarefaction_ASV_uncut.pdf")
v2 <- vegan::rarecurve(ASV, sample = min(rowSums(ASV)), step = 200, label = FALSE,
                       main = "Rarefaction ASV")
v2 <- vegan::rarecurve(ASV, sample = median(rowSums(ASV)), step = 200, label = FALSE,
                       main = "Rarefaction ASV median")
dev.off()
rs <- vegan::rareslope(ASV, median(rowSums(ASV)))
rs2 <- vegan::rareslope(ASV, min(rowSums(ASV)))


counts_ASV <- ASV
colnames(counts_ASV) <- NULL

ASV_counts <- t(counts_ASV)

# Load and prepare metadata for the samples:
meta <- readRDS("output/refined_meta.RDS")
# The missing values of Exact location
meta$ileum <- ifelse(meta$Exact_location == "ileum", "Ileum", "Colon")
meta$ileum[is.na(meta$ileum )] <- "Colon"
meta$SEX[is.na(meta$SEX) | meta$SEX == ""] <- "female"
meta$Time[is.na(meta$Time)] <- "C"
meta$IBD <- as.character(meta$IBD)
meta$IBD[is.na(meta$IBD)] <- "CONTROL"
meta <- meta[!startsWith(meta$Original, "052"), ] # Remove outliers



# Remove unnecessary samples to fit just the relevant
colnames(ASV_counts) <- gsub("_S.*", "", colnames(ASV_counts))
colnames(ASV_counts) <- gsub("_p.*", "", colnames(ASV_counts))
stopifnot(any(meta$Name %in% colnames(ASV_counts))) # Check there aren't any missing samples

ASV_counts <- ASV_counts[ , colnames(ASV_counts) %in% meta$Original]
ASV_counts <- ASV_counts[, match(meta$Original, colnames(ASV_counts))]
saveRDS(ASV_counts, file = "output/refined_ASV.RDS")
colnames(ASV_counts) <- meta$Original

# Filter low abundance taxa
# abundance_ASV <- apply(ASV_counts, 2, function(x){x/sum(x) < 0.5/100})
# keep_ASV <- rowSums(abundance_ASV) < 120
# Filter low mapped samples
ASV_counts_filtered <- ASV_counts[, colSums(ASV_counts) > 3000]
ASV_counts_filtered <- ASV_counts_filtered[rowSums(ASV_counts_filtered) != 0, ]

# Diversity ####
ASV_counts_filtered <- ASV_counts[, !startsWith(colnames(ASV_counts), "C")]
rownames(ASV_counts_filtered) <- paste0("sp", seq_len(nrow(ASV_counts_filtered)))

rownames(meta) <- meta$Original
meta <- mutate(meta,
             active = case_when(
               IBD == "CONTROL" ~ "INACTIVE",
               IBD == "CD" & `CDEIS_partial` <= 4 ~ "INACTIVE",
               IBD == "UC" & `partial MAYO UC` <= 1 ~ "INACTIVE",
               TRUE ~ "ACTIVE"))

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

date <- "20210612"
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
  ggplot(aes(active, effective, col = IBD)) +
  geom_boxplot(alpha = 0, outlier.size = 0) +
  geom_point(position = position_jitterdodge( jitter.height = 0)) +
  facet_grid(~ ileum, scales = "free") +
  labs(y = "Shannon diversity", x = element_blank(),
       title = "Activity ASV uncut")
ggsave("Figures/20210525_uncut_ASV_diversity.png")
ggplot(richness_rel, aes(active, effective, col = IBD)) +
  geom_boxplot(alpha = 0, outlier.size = 0) +
  geom_point(position = position_jitterdodge( jitter.height = 0)) +
  facet_grid(`Alpha diversity` ~ ileum, scales = "free") +
  labs(y = "Alpha diversity", x = element_blank(), title = "Activity")
ggsave(paste0("Figures/", date, "_ASV_alpha_active.png"))

ggplot(richness_rel, aes(Time, col = IBD, effective, fill = IBD, shape = active)) +
  geom_boxplot(alpha = 0, outlier.size = 0) +
  geom_point(position = position_jitterdodge( jitter.height = 0)) +
  facet_grid(`Alpha diversity` ~ ileum, scales = "free") +
  labs(y = "Alpha diversity", x = element_blank(), title = "Disease")
ggsave(paste0("Figures/", date, "_ASV_alpha_time_active.png"))


richness_rel %>%
  filter(`Alpha diversity` == "Shannon") %>%
  ggplot(aes(IBD, effective, col = ileum, shape = ileum)) +
  geom_boxplot(alpha = 0, outlier.size = 0) +
  geom_point(position = position_jitterdodge( jitter.height = 0)) +
  # facet_grid(`Alpha diversity` ~ ileum, scales = "free") +
  labs(y = "Shannon Effective", x = element_blank(), title = "Alpha diversity")
ggsave(paste0("Figures/", date, "_ASV_alpha_ibd_location.png"))

# Comparisons #####


controls_trim <- readRDS("../TRIM/alpha_diversity_aTNF_controls.RDS") %>%
  mutate(Sequenced = "TRIM",
          Exact_location = tolower(Exact_location),
          Exact_location = gsub(" colon", "", Exact_location))
both_studies <- richness_rel %>%
  filter(IBD == "CONTROL") %>%
  mutate(Sequenced = "aTNF") %>%
  droplevels() %>%
  full_join(controls_trim, by = c("Original" = "Sample_Code_uDNA", "SEX" = "SEX",
                                  "ID" = "ID", "Sequenced" = "Sequenced",
                                  "effective" = "effective",
                                  "value" = "value", "Exact_location" = "Exact_location",
                                  "Alpha diversity" = "Alpha diversity"),
            suffix = c(".aTNF", ".TRIM"))

both_studies %>%
  filter(`Alpha diversity` == "Shannon") %>%
  mutate(Sequenced = forcats::fct_relevel(Sequenced, "TRIM", "aTNF")) %>%
  ggplot() +
  geom_point(aes(Sequenced, effective, col = Original)) +
  geom_line(aes(Sequenced, effective, col = Original, group = Original)) +
  labs(y = "Shannon Effective", title = "Differences in alpha diversity",
       x = "Study", col = "Sample")
ggsave("Figures/controls_at_aTNF_and_TRIM_diversity.png")
both_studies %>%
  filter(`Alpha diversity` == "Shannon") %>%
  mutate(Sequenced = forcats::fct_relevel(Sequenced, "TRIM", "aTNF")) %>%
  ggplot() +
  geom_point(aes(Sequenced, effective, col = Original)) +
  geom_line(aes(Sequenced, effective, col = Original, group = Original)) +
  facet_wrap(~Exact_location) +
  labs(y = "Shannon Effective", title = "Differences in alpha diversity",
       x = "Study", col = "Sample")
