# load data ####
library("readxl")
library("tidyr")
library("metagenomeSeq")
library("integration")
library("stringr")
library("purrr")
library("ggplot2")
library("lubridate")
library("org.Hs.eg.db")
library("forcats")
library("ggh4x") # from teunbrand/ggh4x Now on CRAN
library("ggpubr")
library("dplyr")
library("rstatix") # To compare several groups the tidy way

# microbiome data ####
{
  counts_ASV <- readRDS("data/ASV.RDS")
  ASV <- colnames(counts_ASV)
  colnames(counts_ASV) <- NULL

  tab <- t(counts_ASV)
  # Remove trailing numbers
  colnames(tab) <- gsub("_S.*", "", colnames(tab))
  colnames(tab) <- gsub("_p.*", "", colnames(tab))

  meta <- readRDS("output/info_samples.RDS")
  depth <- colSums(tab)
  # Remove duplicate samples and keep the ones with more sequencing depth
  replicates <- table(colnames(tab))
  replicate_samples <- meta[meta$Original %in% names(replicates[replicates > 1]), ]
  depth[names(depth) %in% names(replicates)[replicates > 1]]
  w <- which(names(depth) %in% names(replicates)[replicates > 1])
  df <- data.frame( w = w,
              sample = names(depth)[w],
              depth = depth[w])
  remove_samples <- df %>%
    group_by(sample) %>%
    filter(depth != max(depth)) %>%
    pull(w)
  micros <- tab[, -remove_samples]
  # Remove the two outliers:
  micros <- micros[, !grepl("^052", colnames(micros))]


}

# Group by tax level ####
group_taxa <- function(taxonomy, data, groups) {
  microorganism <-  readRDS("output/taxonomy_ASV.RDS")$tax
  microorganism <- cbind(microorganism, "ASV" = rownames(microorganism))
  microorganism <- as.data.frame(microorganism, stringsAsFactors = FALSE)
  microorganism$rowname <- as.character(seq_len(nrow(microorganism)))
  rownames(microorganism) <- as.character(seq_len(nrow(microorganism)))
  # Remove Mitochondria
  mit <- microorganism[, "Family"] != "Mitochondria"
  microorganism <- microorganism[mit, ]

  tab2 <- data[mit, ]
  colnames(tab2) <- seq_along(colnames(data))
  tab2 <- cbind(tab2, microorganism)

  # Collapse the data so that there is one Genus per sample, not multiple
  rm_micro <- tab2 %>%
    group_by(!!!syms(groups)) %>%
    summarize(across(where(is.numeric), function(x){sum(x, na.rm = TRUE)}),
              .groups = "keep") %>%
    as.data.frame()
  # Change the first number according to the number of columns on group_by +1
  colnames(rm_micro)[(length(groups) + 1):ncol(rm_micro)] <- colnames(data)
  genus <- rm_micro[, 1:(length(groups))] # And here the same number
  list(genus = genus, counts = rm_micro[, (length(groups) + 1):ncol(rm_micro)])
}

# Prepare for GETS ####
A <- readRDS("data/RGCCA_data_wo_out.RDS")
meta_f <- A$Meta

# * Genus #####
genus <- group_taxa(taxonomy = readRDS("output/taxonomy_ASV.RDS")$tax,
                    data = micros,
                    groups = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"))

rS <- rowSums(genus$counts)
genus_prev <- apply(genus$counts[rS != 0, ], 2, function(x){x/sum(x, na.rm = TRUE)})
meta2 <- meta_f[match(colnames(genus_prev), meta_f$Original), c("Original", "IBD", "SEX", "Exact_location", "Ulcers", "sample_location", "diagTime", "Age", "AgeDiag", "Time")] %>%
  arrange(IBD, sample_location, SEX, Age)

write.table(meta2, file = "output/GETS_ASV_samples.tsv", quote = FALSE, sep = "\t", row.names = FALSE, na = " ")
genus_prev <- genus_prev[, match(meta2$Original, colnames(genus_prev))]
stopifnot("Some Nas present" = !anyNA(genus_prev))
micro_norm <- norm_RNAseq(genus$counts[rS != 0, ])
micro_norm <- micro_norm[, meta2$Original]
matrix_norm <- cbind(rownames = rownames(micro_norm), micro_norm)
write.table(matrix_norm, file = "output/GETS_matrix_ASV_norm.tsv", quote = FALSE, sep = "\t", row.names = FALSE)
system2("gzip", args = "-kf output/GETS_matrix_ASV_norm.tsv") # Compress it to upload to website

tax <- genus$genus[rS != 0, ]
matrix_out <- cbind(rownames = rownames(genus_prev), genus_prev)
write.table(matrix_out, file = "output/GETS_matrix_ASV_prev.tsv", quote = FALSE, sep = "\t", row.names = FALSE)
system2("gzip", args = "-kf output/GETS_matrix_ASV_prev.tsv") # Compress it to upload to website
write.table(tax, file = "output/GETS_ASV_genus.tsv", quote = FALSE, sep = "\t", row.names = FALSE, na = "")
system2("gzip", args = "-kf output/GETS_ASV_genus.tsv") # Compress it to upload to website

meta2$sample_location[is.na(meta2$sample_location)] <- ifelse(meta2$Exact_location[is.na(meta2$sample_location)] == "ileum", "ileum", "colon")

# ** Plots of interesting results ####
gen <- genus$counts
rs <- rowSums(gen)
gen <- apply(gen[rs != 0, ], 2, function(x){x/sum(x, na.rm = TRUE)})
n <- cbind(gen, genus$genus[rs != 0, ])
m <- pivot_longer(n, cols = where(is.numeric)) %>%
  # filter(value != 0) %>%
  rename(Original = name, percentage = value) %>%
  left_join(meta2, by = "Original") %>%
  mutate(Ulcers = case_when(IBD == "CONTROL" ~ "no", TRUE ~ Ulcers),
         IBD = fct_relevel(IBD, c("CONTROL", "UC", "CD")))
# Duplicate controls with "ulcers" to be able to make the tests easily
dup_controls <- m %>% filter(IBD == "CONTROL") %>% mutate(Ulcers = "yes")
m <- rbind(m, dup_controls)

m %>% group_by(sample_location) %>%
  summarize(across(starts_with("percentage"), .fns = list(min = min , median = median, mean = mean, max = max)))
m %>% group_by(sample_location, IBD) %>%
  summarize(across(starts_with("percentage"), .fns = list(min = min , median = median, mean = mean, max = max)))

sps <- c("[Ruminococcus] gnavus group", "Lachnospira",
         "[Eubacterium] eligens group", "Monoglobus", "UCG-005", "Sellimonas",
         "CAG-56", "Ruminococcus")
m %>%
  filter(Genus %in% sps ) %>%
  group_by(sample_location, IBD, Genus) %>%
  summarize(mean_se(percentage)) %>%
  ungroup() %>%
  ggplot() +
  geom_point(aes(IBD, y, group = Genus, col = Genus)) +
  geom_errorbar(aes(IBD, y, group = Genus, col = Genus, ymin = ymin, ymax = ymax)) +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~sample_location) + labs(y = "%") +
  theme_minimal()
ggsave("Figures/species_location.png")

m %>%
  filter(Genus %in% sps ) %>%
  group_by(sample_location, IBD, Genus, Ulcers) %>%
  summarize(mean_se(percentage)) %>%
  ungroup() %>%
  mutate(IBD = fct_relevel(IBD, c("CONTROL", "UC", "CD"))) %>%
  ggplot() +
  geom_point(aes(IBD, y, group = Genus, col = Genus)) +
  geom_errorbar(aes(IBD, y, group = Genus, col = Genus, ymin = ymin, ymax = ymax)) +
  scale_y_continuous(labels = scales::percent) +
  facet_grid(Ulcers~sample_location, drop = TRUE) +
  labs(y = "%") +
  theme_minimal()  +
  theme(strip.text.y = element_text(angle = 0))
ggsave("Figures/species_location_ulcers.png")

fms <- c("Christensenellaceae", "Ruminococcaceae")
m %>%
  filter(Family %in% fms ) %>%
  group_by(sample_location, IBD, Family) %>%
  summarize(mean_se(percentage)) %>%
  ungroup() %>%
  mutate(IBD = fct_relevel(IBD, c("CONTROL", "UC", "CD"))) %>%
  ggplot() +
  geom_point(aes(IBD, y, group = Family, col = Family)) +
  geom_errorbar(aes(IBD, y, group = Family, col = Family, ymin = ymin, ymax = ymax)) +
  facet_wrap(~sample_location, drop = TRUE) +
  scale_y_continuous(labels = scales::percent) +
  labs(y = "%") +
  theme_minimal()
ggsave("Figures/family_location.png")
m %>%
  filter(Family %in% fms ) %>%
  group_by(sample_location, IBD, Family, Ulcers) %>%
  summarize(mean_se(percentage)) %>%
  ungroup() %>%
  ggplot() +
  geom_point(aes(IBD, y, group = Family, col = Family)) +
  geom_errorbar(aes(IBD, y, group = Family, col = Family, ymin = ymin, ymax = ymax)) +
  scale_y_continuous(labels = scales::percent) +
  facet_grid(Ulcers~sample_location, drop = TRUE) +
  labs(y = "%") +
  theme_minimal()  +
  theme(strip.text.y = element_text(angle = 0))
ggsave("Figures/family_location_ulcers.png")
# ** Statistics ####
# *** Family ####
# Compare location and family vs control
m %>%
  filter(Family %in% fms ) %>%
  group_by(sample_location, Family) %>% # Can't use ulcers here
  t_test(percentage~IBD, ref.group = "CONTROL") %>%
  select(Family, sample_location, group1, group2, p, p.adj) %>%
  xlsx::write.xlsx(file = "output/statistics_family_sample_location_controls.xlsx")
# Compare no ulcers, on colon and ileum vs controls
m %>%
  filter(Family %in% fms ) %>%
  group_by(sample_location, Family, Ulcers) %>% # Can't use ulcers here
  t_test(percentage~IBD, ref.group = "CONTROL") %>%
  select(Family, Ulcers, sample_location, group1, group2, p, p.adj) %>%
  xlsx::write.xlsx(file = "output/statistics_family_sample_location_ulcers.xlsx")
# Compare between families
m %>%
  filter(Family %in% fms ) %>%
  group_by(sample_location, Ulcers, IBD) %>% # Can't use ulcers here
  t_test(percentage~Family, ref.group = "Christensenellaceae") %>%
  mutate(p.adj = p.adjust(p)) %>%
  select(IBD, Ulcers, sample_location, group1, group2, p, p.adj) %>%
  xlsx::write.xlsx(file = "output/statistics_family_sample_location_ulcers_IBD.xlsx")

m %>%
  filter(Family %in% fms ) %>%
  group_by(sample_location, Family, Ulcers) %>%
  t_test(percentage~IBD, ref.group = "all")
# *** Genus ####

m %>%
  filter(Genus %in% sps ) %>%
  group_by(sample_location, Genus, Ulcers) %>%
  t_test(percentage~IBD, ref.group = "CONTROL") %>%
  select(Genus, Ulcers, sample_location, group1, group2, p, p.adj) %>%
  xlsx::write.xlsx(file = "output/statistics_genus_sample_location_ulcers.xlsx")

m %>%
  filter(Genus %in% sps ) %>%
  group_by(sample_location, Genus) %>% # Can't use ulcers here
  t_test(percentage~IBD, ref.group = "CONTROL") %>%
  select(Genus, sample_location, group1, group2, p, p.adj) %>%
  xlsx::write.xlsx(file = "output/statistics_genus_sample_location_controls.xlsx")

# * Classes #####
class <- group_taxa(taxonomy = readRDS("output/taxonomy_ASV.RDS")$tax,
                    data = micros,
                    groups = c("Kingdom", "Phylum", "Class"))

rS <- rowSums(class$counts)
genus_prev <- apply(class$counts[rS != 0, ], 2, function(x){x/sum(x, na.rm = TRUE)})
genus_prev <- genus_prev[, match(meta2$Original, colnames(genus_prev))]
stopifnot("Some Nas present" = !anyNA(genus_prev))
tax <- class$genus[rS != 0, ]
matrix_out <- cbind(rownames = rownames(genus_prev), genus_prev)
write.table(matrix_out, file = "output/GETS_matrix_ASV_prev_fam.tsv", quote = FALSE, sep = "\t", row.names = FALSE)
system2("gzip", args = "-kf output/GETS_matrix_ASV_prev_fam.tsv") # Compress it to upload to website
write.table(tax, file = "output/GETS_ASV_fam.tsv", quote = FALSE, sep = "\t", row.names = FALSE, na = "")
system2("gzip", args = "-kf output/GETS_ASV_fam.tsv") # Compress it to upload to website
