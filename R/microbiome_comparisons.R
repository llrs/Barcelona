# load data ####
library("readxl")
library("dplyr")
library("tidyr")
library("metagenomeSeq")
library("integration")
library("stringr")
library("purrr")
library("ggplot2")
library("lubridate")
library("org.Hs.eg.db")
library("forcats")

tab <- read.delim("data/Partek_Michigan3_Kraken_Classified_genus.tsv", check.names = FALSE)
colnames(tab) <- gsub("_S.*", "", colnames(tab)) # Remove trailing numbers
counts <- tab[, -1]
genus <- tab[, 1, FALSE]

# From the QC step
meta <- readRDS("data_out/refined_meta.RDS")
ccounts <- colSums(counts)
meta <- meta[meta$Name %in% colnames(counts), ]

# Subset and change the names
otus <- counts[, colnames(counts) %in% meta$Name]
colnames(otus) <- meta$Original[match(colnames(otus), meta$Name)]

# Reorder samples to match!
OTUs2 <- otus[, match(meta$Original, colnames(otus))]
stopifnot(all(colnames(OTUs2) == meta$Original))
OTUs2 <- as.matrix(OTUs2)

rownames(meta) <- colnames(OTUs2)

meta$ileum <- ifelse(meta$Exact_location == "ileum", "Ileum", "colon")
# The missing values of Exact location
meta$ileum[is.na(meta$ileum )] <- "colon"
meta$SEX[is.na(meta$SEX) | meta$SEX == ""] <- "female"
meta$Time[is.na(meta$Time)] <- "C"


# Abundance ####
# * Family level ####
family <- read.delim("data/Partek_Michigan3_Kraken_Classified_family.tsv",
                     check.names = FALSE)
family_tax <- family[, 1]
fam <- family[ , -1]
colnames(fam) <- gsub("_S.*", "", colnames(fam)) # Remove trailing numbers
fam2 <- fam[, match(meta$Original, colnames(fam))]
fam2 <- as.matrix(fam2)

tidy_family <- family %>%
  gather(Sample, Count, -'Sample name') %>%
  mutate(Sample = str_remove(string = Sample, pattern = "_S.*"),
         # Sample = str_replace_all(Sample, "-T-DM-", "-T-TTR-"),
         ) %>%
  filter(Sample %in% meta$Name) %>%
  dplyr::rename(Family = "Sample name") %>%
  filter(!str_detect(Sample, "^500") &
           str_detect(Sample, "^C|^([0-9]+-w)")) %>%
  group_by(Sample) %>%
  mutate(ratio = Count/sum(Count)) %>%
  ungroup() %>%
  group_by(Family) %>%
  ungroup() %>%
  left_join(meta, by = c("Sample" = "Name")) %>%
  mutate_if(is.factor, as.character) %>%
  # mutate(IBD = case_when(grepl(x = Name, pattern = "^C") & is.na(IBD) ~ "CONTROL",
  #                        TRUE ~ IBD)) %>%
  arrange(IBD, Time, SEX) %>%
  mutate(Sample = Original)

tidy_family %>%
  ggplot() +
  geom_point(aes(Sample, Family, size = ratio, col = IBD)) +
  theme(axis.text.x = element_blank())
tidy_family %>%
  arrange(IBD) %>%
  ggplot() +
  geom_col(aes(Sample, ratio, col = IBD, fill = Family)) +
  theme(axis.text.x = element_blank()) +
  labs(fill = element_blank()) +
  guides(fill = FALSE) +
  scale_y_continuous(labels = scales::percent, expand = expansion())
tidy_family <- tidy_family %>%
  mutate(IBD = fct_relevel(IBD, c("CONTROL", "UC", "CD")),
         Time = fct_relevel(Time, c("C", "0", "14", "46")))

tidy_family %>%
  filter(Family == "Enterobacteriaceae") %>%
  ggplot() +
  geom_jitter(aes(IBD, ratio, shape = IBD)) +
  labs(x = element_blank(), y = "Beta diversity") +
  theme_minimal()
tidy_family %>%
  filter(Family == "Streptococcaceae") %>%
  ggplot() +
  geom_jitter(aes(IBD,ratio, shape = IBD)) +
  labs(x = element_blank(), y = "Beta diversity") +
  theme_minimal()
library("ggh4x") # from teunbrand/ggh4x
tidy_family %>%
  filter(Family == "Streptococcaceae") %>%
  ggplot() +
  geom_jitter(aes(IBD, ratio, col = Time, shape = IBD)) +
  labs(x = element_blank(), y = "Beta diversity", title = "Streptococcaceae") +
  theme_minimal() +
  facet_nested(~Time + IBD,
               scales = "free_x", switch = "x", nest_line = TRUE) +
  theme(axis.text.x = element_blank())
tidy_family %>%
  filter(Family == "Streptococcaceae") %>%
  ggplot() +
  geom_jitter(aes(IBD, ratio, shape = IBD)) +
  labs(x = element_blank(), y = "Beta diversity", title = "Streptococcaceae") +
  theme_minimal() +
  facet_nested(~ IBD + Activity,
               scales = "free_x", switch = "x", nest_line = TRUE) +
  theme(axis.text.x = element_blank())
tidy_family %>%
  filter(Family == "Enterobacteriaceae") %>%
  ggplot() +
  geom_jitter(aes(IBD, ratio, shape = IBD)) +
  labs(x = element_blank(), y = "Beta diversity", title = "Enterobacteriaceae") +
  theme_minimal() +
  facet_nested(~ IBD + Activity,
               scales = "free_x", switch = "x", nest_line = TRUE) +
  theme(axis.text.x = element_blank())
tidy_family %>%
  filter(Family == "Enterobacteriaceae") %>%
  mutate(Ileum = ifelse(Exact_location == "ileum", "ileum", "colon")) %>%
  ggplot() +
  geom_jitter(aes(Ileum, ratio, shape = Ileum)) +
  labs(x = element_blank(), y = "Beta diversity", title = "Enterobacteriaceae") +
  theme_minimal() +
  facet_nested(~ IBD + Activity + Ileum,
               scales = "free_x", switch = "x", nest_line = TRUE) +
  theme(axis.text.x = element_blank())

tidy_family %>%
  arrange(IBD) %>%
  mutate(f2 = fct_lump_prop(as.factor(Family), prop = 0.01, w = ratio),
         # Reorder such that the levels are according to the ratio on the plot
         f2 = fct_reorder2(f2, Sample, ratio)) %>%
  ggplot() +
  geom_col(aes(Sample, ratio, fill = f2)) +
  theme(axis.text.x = element_blank()) +
  labs(fill = "Family", title = "Family abundance of the samples",
       y = "Abundance") +
  scale_fill_viridis_d() +
  # guides(fill = FALSE) +
  scale_y_continuous(labels = scales::percent, expand = expansion())

# * genus level ####

tidy_genus <- tab %>%
  gather(Sample, Count, -"Sample name") %>%
  mutate(Sample = str_remove(string = Sample, pattern = "_S.*"),
         # Sample = str_replace_all(Sample, "-T-DM-", "-T-TTR-"),
  ) %>%
  filter(Sample %in% meta$Name) %>%
  dplyr::rename(Genus = "Sample name") %>%
  filter(!str_detect(Sample, "^500") &
           str_detect(Sample, "^C|^([0-9]+-w)")) %>%
  group_by(Sample) %>%
  mutate(ratio = Count/sum(Count)) %>%
  ungroup() %>%
  group_by(Genus) %>%
  ungroup() %>%
  left_join(meta, by = c("Sample" = "Name")) %>%
  mutate_if(is.factor, as.character) %>%
  # mutate(IBD = case_when(grepl(x = Name, pattern = "^C") & is.na(IBD) ~ "CONTROL",
  #                        TRUE ~ IBD)) %>%
  arrange(IBD, Time, SEX) %>%
  mutate(Sample = Original)


fd <- AnnotatedDataFrame(genus)
# Create the object
MR <- newMRexperiment(
  OTUs2,
  phenoData = AnnotatedDataFrame(meta),
  featureData = fd
)

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
# Fill some gaps
# Checked manually on the database
pData(MRtrim)$ileum[is.na(pData(MRtrim)$ileum )] <- "colon"
pData(MRtrim)$SEX[is.na(pData(MRtrim)$SEX) | pData(MRtrim)$SEX == ""] <- "female"

ibd <- as.character(pData(MRtrim)$IBD)
ibd[is.na(ibd)] <- "C"
ibd[ibd == "CONTROL"] <- "C"
ibd <- factor(ibd, levels = c("C", "CD", "UC"))
mod <- model.matrix(~ts+pData(MRtrim)$ileum+ibd)
colnames(mod) <- c("(Intercept)", "ts14", "ts46", "Ileum", "CD","UC")


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

tt <- topTable(fit2, coef = "ileum_vs_colon", number = Inf)
rownames(tt) <- genus[rownames(tt), 1]
write.csv(tt, "data_out/colon_vs_ileum.csv", row.names = TRUE)

# Prevalence ####

filter_prev <- function(x) {
  x <- x[apply(x, 1, function(y){any(y != 1)}), , drop = FALSE]
  x[, apply(x, 2, function(y){any(y != 1)}), drop = FALSE]
}

extract_genus <- function(x) {
  unique(rownames(which(x < 0.05, arr.ind = TRUE)))
}

count_prevalence <- function(tidy_data, selected_genus, ...) {
  tidy_data %>%
    mutate(presence = Count != 0,
           Time = if_else(is.na(Time), "C", Time)) %>%
    group_by(Genus, presence, ...) %>%
    count() %>%
    ungroup() %>%
    group_by(Genus, ...) %>%
    mutate(pos = sum(n), label = paste(n, collapse = "/")) %>%
    ungroup() %>%
    filter(Genus %in% selected_genus)
}

# * genus level ####
rownames(OTUs2) <- genus[, 1]

gen <- filter_prev(comb_prevalence(OTUs2, meta, c("Time")))
gen_se <- comb_prevalence(OTUs2, meta, c("Time", "SEX")) %>% filter_prev()
tre <- comb_prevalence(OTUs2, meta, c("treatment")) %>% filter_prev()
se <- comb_prevalence(OTUs2, meta, c("SEX")) %>% filter_prev()
std <- comb_prevalence(OTUs2, meta, c("Study")) %>% filter_prev()
ibdm <- comb_prevalence(OTUs2, meta, c("IBD", "Time", "ileum")) %>% filter_prev()
ibd <- comb_prevalence(OTUs2, meta, c("IBD", "ileum")) %>% filter_prev()
loc <- comb_prevalence(OTUs2, meta, c("ileum")) %>% filter_prev()
write.csv(ibd, "data_out/prevalence_disease_location.csv", row.names = TRUE)

p <- full_prevalence(OTUs2, meta, "Time")
p <- full_prevalence(OTUs2, meta, "ileum")
genus_ibd <- extract_genus(full_prevalence(OTUs2, meta, "IBD"))
genus_study <- extract_genus(full_prevalence(OTUs2, meta, "Study"))

meta2 <- meta
meta2$ti <- paste(meta$Time, meta$ileum, sep = " & ")
meta2$tI <- paste(meta$Time, meta$IBD, sep = " & ")
meta2$t3 <- paste(meta$Time, meta$Study, sep = " & ")
meta2$t4 <- paste(meta$Study, meta$IBD, sep = " & ")
meta2$t5 <- paste(meta$ileum, meta$IBD, sep = " & ")
meta2$t6 <- paste(meta$ileum, meta$Study, sep = " & ")
genus_time_ileum <- extract_genus(full_prevalence(OTUs2, meta2, "ti"))
genus_time_ibd <- extract_genus(full_prevalence(OTUs2, meta2, "tI"))
genus_time_t3 <- extract_genus(full_prevalence(OTUs2, meta2, "t3"))
genus_time_t4 <- extract_genus(full_prevalence(OTUs2, meta2, "t4"))
genus_time_t5 <- extract_genus(full_prevalence(OTUs2, meta2, "t5"))
genus_time_t6 <- extract_genus(full_prevalence(OTUs2, meta2, "t6"))



(time_genus <- extract_genus(gen))
(time_sex_genus <- extract_genus(gen_se))
(tre_genus <- extract_genus(tre))
(se_genus <- extract_genus(se))
(std_genus <- extract_genus(std))
(ibdm_genus <- extract_genus(ibdm))
(ibd_genus <- extract_genus(ibd))
(loc_genus <- extract_genus(loc))

# So basically I need to plot for ibdm and ibd_genus



theme_set(theme_minimal())

count_prevalence(tidy_genus, ibdm_genus, ileum, Time) %>%
  ggplot() +
  geom_col(aes(fct_relevel(Time, c("C", "0", "14", "46")), n, fill = presence)) +
  geom_text(aes(y = pos * 1.05, x = Time, label = label)) +
  facet_grid(ileum~Genus, scales = "free") +
  labs(x = element_blank(), y = "Samples")
ggsave("Figures/time_location_genus.png")

count_prevalence(tidy_genus, ibd_genus, IBD, ileum) %>%
  ggplot() +
  geom_col(aes(IBD, n, fill = presence)) +
  geom_text(aes(y = pos + 3, x = IBD, label = label)) +
  facet_grid(ileum~Genus) +
  labs(x = element_blank(), y = "Samples")
ggsave("Figures/coprothermobacter_prevotella_disease_location.png")

count_prevalence(tidy_genus, std_genus, Study) %>%
  ggplot() +
  geom_col(aes(Study, n, fill = presence)) +
  geom_text(aes(y = pos + 3, x = Study, label = label)) +
  labs(x = "Disease", y = "Samples", title = std_genus)
ggsave("Figures/ralstonia_study.png")

count_prevalence(tidy_genus, genus_ibd, IBD) %>%
  ggplot() +
  geom_col(aes(fct_relevel(IBD, c("CONTROL", "UC", "CD")), n, fill = presence)) +
  geom_text(aes(y = pos + 3, x = IBD, label = label)) +
  facet_grid(~Genus) +
  labs(x = element_blank(), y = "Samples")
ggsave("Figures/Prevotella_Ralstonia_disease.png")
count_prevalence(tidy_genus, genus_study, Study) %>%
  ggplot() +
  geom_col(aes(Study, n, fill = presence)) +
  geom_text(aes(y = pos + 3, x = Study, label = label)) +
  labs(x = element_blank(), y = "Samples", title = genus_study)
ggsave("Figures/Ralstonia_study.png")



# * family level ####
rownames(fam2) <- family_tax

gen <- filter_prev(comb_prevalence(fam2, meta, c("Time")))
gen_se <- comb_prevalence(fam2, meta, c("Time", "SEX")) %>% filter_prev()
tre <- comb_prevalence(fam2, meta, c("treatment")) %>% filter_prev()
se <- comb_prevalence(fam2, meta, c("SEX")) %>% filter_prev()
std <- comb_prevalence(fam2, meta, c("Study")) %>% filter_prev()
ibdm <- comb_prevalence(fam2, meta, c("IBD", "Time", "ileum")) %>% filter_prev()
ibd <- comb_prevalence(fam2, meta, c("IBD", "ileum")) %>% filter_prev()
loc <- comb_prevalence(fam2, meta, c("ileum")) %>% filter_prev()
write.csv(ibd, "data_out/prevalence_disease_family_location.csv", row.names = TRUE)


(p <- extract_genus(full_prevalence(fam2, meta, "Time")))
(p <- extract_genus(full_prevalence(fam2, meta, "ileum")))
(genus_ibd <- extract_genus(full_prevalence(fam2, meta, "IBD")))
(genus_study <- extract_genus(full_prevalence(fam2, meta, "Study")))

meta2 <- meta
meta2$ti <- paste(meta$Time, meta$ileum, sep = " & ")
meta2$tI <- paste(meta$Time, meta$IBD, sep = " & ")
meta2$t3 <- paste(meta$Time, meta$Study, sep = " & ")
meta2$t4 <- paste(meta$Study, meta$IBD, sep = " & ")
meta2$t5 <- paste(meta$ileum, meta$IBD, sep = " & ")
meta2$t6 <- paste(meta$ileum, meta$Study, sep = " & ")
(family_time_ileum <- extract_genus(full_prevalence(fam2, meta2, "ti")))
(family_time_ibd <- extract_genus(full_prevalence(fam2, meta2, "tI")))
(family_time_t3 <- extract_genus(full_prevalence(fam2, meta2, "t3")))
(family_time_t4 <- extract_genus(full_prevalence(fam2, meta2, "t4")))
(family_time_t5 <- extract_genus(full_prevalence(fam2, meta2, "t5")))
(family_time_t6 <- extract_genus(full_prevalence(fam2, meta2, "t6")))


(time_fam <- extract_genus(gen))
(time_sex_fam <- extract_genus(gen_se))
(tre_fam <- extract_genus(tre))
(se_fam <- extract_genus(se))
(std_fam <- extract_genus(std))
(ibdm_fam <- extract_genus(ibdm))
(ibd_fam <- extract_genus(ibd))
(loc_fam <- extract_genus(loc))
