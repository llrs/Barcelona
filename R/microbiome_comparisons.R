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


seqtab.nochim <- readRDS("data/ASV.RDS")

ASV <- colnames(seqtab.nochim)
counts_ASV <- seqtab.nochim
colnames(counts_ASV) <- NULL

tab <- t(counts_ASV)
# Remove trailing numbers
colnames(tab) <- gsub("_S.*", "", colnames(tab))
colnames(tab) <- gsub("_p.*", "", colnames(tab))

# Group by tax level ####
{
  microorganism <-  readRDS("data_out/taxonomy_ASV.RDS")$tax
  microorganism <- cbind(microorganism, "ASV" = rownames(microorganism))
  microorganism <- as.data.frame(microorganism, stringsAsFactors = FALSE)
  microorganism$rowname <- as.character(seq_len(nrow(microorganism)))
  rownames(microorganism) <- as.character(seq_len(nrow(microorganism)))
  mit <- rownames(microorganism)[microorganism[, "Family"] == "Mitochondria"]

  # Collapse the data so that there is one Genus per sample, not multiple
  tab2 <- tab
  colnames(tab2) <- seq_along(colnames(tab))
  tab2 <- cbind(tab2, microorganism)

  groups <- syms(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"))
  rm_micro <- tab2 %>%
    group_by(!!!groups) %>%
    summarize(across(where(is.numeric), sum)) %>%
    as.data.frame()
  # Change the first number according to the number of columns on group_by +1
  colnames(rm_micro)[(length(groups) + 1):ncol(rm_micro)] <- colnames(tab)
  counts <- rm_micro[(length(groups) + 1):ncol(rm_micro)]
  microorganism <- rm_micro[, 1:(length(groups))] # And here the same number

  # From the QC step
  meta <- readRDS("data_out/info_samples.RDS")
  colnames(counts) <- gsub("\\.[0-9]", "", colnames(counts))
  meta <- meta[match(colnames(counts), meta$Name), ]

  # filter (keep in mind that it should be on the same order)
  if (!all(colnames(counts) == meta$Name)) {
    stop("Reorder the samples to match the condition!!")
  }
  bcn <- counts[, meta$Study %in% c("BCN", "Controls")]
  meta <- meta[meta$Study %in% c("BCN", "Controls"), ]
  meta$Counts <- colSums(tab)

  # Remove duplicate samples
  replicates <- table(meta$Original)
  replicate_samples <- meta[meta$Original %in% names(replicates[replicates > 1]), ]

  ## But keeping those more sequenced
  keepDup <- replicate_samples %>%
    group_by(Original) %>%
    filter(Counts == max(Counts)) %>%
    arrange(desc(abs(Counts))) %>%
    ungroup()

  nam <- c(names(replicates[replicates == 1]), keepDup$Name)
  otus <- bcn[, colnames(bcn) %in% nam]

  # Working with RNAseq
  # From Juanjo: The original counts are ok, but I need to remove the reseq samples as they
  # have different length and bias the PCA
  conn <- gzfile("data/TNF.all.samples.original.counts.tsv.gz")
  rna <- read.table(conn, sep = "\t", check.names = FALSE)

  rna <- rna[ , !grepl(" reseq$", colnames(rna))] # Remove as said
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


  # Filter the samples ####
  rna2 <- rna[, colnames(rna) %in% meta$Original]
  meta2 <- droplevels(meta[meta$Original %in% colnames(rna2), ])
  OTUs2 <- otus[, colnames(otus) %in% meta2$Name]

  colnames(OTUs2) <- meta2$Original[match(colnames(OTUs2), meta2$Name)]

  # Reorder samples to match!
  meta2 <- meta2[match(colnames(rna2), meta2$Original), ]
  meta3 <- readRDS("data_out/refined_meta_wo_out.RDS")
  meta3$IBD <- as.character(meta3$IBD)
  meta3$IBD[meta3$IBD == "CONTROL"] <- "C"
  # If using wo outliers we removed 2 samples
  # stopifnot(all(meta3$Original == meta2$Original))
  meta2 <- meta3
  meta <- meta2
  rna2 <- rna2[, match(meta$Original, colnames(rna2))]
  OTUs2 <- OTUs2[, match(colnames(rna2), colnames(OTUs2))]
  rownames(meta) <- meta$Original
  OTUs3 <- OTUs2
  OTUs3$`Sample name` <- microorganism$Genus
}
if (FALSE) {
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
}
# Abundance ####
# * Family level ####
if (FALSE) {
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
           str_detect(Sample, "^C|^([0-9]+-w)"),
         Count != 0) %>%
  group_by(Sample) %>%
  mutate(ratio = Count/sum(Count)) %>%
  ungroup() %>%
  left_join(meta, by = c("Sample" = "Name")) %>%
  mutate_if(is.factor, as.character) %>%
  # mutate(IBD = case_when(grepl(x = Name, pattern = "^C") & is.na(IBD) ~ "CONTROL",
  #                        TRUE ~ IBD)) %>%
  arrange(IBD, Time, SEX) %>%
  mutate(Sample = Original,
         Activity = ifelse(is.na(Activity), "INACTIVE", Activity),
         Ileum = ifelse(Exact_location == "ileum", "ileum", "colon"),
         IBD = fct_recode(IBD, C = "CONTROL")) %>%
  mutate(IBD = fct_relevel(IBD, c("C", "UC", "CD")),
         Time = fct_relevel(Time, c("C", "0", "14", "46")),
         Activity = fct_relevel(Activity, c("INACTIVE", "ACTIVE")))
theme_set(theme_minimal())
}
# Several plots
if (FALSE) {
p1 <- tidy_family %>%
  ggplot() +
  geom_point(aes(Sample, Family, size = ratio, col = IBD)) +
  theme(axis.text.x = element_blank())
print(p1)
p2 <- tidy_family %>%
  arrange(IBD) %>%
  ggplot() +
  geom_col(aes(Sample, ratio, fill = Family)) +
  theme(axis.text.x = element_blank()) +
  labs(fill = element_blank()) +
  guides(fill = FALSE) +
  scale_y_continuous(labels = scales::percent, expand = expansion())
print(p2)

p3 <- tidy_family %>%
  filter(Family == "Enterobacteriaceae") %>%
  ggplot() +
  geom_jitter(aes(IBD, ratio, shape = IBD)) +
  labs(x = element_blank(), y = "Beta diversity") +
  theme_minimal()
print(p3)

p4 <- tidy_family %>%
  filter(Family == "Streptococcaceae") %>%
  ggplot() +
  geom_jitter(aes(IBD,ratio, shape = IBD)) +
  labs(x = element_blank(), y = "Beta diversity") +
  theme_minimal()
print(p4)
p5 <- tidy_family %>%
  filter(Family == "Streptococcaceae") %>%
  ggplot() +
  geom_jitter(aes(IBD, ratio, col = Time, shape = IBD)) +
  labs(x = element_blank(), y = "Beta diversity", title = "Streptococcaceae") +
  theme_minimal() +
  facet_nested(~Time + IBD,
               scales = "free_x", switch = "x", nest_line = TRUE) +
  theme(axis.text.x = element_blank())
print(p5)
p6 <- tidy_family %>%
  filter(Family == "Streptococcaceae") %>%
  ggplot() +
  geom_jitter(aes(IBD, ratio, shape = IBD)) +
  labs(x = element_blank(), y = "Beta diversity", title = "Streptococcaceae") +
  theme_minimal() +
  facet_nested(~ IBD + Activity,
               scales = "free_x", switch = "x", nest_line = TRUE) +
  theme(axis.text.x = element_blank())
print(p6)
p7 <- tidy_family %>%
  filter(Family == "Enterobacteriaceae") %>%
  ggplot() +
  geom_jitter(aes(IBD, ratio, col = Time, shape = IBD)) +
  labs(x = element_blank(), y = "Beta diversity", title = "Enterobacteriaceae") +
  theme_minimal() +
  facet_nested(~Time + IBD,
               scales = "free_x", switch = "x", nest_line = TRUE) +
  theme(axis.text.x = element_blank())
print(p7)
}
# PDFs prevalences
if (FALSE) {
tidy_family$IBD <- as.character(tidy_family$IBD)
ts <- tidy_family %>%
  mutate(IBD = case_when(is.na(IBD) ~ "C",
                         TRUE ~ IBD)) %>%
  nest_by(Family, .key = "data_plots") %>%
    mutate(
      Time_IBD = list(
        ggplot(data = data_plots) +
          geom_boxplot(aes(IBD, ratio), alpha = 0, outlier.size = -1) +
          geom_jitter(aes(IBD, ratio, shape = IBD, col = IBD), height = 0) +
          labs(x = element_blank(), y = "Beta diversity", title = Family) +
          theme_minimal() +
          facet_nested(~ Time + IBD,
                       scales = "free_x", switch = "x", nest_line = TRUE) +
          theme(axis.text.x = element_blank()) +
          scale_y_continuous(labels = scales::percent, limits = c(0, NA))),
      IBD_activity = list(
        ggplot(data = data_plots) +
          geom_boxplot(aes(IBD, ratio), alpha = 0, outlier.size = -1) +
          geom_jitter(aes(IBD, ratio, shape = IBD), height = 0) +
          labs(x = element_blank(), y = "Beta diversity", title = Family) +
        theme_minimal() +
        facet_nested(~ IBD + Activity,
                     scales = "free_x", switch = "x", nest_line = TRUE) +
        theme(axis.text.x = element_blank()) +
        scale_y_continuous(labels = scales::percent, limits = c(0, NA))),
    Ileum_IBD_activity = list(
      ggplot(data = data_plots) +
        geom_boxplot(aes(IBD, ratio), alpha = 0, outlier.size = -1) +
        geom_jitter(aes(IBD, ratio, shape = Ileum), height = 0) +
        labs(x = element_blank(), y = "Beta diversity", title = Family) +
        theme_minimal() +
        facet_nested(~ Ileum + IBD + Activity,
                     scales = "free_x", switch = "x", nest_line = TRUE) +
        theme(axis.text.x = element_blank()) +
        scale_y_continuous(labels = scales::percent, limits = c(0, NA)))
  )

ts2 <- ts %>%
  # ungroup() %>%
  # filter(n_distinct(data_plots$IBD)>=2 &
  #          n_distinct(data_plots$Time) >= 2 &
  #          n_distinct(data_plots$Activity) >= 2) %>%
  mutate(IBD_test = ifelse(n_distinct(data_plots$IBD) >= 2,
                kruskal.test(ratio ~ IBD, data = data_plots)$p.value,
              NA),
    time_test =
      ifelse(n_distinct(data_plots$Time) >= 2,
             cor.test( ~ ratio + as.numeric(Time), data = data_plots,
                       method = "spearman")$p.value,
                NA),
    activity_test =
      ifelse(n_distinct(data_plots$Activity) >= 2,
                  kruskal.test(ratio ~ Activity, data = data_plots)$p.value,
                NA),
    Location_test =
      ifelse(n_distinct(data_plots$Ileum) >= 2,
             kruskal.test(ratio ~ Ileum, data = data_plots)$p.value,
             NA),
    IBD_activity_test =
      ifelse(all(group_by(data_plots, Activity, IBD) %>% count() %>% pull(n) >= 2) &
               length(unique(paste0(data_plots$IBD, data_plots$Activity))) != 1,
             kruskal.test(ratio ~ paste0(IBD, Activity), data = data_plots)$p.value,
             NA),
    IBD_loc_test =
      ifelse(all(group_by(data_plots, Ileum, IBD) %>% count() %>% pull(n) >= 2) &
               length(unique(paste0(data_plots$IBD, data_plots$Ileum))) != 1,
             kruskal.test(ratio ~ paste0(IBD, Ileum), data = data_plots)$p.value,
             NA))
df <- ts2 %>%
  ungroup() %>%
  dplyr::select(Family, ends_with("_test"))
colnames(df)[2:ncol(df)] <- c("IBD (p.val)", "Time (p.val)", "Activity (p.val)",
                       "Location (p.val)",  "IBD & Activity (p.val)",
                       "IBD & Location (p.val)")
df2 <- apply(df[, 2:ncol(df)], 2, p.adjust)
colnames(df2) <- gsub(pattern = "p.val", replacement = "adj.p.val", x = colnames(df2))
df3 <- cbind(df, df2)
write.csv(x = df3, file = "data_out/abundance_p_values.csv",
          row.names = FALSE)
pdf("Figures/Time_IBD_abundance.pdf")
for (i in seq_len(nrow(ts))) {
  print(ts$Time_IBD[[i]])
}
dev.off()
pdf("Figures/IBD_activity_abundance.pdf")
for (i in seq_len(nrow(ts))) {
  print(ts$IBD_activity[[i]])
}
dev.off()
pdf("Figures/Ileum_IBD_activity_abundance.pdf")
for (i in seq_len(nrow(ts))) {
  print(ts$Ileum_IBD_activity[[i]])
}
dev.off()


p <- tidy_family %>%
  filter(Family == "Acidaminococcaceae") %>%
  mutate(cat = paste(Ileum, IBD, Activity, sep = "_")) %>%
  mutate(cat = fct_relevel(cat,
                           c("colon_CONTROL_INACTIVE", "colon_UC_INACTIVE",
                             "colon_UC_ACTIVE", "colon_CD_INACTIVE",
                             "colon_CD_ACTIVE", "ileum_CONTROL_INACTIVE",
                             "ileum_CD_INACTIVE", "ileum_CD_ACTIVE"))) %>%
  ggplot(aes(cat, ratio)) +
  geom_boxplot() +
  stat_compare_means() +
  labs(x = element_blank(), y = "Beta diversity", title = "Acidaminococcaceae") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  facet_nested(~ Ileum + IBD + Activity,
               scales = "free", switch = "x", nest_line = TRUE) +
  theme(axis.text.x = element_blank())
print(p)
p <- tidy_family %>%
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
print(p)
}
# * genus level ####
if (FALSE) {
tidy_genus <- OTUs3 %>%
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
genus <- microorganism[, "Genus"]

}

fd <- AnnotatedDataFrame(microorganism)
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
write.csv(tt, "data_out/colon_vs_ileum_ASV.csv", row.names = TRUE)

# Prevalence ####
# Functions
{
filter_prev <- function(x) {
  x <- x[apply(x, 1, function(y){any(y != 1)}), , drop = FALSE]
  x[, apply(x, 2, function(y){any(y != 1)}), drop = FALSE]
}


extract_rownames <- function(x) {
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
}
#### * Family level ####
{
rownames(family) <- family$`Sample name`
fam <- family[, -1]
colnames(fam) <- gsub("_S.+$", "", colnames(fam))
fam <- as.matrix(fam[, meta$Original])

fam_t <- fam %>%
  comb_prevalence(meta, c("Time")) %>%
  filter_prev()
fam_ts <- fam %>%
  comb_prevalence(meta, c("Time", "SEX")) %>%
  filter_prev()
fam_tre <- fam %>%
  comb_prevalence(meta, c("treatment")) %>%
  filter_prev()
fam_se <- fam %>%
  comb_prevalence(meta, c("SEX")) %>%
  filter_prev()
fam_std <- fam %>%
  comb_prevalence(meta, c("Study")) %>%
  filter_prev()
fam_ibdm <- fam %>%
  comb_prevalence(meta, c("IBD", "Time", "ileum")) %>%
  filter_prev()
fam %>%
  comb_prevalence(meta, c("IBD", "ileum")) %>%
  filter_prev() %>%
  extract_rownames()
fam %>%
  comb_prevalence(meta, c("ileum")) %>%
  filter_prev() %>%
  extract_rownames()
fam %>%
  full_prevalence(meta, c("Time")) %>%
  extract_rownames()
fam %>%
  full_prevalence(meta, c("treatment")) %>%
  extract_rownames()
fam %>%
  full_prevalence(meta, c("SEX")) %>%
  extract_rownames()
fam %>%
  full_prevalence(meta, c("Study")) %>%
  extract_rownames()
fam %>%
  full_prevalence(meta, "IBD") %>%
  extract_rownames()
fam_ibd <- fam %>%
  full_prevalence(meta, "ileum") %>%
  filter_prev()
fam_loc <- fam %>%
  full_prevalence(meta, c("ileum")) %>%
  filter_prev()
}


#### * Genus level ####
{
rownames(OTUs2) <- genus

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
genus_ibd <- extract_rownames(full_prevalence(OTUs2, meta, "IBD"))
genus_study <- extract_rownames(full_prevalence(OTUs2, meta, "Study"))

meta2 <- meta
meta2$ti <- paste(meta$Time, meta$ileum, sep = " & ")
meta2$tI <- paste(meta$Time, meta$IBD, sep = " & ")
meta2$t3 <- paste(meta$Time, meta$Study, sep = " & ")
meta2$t4 <- paste(meta$Study, meta$IBD, sep = " & ")
meta2$t5 <- paste(meta$ileum, meta$IBD, sep = " & ")
meta2$t6 <- paste(meta$ileum, meta$Study, sep = " & ")
genus_time_ileum <- extract_rownames(full_prevalence(OTUs2, meta2, "ti"))
genus_time_ibd <- extract_rownames(full_prevalence(OTUs2, meta2, "tI"))
genus_time_t3 <- extract_rownames(full_prevalence(OTUs2, meta2, "t3"))
genus_time_t4 <- extract_rownames(full_prevalence(OTUs2, meta2, "t4"))
genus_time_t5 <- extract_rownames(full_prevalence(OTUs2, meta2, "t5"))
genus_time_t6 <- extract_rownames(full_prevalence(OTUs2, meta2, "t6"))


(time_genus <- extract_rownames(gen))
(time_sex_genus <- extract_rownames(gen_se))
(tre_genus <- extract_rownames(tre))
(se_genus <- extract_rownames(se))
(std_genus <- extract_rownames(std))
(ibdm_genus <- extract_rownames(ibdm))
(ibd_genus <- extract_rownames(ibd))
(loc_genus <- extract_rownames(loc))

# So basically I need to plot for ibdm and ibd_genus
{
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
}
}
# * Family level ####
rownames(fam2) <- family_tax
{
gen <- filter_prev(comb_prevalence(fam2, meta, c("Time")))
gen_se <- comb_prevalence(fam2, meta, c("Time", "SEX")) %>% filter_prev()
tre <- comb_prevalence(fam2, meta, c("treatment")) %>% filter_prev()
se <- comb_prevalence(fam2, meta, c("SEX")) %>% filter_prev()
std <- comb_prevalence(fam2, meta, c("Study")) %>% filter_prev()
ibdm <- comb_prevalence(fam2, meta, c("IBD", "Time", "ileum")) %>% filter_prev()
ibd <- comb_prevalence(fam2, meta, c("IBD", "ileum")) %>% filter_prev()
loc <- comb_prevalence(fam2, meta, c("ileum")) %>% filter_prev()
write.csv(ibd, "data_out/prevalence_disease_family_location.csv", row.names = TRUE)
}
(p <- extract_genus(full_prevalence(fam2, meta, "Time")))
(p <- extract_genus(full_prevalence(fam2, meta, "ileum")))
(genus_ibd <- extract_genus(full_prevalence(fam2, meta, "IBD")))
(genus_study <- extract_genus(full_prevalence(fam2, meta, "Study")))
{
meta2 <- meta
meta2$ti <- paste(meta$Time, meta$ileum, sep = " & ")
meta2$tI <- paste(meta$Time, meta$IBD, sep = " & ")
meta2$t3 <- paste(meta$Time, meta$Study, sep = " & ")
meta2$t4 <- paste(meta$Study, meta$IBD, sep = " & ")
meta2$t5 <- paste(meta$ileum, meta$IBD, sep = " & ")
meta2$t6 <- paste(meta$ileum, meta$Study, sep = " & ")
meta2$t7 <- paste(meta$IBD, meta$Activity, sep = " & ")
(family_time_ileum <- extract_genus(full_prevalence(fam2, meta2, "ti")))
(family_time_ibd <- extract_genus(full_prevalence(fam2, meta2, "tI")))
(family_time_t3 <- extract_genus(full_prevalence(fam2, meta2, "t3")))
(family_time_t4 <- extract_genus(full_prevalence(fam2, meta2, "t4")))
(family_time_t5 <- extract_genus(full_prevalence(fam2, meta2, "t5")))
(family_time_t6 <- extract_genus(full_prevalence(fam2, meta2, "t6")))
(family_time_t7 <- extract_genus(full_prevalence(fam2, meta2, "t7")))
}

(time_fam <- extract_genus(gen))
(time_sex_fam <- extract_genus(gen_se))
(tre_fam <- extract_genus(tre))
(se_fam <- extract_genus(se))
(std_fam <- extract_genus(std))
(ibdm_fam <- extract_genus(ibdm))
(ibd_fam <- extract_genus(ibd))
(loc_fam <- extract_genus(loc))
