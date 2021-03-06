library("dplyr")
library("phyloseq")
library("ggplot2")
library("decontam")
library("tidyr")
library("vegan")
library("forcats")

tab <- read.delim("data/Partek_Michigan3_Kraken_Classified_family.tsv", check.names = FALSE)
colnames(tab) <- gsub("_S.*", "", colnames(tab))
counts <- tab[, -1]
family <- tab[, 1, FALSE]
seqs <- colSums(counts)
not_empty <- function(x) {
  sum(x != 0)
}
sp <- apply(counts, 2, not_empty)
ctrls <- grepl("^500_", colnames(counts))


seqs <- read.delim("output/reads.txt", sep = " ", header = FALSE)
colnames(seqs) <- c("file", "reads")
fs <- strcapture("(([0-9]{3}-w[0-9]{3}|C[0-9]+-T.+[GAIH]|[0-9]+-?[DS]?0?-T.*[GAIH]|500)_?(p[0-9][A-Z]?[0-9]*)?_(S[0-9]+)?_L001_R[12]_001.fastq.gz*)", seqs$file,
                 proto = data.frame(file = character(),
                                    sample = character(),
                                    replicate = character(),
                                    S = character(),
                                    stringsAsFactors = FALSE))
fs2 <- merge(fs, seqs, all = TRUE)
stopifnot(nrow(seqs) == nrow(fs2))
f <- function(x){
  if (any(!is.na(x))) {
    plates <- paste0("p", 1:4)
    x[is.na(x)] <- plates[!plates %in% x]
  }
  x
}
fs3 <- fs2 %>%
  mutate(Study = case_when(
    startsWith(sample, "500") ~ "water",
    startsWith(sample, "C") ~ "Controls",
    grepl("-w", sample) ~ "BCN",
    TRUE ~ "TRIM")) %>%
  mutate(replicate = ifelse(replicate == "", NA, replicate)) %>%
  group_by(sample) %>%
  mutate(replicate= f(replicate))


# For some reason (too big?) this is failing
# library("vegan")
# S <- specnumber(counts)
# # min(rowSums(counts))
# Srar <- rarefy(counts, sample = 1)
# plot(S, Srar, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
# abline(0, 1)
# r <- rarecurve(counts, step = 1, sample = 1, col = "blue")


# Contamination by Escherichia
o <- prop.table(as.matrix(counts), margin = 2)
# o[grep("scheri", family, ignore.case = TRUE), ]*100
# No families from Escherichia coli

sp_ctrls <- as.logical(apply(counts[, ctrls], 1, not_empty))

seqs_b <- colSums(counts[!sp_ctrls, ])
sp_b <- apply(counts[!sp_ctrls, ], 2, not_empty)
info_seq <- cbind(sp, seqs, sp_b, seqs_b)


family.tidy <- gather(tab, Sample, Count, -'Sample name') %>%
  filter(Count != 0) %>%
  dplyr::rename(Microorganism = "Sample name") %>%
  group_by(Sample) %>%
  mutate(ratio = Count/sum(Count))

# Input the disease, state and so on
meta <- read.csv("data/batches.csv")
meta <- meta %>%
  dplyr::select(-Study, -ID, -batch)

# Input about the concentration used and location on the plates
# From the post-sequencing.R files
location <- readRDS("output/Samples_concentration_distribution.RDS")
# Add the study and separate name
location <- location %>%
  mutate(
    Study = case_when(
      grepl("-w", Name) ~ "BCN",
      grepl("^500", Name) ~ "NC",
      grepl("^C", Name) ~ "Controls",
      TRUE ~ "TRIM"
    )
  ) %>%
  separate(Name, c("Original", "Replicate"), sep = "_", remove = FALSE)
# Expected a warning of filling pieces

# Merge the different data set
colnames(info_seq) <- c("Families", "Counts", "Families_wo_ctrls", "Counts_wo_ctrls")
out <- merge(location, info_seq, by.x = "Name", by.y = "row.names", all = TRUE)
out <- merge(out, meta, by.x = "Original", by.y = "Sample.Name_RNA", all.x = TRUE)
out <- out[match(colnames(counts), out$Name), ] # Reorder


all_data <- merge(family.tidy, out, by.x = "Sample", by.y = "Name")
# Plot about the number of counts per sample and the column of origin
all_data %>%
  group_by(Sample) %>%
  summarise(counts = sum(Count),
            study = unique(Study),
            Plate = unique(Plate),
            Row = unique(Row),
            Col = unique(Column)
            ) %>%
  ungroup() %>%
  droplevels() %>%
  ggplot() +
  geom_col(aes(lvls_reorder(Sample, order(counts)), log10(counts), fill = Col)) +
  labs(x = "Samples") +
  geom_hline(yintercept = c(log10(1000), log10(20000)), col = c("red", "green")) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank())

all_data %>% as_tibble() %>%
  filter(Study == "BCN" & Microorganism != "no_match")

topMicro <- all_data %>%
  filter(Study == "BCN") %>%
  filter(Microorganism != "no_match") %>%
  group_by(Microorganism) %>%
  summarise(mA = mean(ratio)) %>%
  arrange(desc(mA)) %>%
  top_n(10, mA)
theme_update(strip.background = element_blank())

# Abundance plot about the top microorganisms on the three times
all_data %>%
  filter(Study == "BCN",
         Microorganism %in% topMicro$Microorganism) %>%
  separate(Original, into = c("Patient", "week"), sep = "-w",
           convert = TRUE, remove = FALSE) %>%
  mutate(week = ifelse(week == 38, 46, week),
         Patient = as.factor(Patient)) %>%
  droplevels() %>%
  ggplot() +
  geom_tile(aes(Patient, Microorganism,
                fill = ratio*100)) +
  scale_fill_viridis_c() +
  facet_wrap(~week) +
  labs(fill = "Abundance (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        axis.ticks.x = element_blank())


ggplot(all_data) +
  geom_point(aes(Sample, Microorganism, size = ratio, col = IBD)) +
  theme(panel.grid = element_blank())

out$Concentr[is.na(out$Concentr)] <- 0.0005
out$Exact_location[out$Exact_location == ""] <- NA
out <- droplevels(out)
saveRDS(out, "output/info_samples.RDS")

# Rename the columns and rows
colnames(counts) <- paste0("a", as.character(seq_len(ncol(counts))))
rownames(out) <- paste0("a", as.character(seq_len(nrow(out))))
rownames(family) <- as.character(seq_len(nrow(family)))
colnames(family) <- c("Family")

phyloseq <- phyloseq(otu_table(counts, taxa_are_rows = TRUE),
                     sample_data(out)
                     # tax_table(as.matrix(family))
                     )
theme_set(theme_bw())
beta <- estimate_richness(phyloseq)
res <- cbind(beta, out)
(p <- plot_richness(phyloseq, "Original", "Plate", measures = "Shannon"))
p + theme(panel.grid.major.x = element_blank()) + labs(x = element_blank())
ggsave("Figures/alpha_diversity_plates.png")
(p <- plot_richness(phyloseq, "Name", "Row", measures = "Shannon"))
(p <- plot_richness(phyloseq, "Name", "Column", measures = "Shannon"))
(p <- plot_richness(phyloseq, "Name", "Concentr", measures = "Shannon"))

ggplot(res) +
  geom_point(aes(Shannon, Concentr, col = Column))
ggplot(res) +
  geom_point(aes(Shannon, Concentr, col = Row))
ggplot(res) +
  geom_point(aes(Shannon, Concentr, col = Plate))
ggplot(res) +
  geom_point(aes(Exact_location, log10(Families), col = Study))
ggplot(res) +
  geom_point(aes(Exact_location, log10(Counts), col = Study))
ggplot(res) +
  geom_point(aes(Plate, log10(Counts), col = Exact_location))
ggplot(res) +
  geom_point(aes(Row, log10(Counts), col = Exact_location))
ggplot(res) +
  geom_point(aes(Column, log10(Counts), col = Exact_location))
ggplot(res) +
  geom_point(aes(log10(Counts), Shannon, col = Exact_location))

# Batch is each plate

count_decontam <- t(as.matrix(counts))
colnames(count_decontam) <- family$Family
decontamination <- isContaminant(count_decontam, neg = out$Study %in% "NC",
                                 method = "either", conc = out$Concentr,
                                 batch = out$Plate)

# I would remove this sequences (also those that are only present in one sample)
contam <- is.na(decontamination$p.prev) | decontamination$contaminant
summary(contam)

# Batch is each study
decontamination2 <- isContaminant(count_decontam, neg = out$Study %in% "NC",
                                 method = "freq", conc = out$Concentr,
                                 batch = out$Study)

# I would remove this sequences (also those that are only present in one sample)
contam2 <- is.na(decontamination2$p.prev) | decontamination2$contaminant
summary(contam2)
# Conclusion, no contaminants if we take into account the different batches

contaminant_families <- as.character(family)[contam]
metagSeq <- phyloseq_to_metagenomeSeq(phyloseq)
