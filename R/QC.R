library("dplyr")
library("phyloseq")
library("ggplot2")
library("decontam")

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

# For some reason this is failing
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
  rename(Microorganism = 'Sample name') %>%
  group_by(Sample) %>%
  mutate(ratio = Count/sum(Count))

# Input the disease, state and so on
meta <- read.csv("data/batches.csv")
meta <- meta %>%
  select(-Study, -ID, -batch)

# Input about the concentration used and location on the plates
location <- readRDS("Samples_concentration_distribution.RDS")
# Add th study and separate name
location <- location %>%
  mutate(
    Study = case_when(
      grepl("-w", Name) ~ "BCN",
      grepl("^500", Name) ~ "NC",
      grepl("^C", Name) ~ "Controls",
      TRUE ~ "TRIM"
    )
  ) %>%
  tidyr::separate(Name, c("Original", "Replicate"), sep = "_", remove = FALSE)


# Merge the different dataset
colnames(info_seq) <- c("Families", "Counts", "Families_wo_ctrls", "Counts_wo_ctrls")
out <- merge(location, info_seq, by.x = "Name", by.y = "row.names", all = TRUE)
out <- merge(out, meta, by.x = "Original", by.y = "Sample.Name_RNA", all.x = TRUE)
out <- out[match(colnames(counts), out$Name), ] # Reorder


all_data <- merge(family.tidy, out, by.x = "Sample", by.y = "Name")
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
  geom_hline(yintercept = c(log10(1000), log10(20000)), col = c("red", "green"))

all_data %>% as_tibble() %>%
  filter(Study == "BCN") %>%
  filter(Microorganism != "no_match")

topMicro <- all_data %>%
  filter(Study == "BCN") %>%
  filter(Microorganism != "no_match") %>%
  group_by(Microorganism) %>%
  summarise(mA = mean(ratio)) %>%
  arrange(desc(mA)) %>%
  top_n(10, mA)
theme_update(strip.background = element_blank())
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
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

out$Concentr[is.na(out$Concentr)] <- 0.0005
out$Exact_location[out$Exact_location == ""] <- NA
out <- droplevels(out)
saveRDS(out, "info_samples.RDS")

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

contaminant_families <- as.character(family)[contam]
metagSeq <- phyloseq_to_metagenomeSeq(phyloseq)
