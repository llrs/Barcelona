library("dada2")
library("dplyr")
library("tidyr")
library("phyloseq")
library("ggplot2")

# Code from http://bioconductor.org/packages/3.10/bioc/vignettes/dada2/inst/doc/dada2-intro.html
# And http://benjjneb.github.io/dada2/tutorial.html
samplesR1 <- list.files(path = "data/fastq_ASV",  pattern = "_R1_", full.names = TRUE)
samplesR2 <- list.files(path = "data/fastq_ASV",  pattern = "_R2_", full.names = TRUE)

files <- read.table("output/reads.txt", sep = " ", row.names = NULL)
keep_files <- files$V1[files$V2 > 3000 & files$V2 < 2*sd(files$V2)+mean(files$V2)]
samplesR1 <- samplesR1[samplesR1 %in% paste0("data/fastq_ASV/", keep_files)]
samplesR2 <- samplesR2[samplesR2 %in% paste0("data/fastq_ASV/", keep_files)]
plotQualityProfile(samplesR1[1]) # Forward
plotQualityProfile(samplesR2[1]) # Forward


tempdir_f1 <- tempdir()
tempdir_r1 <- tempdir()
filtF1 <- paste0(tempdir_f1, "/filt_", basename(samplesR1))
filtR1 <- paste0(tempdir_r1, "/filt_", basename(samplesR2))

# filtF1 <- tempfile(fileext = ".fastq.gz")
# filtR1 <- tempfile(fileext = ".fastq.gz")
out <- filterAndTrim(fwd = samplesR1, filt = filtF1,
                     rev = samplesR2, filt.rev = filtR1,
                     trimLeft = 10, truncLen = c(240, 200),
                     maxN = 0, maxEE = 2, matchIDs = TRUE,
                     compress = TRUE, verbose = FALSE, multithread = 6)

derepF1 <- derepFastq(filtF1, verbose = FALSE)
derepR1 <- derepFastq(filtR1, verbose = FALSE)

errF <- learnErrors(derepF1, multithread = TRUE)
errR <- learnErrors(derepR1, multithread = TRUE)

plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)

dadaF1 <- dada(derepF1, err = errF, multithread = TRUE)
dadaR1 <- dada(derepR1, err = errR, multithread = TRUE)

merger1 <- mergePairs(dadaF1, derepF1, dadaR1, derepR1, verbose = FALSE)

seqtab <- makeSequenceTable(merger1)
seqtab.nochim <- removeBimeraDenovo(seqtab, verbose = FALSE)

# The exact same results as without matchID = TRUE on filterAndTrim
# saveRDS(seqtab.nochim, file = "data/ASV_matchID.RDS")
saveRDS(seqtab.nochim, file = "output/ASV_cut3000_up.RDS")
seqtab.nochim <- readRDS("output/ASV_cut3000_up.RDS")
# seqtab.nochim <- readRDS("data/ASV.RDS")

ASV <- colnames(seqtab.nochim)
counts_ASV <- seqtab.nochim
colnames(counts_ASV) <- NULL

ASV_counts <- t(counts_ASV)
cASV <- sort(colSums(ASV_counts))
barplot(log10(cASV))
abline(h = c(log10(500), log10(median(cASV))), col = c("red", "green"))

sample_names <- gsub("_L.*", "", gsub("filt_", "", colnames(ASV_counts)))
s <- strcapture(pattern = "([0-9]{3}-w[0-9]{3}|C[0-9]-T-DM-.+)_?(p[0-9])?_(S[0-9]+)",
           x = sample_names,
           proto = data.frame(sample = character(),
                              replicate = character(),
                              s = character(),
                              stringsAsFactors = FALSE))
s$file <- colnames(ASV_counts)
s$replicate[!nzchar(s$replicate)] <- NA
s$reads <- colSums(ASV_counts)
s2 <- s %>%
  group_by(sample) %>%
  filter(reads == max(reads)) %>%
  ungroup()
ASV_count_f <- ASV_counts[, s2$file]
colnames(ASV_count_f) <- s2$sample
ASV_count_f <- ASV_count_f[rowSums(ASV_count_f) != 0, ]

meta <- readRDS("output/refined_meta_all.RDS")
meta <- meta[match(colnames(ASV_count_f), meta$Original), ]
rownames(meta) <- colnames(ASV_count_f)
meta$Loc <- ifelse(meta$Exact_location == "ileum", "ileum", "colon")
meta$Activity <- ifelse(is.na(meta$Activity), "INACTIVE", meta$Activity)

phyloseq <- phyloseq(otu_table(ASV_count_f, taxa_are_rows = TRUE),
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
  labs(y = "Alpha diversity", x = element_blank(), title = "ASV diversity",
       subtitle = "cutoff 3000 reads") +
  theme_minimal()


richness_rel %>%
  filter(`Alpha diversity` == "Shannon") %>%
  ggplot(aes(Activity, effective, col = IBD)) +
  geom_boxplot(alpha = 0, outlier.size = 0) +
  geom_point(position = position_jitterdodge(jitter.height = 0, jitter.width = 1/4)) +
  facet_wrap( ~ Loc, drop = TRUE) +
  labs(y = "Shannon alpha diversity", x = element_blank(), title = "Diversity ASV",
       subtitle = "cutoff top and bottom reads") +
  theme_minimal() +
  scale_y_continuous(limits = c(0, NA))
ggsave("Figures/cut3000_up_ASV_diversity.png")
