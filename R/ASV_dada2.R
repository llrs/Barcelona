library("dada2")

# Code from http://bioconductor.org/packages/3.10/bioc/vignettes/dada2/inst/doc/dada2-intro.html
# And http://benjjneb.github.io/dada2/tutorial.html
samplesR1 <- list.files(path = "data/fastq_ASV/",  pattern = "_R1_", full.names = TRUE)
samplesR2 <- list.files(path = "data/fastq_ASV/",  pattern = "_R2_", full.names = TRUE)
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
# saveRDS(seqtab.nochim, file = "data/ASV.RDS")
seqtab.nochim <- readRDS("data/ASV.RDS")

ASV <- colnames(seqtab.nochim)
counts_ASV <- seqtab.nochim
colnames(counts_ASV) <- NULL

ASV_counts <- t(counts_ASV)
cASV <- sort(colSums(ASV_counts))
barplot(log10(cASV))
abline(h = c(log10(500), log10(median(cASV))), col = c("red", "green"))

# out <- assignTaxonomy(ASV, refFasta = "data/silva_nr99_v138_train_set.fa.gz",
#                       outputBootstraps = TRUE,
#                tryRC = TRUE, multithread = TRUE)
# saveRDS(out, "output/taxonomy_ASV.RDS")
taxonomy <- readRDS("output/taxonomy_ASV.RDS")
## Keep track of read along the pipeline:
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaF1, getN), sapply(dadaR1, getN), sapply(merger1, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")

track %>%
  as.data.frame() %>%
  mutate(samples = rownames(.)) %>%
  pivot_longer(input:nonchim) %>%
  mutate(name = as.factor(name)) %>%
  mutate(name = fct_relevel(name, colnames(track))) %>%
  group_by(samples) %>%
  mutate(per = value/max(value)*100) %>%
  ungroup() %>%
  ggplot() +
  # Reorder samples by the initial amount of ASV
  geom_tile(aes(name, fct_reorder(samples, value, .fun = max), fill = per)) +
  labs(x = "Steps", y = "Samples", title = "Amplicons present",
       fill = "Percentage (%)") +
  theme_minimal() +
  scale_x_discrete(position = "top") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

