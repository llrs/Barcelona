library("ggplot2")
# load data ####
library("readxl")
library("dplyr")
library("metagenomeSeq")
library("integration")
library("stringr")
library("purrr")
library("ggplot2")
library("lubridate")
library("org.Hs.eg.db")
library("MultiAssayExperiment")
library("MOFA")

{
tab <- read.delim("data/Partek_Michigan3_Kraken_Classified_genus.tsv", check.names = FALSE)
colnames(tab) <- gsub("_S.*", "", colnames(tab)) # Remove trailing numbers
counts <- tab[, -1]
genus <- tab[, 1, FALSE]

# From the QC step
meta <- readRDS("data_out/info_samples.RDS")
meta$Counts <- colSums(counts)
}
{
# filter (keep in mind that it should be on the same order)
if (!all(colnames(counts) == meta$Name)) {
  stop("Reorder the samples to match the condition!!")
}
bcn <- counts[, meta$Study %in% c("BCN", "Controls")]
meta <- meta[meta$Study %in% c("BCN", "Controls"), ]

# Remove duplicate samples
replicates <- table(meta$Original)
replicate_samples <- meta[meta$Original %in% names(replicates[replicates > 1]), ]

## By keeping those more sequenced
keepDup <- replicate_samples %>%
  group_by(Original) %>%
  filter(Counts == max(Counts)) %>%
  arrange(desc(abs(Counts))) %>%
  ungroup()

nam <- c(names(replicates[replicates == 1]), keepDup$Name)
otus <- bcn[, colnames(bcn) %in% nam]
meta <- meta[meta$Name %in% nam, ]
}
{
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


# Filter the samples
rna2 <- rna[, colnames(rna) %in% meta$Original]
meta2 <- droplevels(meta[meta$Original %in% colnames(rna2), ])
OTUs2 <- otus[, colnames(otus) %in% meta2$Name]

colnames(OTUs2) <- meta2$Original[match(colnames(OTUs2), meta2$Name)]

# Reorder samples to match!
meta2 <- meta2[match(colnames(rna2), meta2$Original), ]
OTUs2 <- OTUs2[match(colnames(rna2), colnames(OTUs2))]


OTUs2 <- norm_RNAseq(OTUs2)
rna2 <- norm_RNAseq(rna2)
rna2 <- filter_RNAseq(rna2)

experiments <- list(RNA = rna2, OTUS = OTUs2)
}
mae <- MultiAssayExperiment(
  experiments = experiments,
  # colData = meta2
)
sm <- data.frame(assay = rep(names(experiments), each = nrow(meta2)),
                 primary = rep(meta2$Original, 2),
                 colname = c(colnames(rna2), colnames(OTUs2)))
sampleMap(mae) <- as(sm, "DataFrame")
colData(mae) <- as(meta2[, c("IBD", "SEX", "Exact_location")], "DataFrame")
MOFAobject <- createMOFAobject(mae)
plotDataOverview(MOFAobject)

mofa <- prepareMOFA(MOFAobject)
n_inits <- 3
MOFAobject0 <- MOFAobject

TrainOptions <- getDefaultTrainOptions()
ModelOptions <- getDefaultModelOptions(MOFAobject)
DataOptions <- getDefaultDataOptions()

TrainOptions$DropFactorThreshold <- 0.01

MOFAlist <- lapply(seq_len(n_inits), function(it) {

  TrainOptions$seed <- 2020 + it

  MOFAobject <- prepareMOFA(
    MOFAobject0,
    DataOptions = DataOptions,
    ModelOptions = ModelOptions,
    TrainOptions = TrainOptions
  )

  runMOFA(MOFAobject)
})
