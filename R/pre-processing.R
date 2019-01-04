library("readxl")
library("dplyr")
library("metagenomeSeq")
library("integration")
library("stringr")
library("purrr")

tab <- read.delim("data/Partek_Michigan3_Kraken_Classified_genus.tsv", check.names = FALSE)
colnames(tab) <- gsub("_S.*", "", colnames(tab)) # Remove trailing numbers
counts <- tab[, -1]
genus <- tab[, 1, FALSE]

# From the QC step
meta <- readRDS("info_samples.RDS")
meta$Counts <- colSums(counts)

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

# normalize 16S
OTUs <- norm_otus(otus, genus)


# Working with RNAseq
conn <- gzfile("data/voom.RNAseq.data.all.cal.noduplications.tsv.gz")
rna <- read.table(conn, sep = "\t", check.names = FALSE)
close(conn)

colnames(rna) <- gsub(" reseq$", "", colnames(rna))
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

# Filter them
rna2 <- rna[, colnames(rna) %in% meta$Original]
meta2 <- droplevels(meta[meta$Original %in% colnames(rna2), ])
OTUs2 <- OTUs[, colnames(OTUs) %in% meta2$Name]


# Meta ####

db <- data.table::fread(
  "data/db_biopsies_bcn_seq16S_noTRIM.txt", sep = "\t",
  stringsAsFactors = FALSE
)

meta3 <- merge(meta2, db,
               by.x = "Original", by.y = "Sample_Code",
               all.x = TRUE, all.y = FALSE)
meta3 <- droplevels(meta3)

# From the copy of the access database I have 04/01/2019
dates <- read_xlsx("data/fechas.xlsx",
                   col_types = c("text", "text", "guess", "guess"), na = "")
dates$Date <- as.Date(dates$Date)
dates$`Date of diagnosis` <- as.Date(dates$`Date of diagnosis`)
colnames(dates) <- c("Visit", "ID", "DATE_SAMPLE", "Date_diagnostic")
meta4 <- merge(meta3, dates, by.x = "Original", by.y = "Visit", all.x = TRUE, all.y = FALSE)

diagTime <- as.Date(meta4$DATE_SAMPLE, "%Y-%m-%d") - meta4$Date_diagnostic
diagTime <- as.numeric(diagTime/365.25)
diagTime[is.na(diagTime)] <- 0 # Replace by date of the sample at least
AgeDiag <- as.numeric(meta4$Date_diagnostic - as.Date(meta4$Birth_date,
                                         "%m/%d/%Y"))/365.25
meta5 <- cbind.data.frame(meta4, diagTime, AgeDiag)

A <- list("RNAseq" = t(rna2), "Micro" = t(OTUs2), "Meta" = meta3)
A[1:2] <- clean_unvariable(A[1:2]) # Just the numeric ones
saveRDS(A, "data/RGCCA_data.RDS")
