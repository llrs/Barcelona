library("readxl")
library("dplyr")
library("metagenomeSeq")
library("integration")
library("stringr")
library("purrr")
library("ggplot2")
library("lubridate")

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
conn <- gzfile("data/TNF.all.samples.original.counts.tsv.gz") # TODO See if this is a good choice
rna <- read.table(conn, sep = "\t", check.names = FALSE)

colnames(rna) <- gsub(" reseq$", "", colnames(rna))
colnames(rna)[grep("[Ww]", colnames(rna))] <- tolower(colnames(rna)[grep("[Ww]", colnames(rna))])

rna <- norm_RNAseq(rna)
rna <- filter_RNAseq(rna)

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
               by.x = c("Original", "IBD"), by.y = c("Sample_Code", "IBD"),
               all.x = TRUE, all.y = FALSE)
meta3 <- droplevels(meta3)

# From the copy of the access database I have 04/01/2019
dates <- read_xlsx("data/fechas.xlsx",
                   col_types = c("text", "text", "guess", "guess"), na = "")
dates$Date <- as.Date(dates$Date)
dates$`Date of diagnosis` <- as.Date(dates$`Date of diagnosis`)
colnames(dates) <- c("Visit", "ID", "DATE_SAMPLE", "Date_diagnostic")
meta4 <- merge(meta3, dates, by.x = "Original", by.y = "Visit", all.x = TRUE, all.y = FALSE)

# Arrange the IDs
meta4$ID[grep("^C", meta4$Original)] <- gsub("^(C[0-9]+)-.+", "\\1", meta4$Original[grep("^C", meta4$Original)])
meta4$ID[grep("^001", meta4$Original)] <- "001"

# Date diagnostic related
meta4$Date_diagnostic[meta4$ID %in% "001"] <- as.Date("26/01/2009", "%d/%m/%Y")
meta4$Date_diagnostic[meta4$ID %in% "017"] <- as.Date("03/06/2013", "%d/%m/%Y")
meta4$Date_diagnostic[meta4$ID %in% "122"] <- as.Date("28/07/2015", "%d/%m/%Y")
meta4$Date_diagnostic[meta4$ID %in% "123"] <- as.Date("10/11/2015", "%d/%m/%Y")

# Use the imported metadata from the TRIM project
# For the controls that are shared
meta_trim <- readRDS("data/metaTRIM.RDS")

meta_c <- meta_trim[grep("^C", meta_trim$ID), ]
meta_c$Sample_Code_uDNA <- gsub("_", "-", meta_c$Sample_Code_uDNA)
order_c <- match(meta4$Original[meta4$Original %in% meta_c$Sample_Code_uDNA], meta_c$Sample_Code_uDNA)
meta_c <- meta_c[order_c, ]
keep_controls <- meta4$Original %in% meta_c$Sample_Code_uDNA
meta4$SEX[keep_controls] <- tolower(meta_c$SEX)

# Different time formats in Birth date!!
meta4$Birth_date <- as.Date(meta4$Birth_date, "%d/%m/%Y")
meta4$Birth_date[keep_controls] <- as.Date(meta_c$Birth_date, "%m/%d/%Y")
meta4$DATE_SAMPLE[keep_controls] <- as.Date(meta_c$DATE_SAMPLE, "%m/%d/%Y")
meta4$Exact_location[keep_controls] <- gsub(" colon", "", tolower(meta_c$Exact_location))

# Patient 017 has a sample without birth date
p17 <- meta4[meta4$ID == "017", c("Original", "ID", "DATE_SAMPLE", "Birth_date", "Age")]
meta4$Birth_date[meta4$ID == "017"] <- p17$Birth_date[1]

ggplot(meta4[grep("^[0-9]", meta4$ID), ]) +
  geom_point(aes(ID, DATE_SAMPLE, color = factor(substr(Original, 6, 8)),
                 size = Age)) +
  labs(x = "ID", y = "Date sample", color = "Sample id (weeks)")
# See the plot!!
# The patient ID 113 has the week 46 before the week 000??
# Also 2 missing DATE_SAMPLES warning message

# Using the other database we can complete and correct the record:
meta4$DATE_SAMPLE[meta4$Original == "001-w000"] <- as.Date("11/14/2012", "%m/%d/%Y")
meta4$DATE_SAMPLE[meta4$Original == "111-w038"] <- as.Date("06/07/2016", "%m/%d/%Y")
meta4$DATE_SAMPLE[meta4$Original == "113-w046"] <- as.Date("08/09/2016", "%m/%d/%Y")

ggplot(meta4[grep("^[0-9]", meta4$ID), ]) +
  geom_point(aes(ID, DATE_SAMPLE, color = factor(substr(Original, 6, 8)),
                 size = Age)) +
  labs(x = "ID", y = "Date sample", color = "Sample id (weeks)",
       title = "After corrections")

# Calculate years since diagnostic and related problems
diagTime <- as.Date(meta4$DATE_SAMPLE, "%Y-%m-%d") - meta4$Date_diagnostic
diagTime <- as.numeric(diagTime/365.25)
diagTime[is.na(diagTime)] <- 0 # Replace by date of the sample at least
AgeDiag <- as.numeric(meta4$Date_diagnostic - as.Date(meta4$Birth_date,
                                         "%d/%m/%Y"))/365.25

# Calculate the age because they are off up to 5 years up
hist(meta4$Age - as.numeric(meta4$DATE_SAMPLE - meta4$Birth_date)/365.25)
meta4$Age <- as.numeric(meta4$DATE_SAMPLE - meta4$Birth_date)/365.25

AgeDiag[is.na(AgeDiag)] <- 0
meta4$diagTime <- diagTime
meta4$AgeDiag <- AgeDiag

db2 <- readxl::read_xls("data/bd_BCN_tnf_biopsies_110119.xls", na = "n.a.")
db2 <- db2[!is.na(db2$RNA_seq_batch), ]
db2$NHC <- as.character(db2$NHC)
db2$Sample_id <- gsub(" reseq", "", db2$Sample_id)
db2$Sample_id <- gsub("-T-TR-", "-T-DM-", db2$Sample_id)
db2$Sample_id[grep("-[wW]", db2$Sample_id)] <- tolower(db2$Sample_id)[grep("-[wW]", db2$Sample_id)]
db2$Sample_id <- str_split(db2$Sample_id, "-w") %>% # Ready for BCN
  map(correct_bcn) %>%
  unlist()
# Based that on the resequenced is after the original
db3 <- db2[-c(which(duplicated(db2$Sample_id) == TRUE) - 1), ]
db3 <- db3[db3$Sample_id %in% colnames(rna2), ]
meta5 <- merge(meta4, db3,
               by.x = c("Original", "IBD", "SEX"),
               by.y = c("Sample_id", "IBD", "Gender"),
               all.x = TRUE)
meta5$sample_date <- as.Date(meta5$sample_date, "%m/%d/%Y")

# Different sample date?? The right one is DATE_SAMPLE
meta5[meta5$sample_date != meta5$DATE_SAMPLE, c("sample_date", "DATE_SAMPLE", "Original")]
meta5 <- meta5[, -grep("sample_date", colnames(meta5))]

duplicates <- group_by(meta5, NHC) %>%
  summarise(diff = n_distinct(Patient_ID))
dupli <- duplicates$NHC[duplicates$diff > 1]

# How do we codify it ?? The one with less  samples loses its Patient_ID
meta5[meta5$NHC %in% dupli, "Patient_ID"] <- c("17", "17", "17", "86", "92", "86", "86", "86", "92", "92")

# Treatment. check with the database? Better with the database
meta5 <- meta5[, -grep("Treatment.x", colnames(meta5))] # Only signaling those that are not controls
treat <- readxl::read_xlsx("data/Treatment.xlsx") # On 17/01/2019 on %D/%M/%Y format
treat <- treat[, c("Visit", "Drug")]
treat <- treat[!is.na(treat$Visit), ]
treat <- treat[treat$Visit %in% meta5$Original, ]
Drugs <- unique(treat$Drug)
visits <- unique(treat$Visit)
incidence <- sapply(visits, function(x){
  drugs <- treat$Drug[treat$Visit == x]
  keep <- as.numeric(Drugs %in% drugs)
  names(keep) <- Drugs
  keep
})
incidence <- t(incidence)
# We don't know the treatment of many samples
meta5$Original[grep("-w", meta5$Original)][!(meta5$Original[grep("-w", meta5$Original)] %in% rownames(incidence))]


meta6 <- merge(meta5, treat, by.x = "Original", by.y = "Visit",
               all.x = TRUE, all.y = FALSE)

A <- list("RNAseq" = t(rna2), "Micro" = t(OTUs2), "Meta" = meta5)
A[1:2] <- clean_unvariable(A[1:2]) # Just the numeric ones
saveRDS(A, "data/RGCCA_data.RDS")
