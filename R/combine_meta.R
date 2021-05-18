library("readr")
library("dplyr")
library("data.table")
library("readxl")
library("ggplot2")
library("stringr")
library("purrr")
library("tidyr")

# From the QC step
meta <- readRDS("output/info_samples.RDS")

# Remove duplicate samples
replicates <- table(meta$Original)
replicate_samples <- meta[meta$Original %in% names(replicates[replicates > 1]), ]

correct_bcn <- function(x) {
  if (length(x) > 1) {
    a <- str_pad(x[1], width = 3, pad = "0")
    b <- str_pad(x[2], width = 3, pad = "0")
    x <- paste(a, b, sep = "-w")
  }
  x
}


# Meta ####
db <- data.table::fread(
  "data/db_biopsies_bcn_seq16S_noTRIM.txt", sep = "\t",
  stringsAsFactors = FALSE)

meta3 <- merge(meta, db,
               by.x = "Original", by.y = "Sample_Code", all = TRUE) %>%
  droplevels() %>%
  mutate(IBD = NA) %>%
  mutate_if(is.factor, as.character)

both_equal <- which(meta3$IBD.x == meta3$IBD.y)
meta3$IBD[both_equal] <- meta3$IBD.x[both_equal]
x_empty <- is.na(meta3$IBD.x) & !is.na(meta3$IBD.y)
meta3$IBD[x_empty] <- meta3$IBD.y[x_empty]
y_empty <- is.na(meta3$IBD.y) & !is.na(meta3$IBD.x)
meta3$IBD[y_empty] <- meta3$IBD.x[y_empty]

contr_original <- grepl("^C", meta3$Original)
meta3[contr_original, "Patient_ID"] <- gsub("-.*", "", meta3$Original)[contr_original]

# From the copy of the access database I have 04/01/2019
dates <- read_xlsx("data/fechas.xlsx",
                   col_types = c("text", "text", "guess", "guess"), na = "")
dates$Date <- as.Date(dates$Date)
dates$`Date of diagnosis` <- as.Date(dates$`Date of diagnosis`)
colnames(dates) <- c("Visit", "ID", "DATE_SAMPLE", "Date_diagnostic")
meta4 <- merge(meta3, dates, by.x = "Original", by.y = "Visit", all.x = TRUE, all.y = FALSE)

# Arrange the IDs
contr4_original <- grepl("^C", meta4$Original)
meta4$ID[contr4_original] <- gsub("^(C[0-9]+)-.+", "\\1", meta4$Original[contr4_original])
meta4$ID[grepl("^001", meta4$Original) & is.na(meta4$ID)] <- "001" # Fill an NA

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
p17 <- meta4[meta4$ID == "017", c("Original", "ID", "DATE_SAMPLE", "Birth_date", "Age", "Patient_ID")]
meta4$Birth_date[meta4$ID == "017"] <- p17$Birth_date[1]
meta4$Patient_ID[meta4$ID == "017"] <- p17$Patient_ID[1]


ggplot(meta4[grep("^[0-9]", meta4$ID), ]) +
  geom_point(aes(as.numeric(ID), DATE_SAMPLE, color = factor(substr(Original, 6, 8)))) +
  labs(x = "ID", y = "Date sample", color = "Sample id (weeks)")
# See the plot!!
# The patient ID 113 has the week 46 before the week 000??
# Also 2 missing DATE_SAMPLES warning message

# Using the other database we can complete and correct the record:
meta4$DATE_SAMPLE[meta4$Original == "001-w000"] <- as.Date("11/14/2012", "%m/%d/%Y")
meta4$DATE_SAMPLE[meta4$Original == "111-w038"] <- as.Date("06/07/2016", "%m/%d/%Y")
meta4$DATE_SAMPLE[meta4$Original == "113-w046"] <- as.Date("08/09/2016", "%m/%d/%Y")

ggplot(meta4[grep("^[0-9]", meta4$ID), ]) +
  geom_point(aes(DATE_SAMPLE, as.numeric(ID), color = factor(substr(Original, 6, 8)))) +
  labs(y = "ID", x = "Date sample", color = "Phase (weeks)",
       title = "Processed samples") +
  theme_minimal()

# Calculate years since diagnostic and related problems
diagTime <- as.Date(meta4$DATE_SAMPLE, "%Y-%m-%d") - meta4$Date_diagnostic
diagTime <- as.numeric(diagTime/365.25)
diagTime[is.na(diagTime)] <- 0 # Replace by date of the sample at least
AgeDiag <- as.numeric(meta4$Date_diagnostic - as.Date(meta4$Birth_date,
                                                      "%d/%m/%Y"))/365.25

# Calculate the age because they are off up to 5 years up
# hist(meta4$Age - as.numeric(meta4$DATE_SAMPLE - meta4$Birth_date)/365.25)
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
  purrr::map(correct_bcn) %>%
  unlist()
db2$IBD[db2$IBD == "ctrl"] <- "CONTROL"

# Based that on the resequenced is after the original
db3 <- db2[-c(which(duplicated(db2$Sample_id) == TRUE) - 1), ]
meta5 <- merge(meta4, db3,
               by.x = "Original",
               by.y = "Sample_id",
               all = TRUE, suffixes = c(".z", ".w"))

meta5$IBD <- NA
both_equal <- which(meta5$IBD.z == meta5$IBD.w)
meta5$IBD[both_equal] <- meta5$IBD.z[both_equal]
x_empty <- is.na(meta5$IBD.z) & !is.na(meta5$IBD.w)
meta5$IBD[x_empty] <- meta5$IBD.w[x_empty]
y_empty <- is.na(meta5$IBD.w) & !is.na(meta5$IBD.z)
meta5$IBD[y_empty] <- meta5$IBD.z[y_empty]

meta5$sample_date <- as.Date(meta5$sample_date, "%m/%d/%Y")
meta5$Time[meta5$Name == "017-w014"] <- "14"

# Different sample date?? The right one is DATE_SAMPLE
# meta5[meta5$sample_date != meta5$DATE_SAMPLE & !is.na(meta5$sample_date),
#       c("sample_date", "DATE_SAMPLE", "Original")]
meta5 <- meta5[, -grep("sample_date", colnames(meta5))]

duplicates <- group_by(meta5, NHC) %>%
  summarise(diff = n_distinct(Patient_ID))
dupli <- duplicates$NHC[duplicates$diff > 1]
dupli <- !is.na(dupli)

# How do we codify it ?? The one with less samples loses its Patient_ID
meta5[meta5$NHC %in% dupli, "Patient_ID"] <- c("17", "17", "17", "86", "92",
                                               "86", "86", "86", "92", "92")

# Merging data of exact localization and remove gender
consulta <- read_xlsx("data/20200306_consulta_BBDD.xlsx")
meta6 <- merge(meta5, consulta[, -4], all.x = TRUE, all.y = FALSE,
               by.x = "Original", by.y = "muestra")



# Check which samples have RNA ####
conn <- gzfile("data/TNF.all.samples.original.counts.tsv.gz") # TODO See if this is a good choice
rna <- read.table(conn, sep = "\t", check.names = FALSE, nrows = 2)

samples <- colnames(rna)
samples <- gsub(" reseq$", "", samples)
samples[grep("[Ww]", samples)] <- tolower(samples[grep("[Ww]", samples)])
samples <- samples %>%
  str_split("-w") %>% # Ready for BCN
  map(correct_bcn) %>%
  unlist() %>%
  gsub("-T-TR-", "-T-DM-", .)

# Check which samples have 16S ####
tab <- read.delim("data/Partek_Michigan3_Kraken_Classified_genus.tsv",
                  check.names = FALSE, nrows = 2)
dna_samples <- colnames(tab)[-1]
dna_samples <- gsub("_S.+", "", rna_samples)



inside <- samples %in% meta6$Original
stopifnot(all(inside))
meta6$Transcriptome <- ifelse(meta6$Original %in% samples, "Yes", "No")
meta6$Microbiome <- ifelse(meta6$Name %in% dna_samples, "Yes", "No")

write.csv(meta6[(meta6$Transcriptome == "Yes" | meta6$Microbiome == "Yes"), ],
          "output/pheno_data.csv")
