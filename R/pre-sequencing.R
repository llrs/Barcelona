#' # Load data

trim <- readRDS("../TRIM/meta_r_full.RDS")
bcn <- readRDS("../BCN/meta_BCN_full.RDS")

#' # Clean the data
#'
#' Removes the unwanted columns for plate organization
keep <- colnames(bcn)[-c(
  2, 3, 4, 5, 6, 8, 15, 16, 17, 22, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35,
  36, 37, 38, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56,
  57, 58)]
bcn <- bcn[, keep]

keep <- colnames(trim)[-c(4, 5, 6, 8, 13, 15, 16, 17, 23)]
trim <- trim[, keep]

#' # Match column names
#'
#' This piece of code should match the content and column names of the two data
#' frames, in order to be able to merge them accordingly
colnames(bcn) <- gsub("gender", "SEX", colnames(bcn), ignore.case = TRUE)
trim$Exact_location <- tolower(gsub(" COLON", "", trim$Exact_location))
trim$ID <- ifelse(grepl("^C", as.character(trim$ID)),
                  as.character(trim$ID),
                  paste0("T-", as.character(trim$ID)))
trim$Study <- "HSCT"
bcn$Study <- "BCN"
pos <- grep("092", as.character(bcn$Pacient_id))
# Looking at the database and filling some missing values
bcn[pos, "Pacient_id"] <- "92"
bcn[pos, "SEX"] <- "female"
bcn[pos, "biopsied_segment"] <- "sigmoid"

bcn[bcn$Pacient_id == "C13", "SEX"] <- "male"
bcn[bcn$Pacient_id == "138", "SEX"] <- "female"
bcn[bcn$Sample_id == "138-w046", "biopsied_segment"] <- "sigmoid"
bcn[bcn$Pacient_id == "145", "SEX"] <- "female"
bcn[bcn$Sample_id == "145-w046", "biopsied_segment"] <- "descending"

bcn$Pacient_id <- gsub("092", "92", bcn$Pacient_id)
bcn$Pacient_id <- ifelse(grepl("^C", as.character(bcn$Pacient_id)),
                         gsub("-.*", "",as.character(bcn$Pacient_id)),
                         paste0("B-", as.character(bcn$Pacient_id)))
levels(bcn$IBD) <- c("CD", "CONTROL", "UC")
bcn <- droplevels(bcn)
trim <- droplevels(trim)

bcn$week <- as.character(bcn$week)
trim$Time <- gsub(".*-(.*)-TT?R-.*", "\\1", trim$`Sample Name_RNA`)
trim$Time[trim$IBD == "CONTROL"] <- NA

common <- c("Treatment", "IBD", "SEX", "Study")
pheno <- merge(trim, bcn, by.x = c("Sample Name_RNA", common, "Exact_location",
                                   "AgeDiag", "ID", "Time"),
               by.y = c("Sample_id", common, "biopsied_segment",
                        "Diagnostic age", "Pacient_id", "week") , all = TRUE)

pheno$SEX <- as.factor(tolower(pheno$SEX))
pheno <- pheno[!duplicated(pheno$`Sample Name_RNA`), ]
# Should be 380 samples
stopifnot(nrow(pheno) == 380)
table(pheno$IBD, pheno$Study)

# Add the phenotype
bcn.responders <- read.csv("bcn_responders.csv")

bcn.responders <- bcn.responders[match(pheno$`Sample Name_RNA`, bcn.responders$muestra),]
bcn.responders$Responder <- as.character(bcn.responders$Responder)
bcn.responders$Responder[bcn.responders$Responder %in% "Responder"] <- "YES"
bcn.responders$Responder[bcn.responders$Responder %in% "Non-responder"] <- "NO"

bcn.subset <- bcn.responders[, c("muestra", "Responder")]

pheno$HSCT_responder[!is.na(bcn.subset$Responder)] <- bcn.subset$Responder[!is.na(bcn.subset$Responder)]





saveRDS(pheno, "samples2sequence.RDS")
