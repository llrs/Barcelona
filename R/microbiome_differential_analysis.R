library("phyloseq")
library("metagenomeSeq")

# OTUs ####
{
  # * Read ####
  tab <- read.delim("data/Partek_Michigan3_Kraken_Classified_family.tsv", check.names = FALSE)
  # tab <- read.delim("data/20200529_Partek_Michigan3_Kraken_Classified_phylum.txt", check.names = FALSE)
  colnames(tab) <- gsub("_S.*", "", colnames(tab)) # Remove trailing numbers
  counts <- tab[, -1]
  microorganism <- tab[, 1, FALSE]

  # From the QC step
  meta <- readRDS("data_out/refined_meta.RDS")
  otus <- tab[, colnames(tab) %in% meta$Name]
  colnames(otus) <- meta$Original[match(colnames(otus), meta$Name)]
  {
    # Check that they match
    otus_new <- meta$Original[match(colnames(otus), meta$Name)]
    names(otus_new) <- colnames(otus)
    otus_new[names(otus_new) != otus_new]
  }

  # Reorder samples to match!
  otus <- otus[, match(meta$Original, colnames(otus))]
  rownames(meta) <- meta$Original
  stopifnot(all(colnames(otus) == meta$Original))
  # microorganism <- read.csv("data/family.csv", row.names = 1)

  meta$ileum <- ifelse(meta$Exact_location == "ileum", "Ileum", "Colon")
  # The missing values of Exact location
  meta$ileum[is.na(meta$ileum )] <- "Colon"
  meta$SEX[is.na(meta$SEX) | meta$SEX == ""] <- "female"
  meta$Time[is.na(meta$Time)] <- "C"

  otus <- as.matrix(otus)
  meta$IBD <- as.character(meta$IBD)
  meta$IBD[is.na(meta$IBD)] <- "C"
  rownames(otus) <- paste0("sp", seq_len(nrow(otus)))
  phyloseq <- phyloseq(otu_table(otus, taxa_are_rows = TRUE),
                       sample_data(meta),
                       tax_table(as.matrix(microorganism)))
}
{
  # * Comparisons ####
  MR <- phyloseq_to_metagenomeSeq(phyloseq) # For testing and comparing data
  filterData(MR, present = 10, depth = 1000)
  p <- cumNormStatFast(MR)
  MR <- cumNorm(MR, p = p)
  assayData(MR)$relative <- assayData(MR)$counts / rowSums(assayData(MR)$counts)
  assayData(MR)$prevalence <- as.matrix(assayData(MR)$counts != 0)

  rareFeatures <- which(rowSums(MRcounts(MR) > 0) < 10)
  MRtrim <- MR[-rareFeatures, ]
  MRp <- cumNormStat(MRtrim, pFlag = TRUE, main = "Trimmed data")
  MRtrim <- cumNorm(MRtrim, p = MRp)

  normFactor <- normFactors(MRtrim)
  assayData(MRtrim)$relative <- assayData(MRtrim)$counts / rowSums(assayData(MRtrim)$counts)
  assayData(MRtrim)$prevalence <- as.matrix(assayData(MRtrim)$counts != 0)

  normFactor <- log2(normFactor / median(normFactor) + 1)
  settings <- zigControl(maxit = 10, verbose = FALSE)
  mod <- model.matrix(~ 0 + IBD + SEX, data = pData(MRtrim))
  Time <- as.numeric(pData(MRtrim)$Time)
  Time[is.na(Time)] <- 0
  Ileum <- ifelse(pData(MRtrim)$Exact_location == "ileum", 1, 0)
  mod <- cbind(mod, Time = Time, Ileum = Ileum)
  fit <- fitZig(
    obj = MRtrim, mod = mod, useCSSoffset = FALSE,
    control = settings,
    block =  pData(MRtrim)$ID
  )
  zigFit <- slot(fit, "fit")
  finalMod <- slot(fit, "fit")$design
  contrast.matrix <- makeContrasts(CD = IBDCD - IBDC,
                                   UC = IBDUC - IBDC,
                                   Ileum_vs_Colon = Ileum,
                                   Male_vs_femal = SEXmale,
                                   levels = finalMod)
  fit2 <- contrasts.fit(zigFit, contrast.matrix)
  fit2 <- eBayes(fit2)
  dt <- decideTests(fit2)
  summary(dt)
}

# ASV ####
{
  # * Read ####
  ASV <- readRDS("data_out/refined_ASV.RDS")

  # From the QC step
  meta <- readRDS("data_out/refined_meta.RDS")
  ASV <- ASV[, colnames(ASV) %in% meta$Name]
  colnames(ASV) <- meta$Original[match(colnames(ASV), meta$Name)]
  {
    # Check that they match
    ASV_new <- meta$Original[match(colnames(ASV), meta$Name)]
    names(ASV_new) <- colnames(ASV)
    ASV_new[names(ASV_new) != ASV_new]
  }

  # Reorder samples to match!
  ASV <- ASV[, match(meta$Original, colnames(ASV))]
  rownames(meta) <- meta$Original
  stopifnot(all(colnames(ASV) == meta$Original))
  # microorganism <- read.csv("data/family.csv", row.names = 1)

  meta$ileum <- ifelse(meta$Exact_location == "ileum", "Ileum", "Colon")
  # The missing values of Exact location
  meta$ileum[is.na(meta$ileum )] <- "Colon"
  meta$SEX[is.na(meta$SEX) | meta$SEX == ""] <- "female"
  meta$Time[is.na(meta$Time)] <- "C"
  meta$IBD <- as.character(meta$IBD)
  meta$IBD[is.na(meta$IBD)] <- "C"

  ASV <- as.matrix(ASV)
  rownames(ASV) <- paste0("sp", seq_len(nrow(ASV)))
  phyloseq_ASV <- phyloseq(otu_table(ASV, taxa_are_rows = TRUE),
                       sample_data(meta))
}

{
  # * Comparisons ####
  MR <- phyloseq_to_metagenomeSeq(phyloseq_ASV) # For testing and comparing data
  filterData(MR, present = 10, depth = 1000)
  p <- cumNormStatFast(MR)
  MR <- cumNorm(MR, p = p)
  assayData(MR)$relative <- assayData(MR)$counts / rowSums(assayData(MR)$counts)
  assayData(MR)$prevalence <- as.matrix(assayData(MR)$counts != 0)

  rareFeatures <- which(rowSums(MRcounts(MR) > 0) < 10)
  MRtrim <- MR[-rareFeatures, ]
  MRp <- cumNormStat(MRtrim, pFlag = TRUE, main = "Trimmed data")
  MRtrim <- cumNorm(MRtrim, p = MRp)

  normFactor <- normFactors(MRtrim)
  assayData(MRtrim)$relative <- assayData(MRtrim)$counts / rowSums(assayData(MRtrim)$counts)
  assayData(MRtrim)$prevalence <- as.matrix(assayData(MRtrim)$counts != 0)

  normFactor <- log2(normFactor / median(normFactor) + 1)
  settings <- zigControl(maxit = 10, verbose = FALSE)
  mod <- model.matrix(~ 0 + IBD + SEX, data = pData(MRtrim))
  Time <- as.numeric(pData(MRtrim)$Time)
  Time[is.na(Time)] <- 0
  Ileum <- ifelse(pData(MRtrim)$Exact_location == "ileum", 1, 0)
  mod <- cbind(mod, Time = Time, Ileum = Ileum)
  fit <- fitZig(
    obj = MRtrim, mod = mod, useCSSoffset = FALSE,
    control = settings,
    block =  pData(MRtrim)$ID
  )
  zigFit <- slot(fit, "fit")
  finalMod <- slot(fit, "fit")$design
  contrast.matrix <- makeContrasts(CD = IBDCD - IBDC,
                                   UC = IBDUC - IBDC,
                                   Ileum_vs_Colon = Ileum,
                                   Male_vs_femal = SEXmale,
                                   levels = finalMod)
  fit2 <- contrasts.fit(zigFit, contrast.matrix)
  fit2 <- eBayes(fit2)
  dt2 <- decideTests(fit2)
  summary(dt2)
}

# In conclusion ASV is filtered more, but we can see more differential abundance
