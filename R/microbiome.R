

library("metagenomeSeq")

# Create the object
MR <- newMRexperiment(
  otus_table_i[, ordSamples],
  phenoData = AnnotatedDataFrame(meta_i[ordSamples, ]),
  featureData = AnnotatedDataFrame(as.data.frame(otus_tax_i))
)
filterData(MR, present = 10, depth = 1000)
p <- cumNormStatFast(MR)
MR <- cumNorm(MR, p = p)
assayData(MR)$relative <- assayData(MR)$counts / rowSums(assayData(MR)$counts)
assayData(MR)$prevalence <- as.matrix(assayData(MR)$counts != 0)

rareFeatures <- which(rowSums(MRcounts(MR) > 0) < 10)
MRtrim <- MR[-rareFeatures, ]
MRp <- cumNormStat(MRtrim, pFlag = TRUE, main = "Trimmed lung data")
MRtrim <- cumNorm(MRtrim, p = MRp)
assayData(MRtrim)



normFactor <- normFactors(MRtrim)
assayData(MRtrim)$relative <- assayData(MRtrim)$counts / rowSums(assayData(MRtrim)$counts)
assayData(MRtrim)$prevalence <- as.matrix(assayData(MRtrim)$counts != 0)

normFactor <- log2(normFactor / median(normFactor) + 1)
settings <- zigControl(maxit = 10, verbose = TRUE)
mod <- model.matrix(~ 0 + Endoscopic_Activity, data = pData(MRtrim))
mod <- cbind(mod, ID = pData(MRtrim)$ID, Time = as.numeric(as.factor(pData(MRtrim)$Time)))
fit <- fitZig(
  obj = MRtrim, mod = mod, useCSSoffset = FALSE,
  control = settings
)

eb <- fit$eb
head(cbind(eb$F, eb$F.p.value))
