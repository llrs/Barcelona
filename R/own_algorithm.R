library("experDesign")
pheno <- readRDS("samples2sequence.RDS")

keep <- c("IBD", "SEX", "Study", "Exact_location", "ID", "Time", "HSCT_responder")
omit <- colnames(pheno)[!colnames(pheno) %in% keep]
omit_wo <- omit[-1]
choose(nrow(pheno), 96)
i <- design(pheno, size.batch = 96, omit = omit, 100000)
saveRDS(i, file = "subsets.RDS")
out <- inspect(i, pheno, omit_wo)
out <- out[order(out$batch), ]
# o <- table(out$ID, out$batch)
sapply(keep, distribution, report = out)
o <- table(out$HSCT_responder, out$Study, out$batch, useNA = "ifany")
