library("integration")
library("RGCCA2")
A <- readRDS("data/RGCCA_data.RDS")
meta <- A$Meta

shrinkage <- sapply(A[1:2], tau.estimate) # 0.11503779803812 0.318145965316924
shrinkage[3] <- 1

C <- matrix(
  0, ncol = length(A), nrow = length(A),
  dimnames = list(names(A), names(A))
)

designs <- weight_design(weights = 11, size = length(A))
keep <- vapply(designs, correct, logical(1L))
designs <- designs[keep]

testing <- function(x, ...) {
  result.sgcca <- RGCCA2::sgcca(C = x,
                                scheme = "centroid",
                                verbose = FALSE,
                                scale = FALSE,
                                ...)
  analyze(result.sgcca)
}

A$Meta <- model_RGCCA(A$Meta, c("ID", "AgeDiag", "diagTime", "Exact_location", "SEX"))
Ab <- lapply(A, function(x) scale2(x, bias = TRUE)/sqrt(NCOL(x)))
# Estimated time of 8 hours
out <- sapply(designs, testing, A = Ab, c1 = shrinkage, USE.NAMES = FALSE)
out2 <- as.data.frame(t(out))
saveRDS(out2, "model2_optimization.RDS")
out2 <- readRDS("model2_optimization.RDS")
model2b <- symm(C, out2[out2$AVE_inner == max(out2$AVE_inner),
                        grep("^var", colnames(out2))])

model2b_sgcca <- sgcca(A = Ab, c1 = shrinkage, scheme = "centroid", C = model2b,
                       verbose = FALSE, scale = FALSE, ncomp = rep(2, length(Ab)))
model2b_sgcca <- improve.sgcca(model2b_sgcca, names(A))
saveRDS(model2b_sgcca, "model2b_sgcca.RDS")
plot(model2b_sgcca$Y$RNAseq[, 1], model2b_sgcca$Y$Micro[, 1], col = meta$Exact_location)







A$Meta <- model_RGCCA(meta, c("ID", "AgeDiag", "diagTime", "Exact_location", "SEX", "IBD"))
Ab <- lapply(A, function(x) scale2(x, bias = TRUE)/sqrt(NCOL(x)))
# Estimated time of 1 hours
out <- sapply(designs, testing, A = Ab, c1 = shrinkage, USE.NAMES = FALSE)
out2 <- as.data.frame(t(out))
saveRDS(out2, "model2_optimization_IBD.RDS")

model2b2 <- symm(C, out2[out2$AVE_inner == max(out2$AVE_inner), grep("^var", colnames(out2))])

model2b2_sgcca <- sgcca(A = Ab, c1 = shrinkage, scheme = "centroid", C = model2b2,
                       verbose = FALSE, scale = FALSE, ncomp = rep(2, length(Ab)))
model2b2_sgcca <- improve.sgcca(model2b2_sgcca, names(A))
saveRDS(model2b2_sgcca, "model2b2_sgcca.RDS")
plot(model2b2_sgcca$Y$RNAseq[, 1], model2b2_sgcca$Y$Micro[, 1], col = meta$IBD)
plot(model2b2_sgcca$Y$RNAseq[, 1], model2b2_sgcca$Y$Micro[, 1], col = as.factor(ifelse(meta$Exact_location == "ileum", "ileum", "colon")))
