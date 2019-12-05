library("integration")
library("RGCCA")
library("inteRmodel")
A <- readRDS("data/RGCCA_data.RDS")
meta <- A$Meta

shrinkage[3] <- 1


A$Meta <- model_RGCCA(A$Meta, c("ID", "AgeDiag", "diagTime", "Exact_location", "SEX"))
Ab <- lapply(A, function(x) scale2(x, bias = TRUE)/sqrt(NCOL(x)))
shrinkage[1:2] <- sapply(Ab[1:2], tau.estimate) # 0.11503779803812 0.318145965316924

out_model <- search_model(A = Ab, c1 = shrinkage, scheme = "centroid",
                          scale = FALSE, verbose = FALSE,
                          bias = TRUE, nWeights = 11)

saveRDS(out_model, "model2_optimization.RDS")
out2 <- readRDS("model2_optimization.RDS")
model2b <- symm(C, out2[out2$AVE_inner == max(out2$AVE_inner),
                        grep("^var", colnames(out2))])

model2b_sgcca <- sgcca(A = Ab, c1 = shrinkage, scheme = "centroid", C = model2b,
                       verbose = FALSE, scale = FALSE, ncomp = rep(2, length(Ab)))

model2b_sgcca <- improve.sgcca(model2b_sgcca, names(A))
saveRDS(model2b_sgcca, "model2b_sgcca.RDS")
plot(model2b_sgcca$Y$RNAseq[, 1], model2b_sgcca$Y$Micro[, 1], col = meta$Exact_location)







A3 <- model_RGCCA(meta, c("ID", "AgeDiag", "diagTime", "Exact_location", "SEX", "IBD"))
Ab[[3]] <- scale2(A3, bias = TRUE)/sqrt(NCOL(A3))

# Estimated time of 1 hours
out_model <- search_model(A = Ab, c1 = shrinkage, scheme = "centroid",
                          scale = FALSE, verbose = FALSE,
                          bias = TRUE, nWeights = 11)

saveRDS(out_model, "model2_optimization_IBD.RDS")
out2 <- readRDS("model2_optimization_IBD.RDS")
model2b2 <- symm(C, out2[out2$AVE_inner == max(out2$AVE_inner),
                        grep("^var", colnames(out2))])

model2b2_sgcca <- sgcca(A = Ab, c1 = shrinkage, scheme = "centroid", C = model2b,
                       verbose = FALSE, scale = FALSE, ncomp = rep(2, length(Ab)))
model2b2_sgcca <- improve.sgcca(model2b2_sgcca, names(A))
saveRDS(model2b2_sgcca, "model2b2_sgcca.RDS")

plot(model2b2_sgcca$Y$RNAseq[, 1], model2b2_sgcca$Y$Micro[, 1], col = meta$IBD)
plot(model2b2_sgcca$Y$RNAseq[, 1], model2b2_sgcca$Y$Micro[, 1], col = as.factor(ifelse(meta$Exact_location == "ileum", "ileum", "colon")))
