library("integration")
library("RGCCA")
library("inteRmodel")
library("BiocParallel")

mcp <- MulticoreParam(workers = 7, progressbar = TRUE)
A <- readRDS("data/RGCCA_data_wo_out.RDS")
meta <- A$Meta

shrinkage <- rep(1, 3)

A$Meta <- model_RGCCA(A$Meta, c("ID", "AgeDiag", "diagTime", "Exact_location", "SEX"))
Ab <- lapply(A, function(x) scale2(x, bias = TRUE)/sqrt(NCOL(x)))
# shrinkage[1:2] <- vapply(Ab[1:2], tau.estimate, numeric(1L)) # c(0.322297910454825, 0.959997006295785) # With outliers

shrinkage[1:2] <- c(0.322020648273615, 0.866155549496009)


out_model <- search_model(A = Ab, c1 = shrinkage, scheme = "centroid",
                          scale = FALSE, verbose = FALSE,
                          bias = TRUE, nWeights = 11,
                          BPPARAM = mcp)

saveRDS(out_model, "output/model2_optimization_b.RDS")
out2 <- readRDS("output/model2_optimization_b.RDS")
C <- matrix(0, nrow = 3, ncol = 3)
model2b <- symm(C, out2[out2$AVE_inner == max(out2$AVE_inner),
                        grep("^var", colnames(out2))])

model2b_sgcca <- sgcca(A = Ab, c1 = shrinkage, scheme = "centroid", C = model2b,
                       verbose = FALSE, scale = FALSE, ncomp = rep(2, length(Ab)))

model2b_sgcca <- improve.sgcca(model2b_sgcca, names(A))
saveRDS(model2b_sgcca, "output/model2b_sgcca_b.RDS")
model2b_sgcca <- readRDS("output/model2b_sgcca_b.RDS")
plot(model2b_sgcca$Y$RNAseq[, 1], model2b_sgcca$Y$Micro[, 1], col = as.factor(meta$Exact_location))

A3 <- model_RGCCA(meta, c("ID", "AgeDiag", "diagTime", "Exact_location", "SEX", "IBD"))
Ab[[3]] <- scale2(A3, bias = TRUE)/sqrt(NCOL(A3))

C <- matrix(0, nrow = 3, ncol = 3)
subOut <- out2[out2$var13 == 1 & out2$var23 == 1 & out2$var12 == 0, ]
vars <- subOut[subOut$AVE_inner == max(subOut$AVE_inner),
               grep("^var", colnames(subOut))]
model1.1 <- symm(C, vars)

model1.1_sgcca <- sgcca(A = Ab, c1 = shrinkage, scheme = "centroid", C = model1.1,
                       verbose = FALSE, scale = FALSE, ncomp = rep(2, length(Ab)))

model1.1_sgcca <- improve.sgcca(model1.1_sgcca, names(A))
saveRDS(model1.1_sgcca, "output/model1.1_sgcca_b.RDS")

# Estimated time of 1 hours
out_model <- search_model(A = Ab, c1 = shrinkage, scheme = "centroid",
                          scale = FALSE, verbose = FALSE,
                          bias = TRUE, nWeights = 11)

saveRDS(out_model, "output/model2_optimization_IBD_b.RDS")
out2 <- readRDS("output/model2_optimization_IBD_b.RDS")
model2b2 <- symm(C, out2[out2$AVE_inner == max(out2$AVE_inner),
                        grep("^var", colnames(out2))])

model2b2_sgcca <- sgcca(A = Ab, c1 = shrinkage, scheme = "centroid", C = model2b,
                       verbose = FALSE, scale = FALSE, ncomp = rep(2, length(Ab)))
model2b2_sgcca <- improve.sgcca(model2b2_sgcca, names(A))
saveRDS(model2b2_sgcca, "output/model2b2_sgcca_b.RDS")

plot(model2b2_sgcca$Y$RNAseq[, 1], model2b2_sgcca$Y$Micro[, 1], col = as.factor(meta$IBD))
plot(model2b2_sgcca$Y$RNAseq[, 1], model2b2_sgcca$Y$Micro[, 1], col = as.factor(ifelse(meta$Exact_location == "ileum", "ileum", "colon")))
