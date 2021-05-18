library("integration")
library("RGCCA")
library("dplyr")
library("inteRmodel")
library("ggplot2")

A <- readRDS("data/RGCCA_ileum_data.RDS")
meta <- A$Meta

stopifnot(all(A$Meta$Original == rownames(A$RNAseq)))
testing <- function(x, ...) {
  tryCatch({
    result.sgcca <- RGCCA::sgcca(C = x,
                                 scheme = "centroid",
                                 verbose = FALSE,
                                 scale = FALSE,
                                 ...)
    analyze(result.sgcca)}
    , error = function(x){NA})
}

meta$`Phenotype CD` <- tolower(meta$`Phenotype CD`)
meta$`Location` <- tolower(meta$`Location`)

Localization <- model_RGCCA(meta, c("Exact_location")) # With SESCD local it increase the AVE_inner
Time <- model_RGCCA(meta, c("AgeDiag", "Age"))
Demographics <- model_RGCCA(meta, c("ID","SEX"))
Time$AgeDiag[is.na(Time$AgeDiag)] <- 0 # Time has NA values
A2 <- A[1:2]
A2$Demographics <- Demographics
A2$Location <- Localization
A2$Time <- Time

A2 <- clean_unvariable(A2)
A2 <- A2[-4] # Location don't have any variable as all are from the ileum
saveRDS(A2, "output/model3_BCN_ileum.RDS")

# shrinkage <- vapply(A2[1:2], tau.estimate, numeric(1L)) # 0.11503779803812 0.318145965316924
shrinkage <- rep(1, 4)
names(shrinkage) <- names(A2)
shrinkage[1:2] <- c(RNAseq = 0.360116796197683, Micro = 0.543260440351768)
shrinkage[3:4] <- 1
names(shrinkage) <- names(A2)
Ab <- lapply(A2, function(x) scale2(x, bias = TRUE)/sqrt(NCOL(x)))



# The design of model 3
C <- matrix(
  0, ncol = length(Ab), nrow = length(Ab),
  dimnames = list(names(Ab), names(Ab))
)

designs <- weight_design(weights = 3, size = length(Ab))
keep <- vapply(designs, correct, logical(1L))
designs <- designs[keep]

# Subset the designs
# set.seed(46726279)
# s <- sample(designs, size = min(length(designs)*.1, 1000))
out <- sapply(designs, testing, A = Ab, c1 = shrinkage, USE.NAMES = FALSE)
out2 <- as.data.frame(t(out))
saveRDS(out2, "output/sample_model3_ileum_boot.RDS")

out1 <- readRDS("output/sample_model3_ileum_boot.RDS")
ggplot(out1, aes(AVE_inner, AVE_outer)) +
  geom_point()


out <- out1
ggplot(out, aes(AVE_inner, AVE_outer, color = cc1)) +
  geom_point()
best3 <- out[out$AVE_inner == max(out$AVE_inner), grep("var", colnames(out))]
best3 <- symm(designs[[1]], best3)
colnames(best3) <- names(Ab)
rownames(best3) <- names(Ab)

model3_best <- sgcca(A = Ab, c1 = shrinkage, C = best3, ncomp = rep(2, 4), scheme = "centroid")
model3_best <- improve.sgcca(model3_best, names(Ab))
saveRDS(model3_best, "model3_ileum_best.RDS")
