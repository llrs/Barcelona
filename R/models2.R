library("integration")
library("RGCCA2")

A <- readRDS("data/RGCCA_data.RDS")

shrinkage <- sapply(A[1:2], tau.estimate)
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
