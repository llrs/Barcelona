library("integration")
library("RGCCA2")

A <- readRDS("data/RGCCA_data.RDS")
A2 <- A[1:2]

# The design
C <- matrix(
  0, ncol = length(A2), nrow = length(A2),
  dimnames = list(names(A2), names(A2))
)
model0 <- subSymm(C, "Micro", "RNAseq", 1)
model0i <- subSymm(model0, 1, 1, 1)

C <- matrix(
  0, ncol = length(A), nrow = length(A),
  dimnames = list(names(A), names(A))
)
model1 <- subSymm(C, "Micro", "Meta", 1)
model1 <- subSymm(model1, "RNAseq", "Meta", 1)
model1i <- subSymm(model1, 1, 1, 1)
model2 <- subSymm(model1, "RNAseq", "Micro", 1)
model2i <- subSymm(model2, 1, 1, 1)



# We cannnot comput eht tau.estimate for A[[1]]
# (shrinkage <- sapply(A, tau.estimate))
shrinkage <- sapply(A2, tau.estimate)
(min_shrinkage <- sapply(A, function(x) {
  1 / sqrt(ncol(x))
}))
# # Don't let the shrinkage go below the threshold  allowed
(shrinkage <- ifelse(shrinkage < min_shrinkage, min_shrinkage, shrinkage))

Ab2 <- lapply(A2, function(x) scale2(x, bias = TRUE)/sqrt(NCOL(x)))
use <- function(...){
  sgcca.centroid <- sgcca(
    scheme = "centroid",
    scale = FALSE,
    verbose = FALSE,
    ...
  )
  args <- list(...)
  improve.sgcca(sgcca.centroid, names(args$A))
}


models0 <- list(model0, model0i)
out <- lapply(models0, use, A = Ab2, c1 = shrinkage[1:2])


# TODO: Merge with the original data.frame with information about the disease.
models2 <- list(model1, model1i, model2, model2i)
Ab <- lapply(A, function(x) scale2(x, bias = TRUE)/sqrt(NCOL(x)))
out <- lapply(models2, use, A = Ab2, c1 = shrinkage[1:2])
