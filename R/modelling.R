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
names(models0) <- c("model0", "model0i")
out <- lapply(models0, use, A = Ab2, c1 = shrinkage[1:2])
saveRDS(out, "models0.RDS")
samples <- sapply(out[[1]]$Y, function(x) {
  x[, 1]
})
comm <- ggplot(as.data.frame(samples), aes(RNAseq, Micro)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  ggtitle("All samples at all times ") +
  xlab("RNAseq (component 1)") +
  ylab("16S (component 1)") +
  theme(plot.title = element_text(hjust = 0.5))
comm +
  geom_text(aes(color = A$Meta$IBD.x, label = A$Meta$Original)) +
  guides(col = guide_legend(title = "Patient"))
# Plot not interesting low AVE and not separating by disease or controls

samples <- sapply(out[[2]]$Y, function(x) {
  x[, 1]
})
comm <- ggplot(as.data.frame(samples), aes(RNAseq, Micro)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  ggtitle("All samples at all times ") +
  xlab("RNAseq (component 1)") +
  ylab("16S (component 1)") +
  theme(plot.title = element_text(hjust = 0.5))
comm +
  geom_text(aes(color = A$Meta$IBD.x, label = A$Meta$Original)) +
  guides(col = guide_legend(title = "Patient"))

# TODO: Merge with the original data.frame with information about the disease.
models2 <- list(model1, model1i, model2, model2i)
Ab <- lapply(A, function(x) scale2(x, bias = TRUE)/sqrt(NCOL(x)))
out <- lapply(models2, use, A = Ab2, c1 = shrinkage[1:2])
