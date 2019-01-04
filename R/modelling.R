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
  0, ncol = length(A2), nrow = length(A2),
  dimnames = list(names(A2), names(A2))
)
model1 <- subSymm(C, "Micro", "meta", 1)
model1 <- subSymm(model1, "RNaseq", "meta", 1)
model1i <- subSymm(model1, 1, 1, 1)
model2 <- subSymm(model1, "RNaseq", "Micro", 1)
model2i <- subSymm(model2, 1, 1, 1)



# We cannnot comput eht tau.estimate for A[[1]]
# (shrinkage <- sapply(A, tau.estimate))
shrinkage <- c(0.249488046688595, 0) # We guess a 0.1
shrinkage[2] <- tau.estimate(A[[2]])
(min_shrinkage <- sapply(A, function(x) {
  1 / sqrt(ncol(x))
}))
# # Don't let the shrinkage go below the threshold  allowed
(shrinkage <- ifelse(shrinkage < min_shrinkage, min_shrinkage, shrinkage))
# shrinkage <- rep(1, length(A))

ncomp <- rep(2, length(A))
sgcca.centroid <- sgcca(
  A, C = model0i, c1 = shrinkage,
  ncomp = ncomp,
  scheme = "centroid",
  scale = TRUE,
  verbose = FALSE
)
sgcca.centroid <- improve.sgcca(sgcca.centroid, names(A))
