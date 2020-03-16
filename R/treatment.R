library("RGCCA")
library("ggplot2")
library("integration")
library("inteRmodel")

A <- readRDS("data/RGCCA_data.RDS")
A2 <- A[1:2]
meta <- A$Meta


Localization <- model_RGCCA(meta, c("Exact_location")) # With SESCD local it increase the AVE_inner
Treatment <- model_RGCCA(meta, "treatment")
Demographics <- model_RGCCA(meta, c("ID","SEX"))


# Basic ####
C <- matrix(
  0, ncol = length(A), nrow = length(A),
  dimnames = list(names(A), names(A))
)
model0 <- subSymm(C, "Micro", "Meta", 1)
model0 <- subSymm(model0, "Meta", "RNAseq", 1)

model0i <- subSymm(model0, 1, 1, 1)

shrinkage <- c(RNAseq = 0.142456045648404, Micro = 0.595614967666157, 1, 1, 1)

B <- A[1:2]
B$Meta <- model_RGCCA(meta, c("Exact_location", "treatment", "ID", "SEX"))

out <- search_model(A = B, bias = TRUE, verbose = FALSE,
                    ncomp = rep(1, length(B)),
                    C = model0, c1 = shrinkage[1:3], scheme = "factorial",
                    scale = TRUE, nWeights = 11)
columns <- grep("var", colnames(out))
model <- symm(C, out[which.max(out$AVE_inner), columns])
# We then look for a variation of the weights of this model
out2 <- iterate_model(A = B, C = model, c1 =shrinkage[1:3], scheme = "factorial",
                     scale = FALSE, verbose = FALSE,
                     ncomp = rep(1, length(A)),
                     bias = TRUE)
