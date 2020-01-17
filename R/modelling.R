library("integration")
library("RGCCA")
library("ggplot2")

A <- readRDS("data/RGCCA_data.RDS")
A2 <- A[1:2]

# Basic ####
C <- matrix(
  0, ncol = length(A2), nrow = length(A2),
  dimnames = list(names(A2), names(A2))
)
model0 <- subSymm(C, "Micro", "RNAseq", 1)
model0i <- subSymm(model0, 1, 1, 1)

# We cannnot comput the tau.estimate for A[[1]]
# (shrinkage <- sapply(A, tau.estimate))
shrinkage <- c(0.285693348851943, 0, 1) #Calculated from the server for the data derived from original data
shrinkage[2] <- tau.estimate(A2[[2]])
(min_shrinkage <- sapply(A, function(x) {
  1 / sqrt(ncol(x))
}))
# # Don't let the shrinkage go below the threshold  allowed
shrinkage[shrinkage < min_shrinkage] <- min_shrinkage[shrinkage < min_shrinkage]

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
  geom_text(aes(color = A$Meta$IBD, label = A$Meta$Original)) +
  guides(col = guide_legend(title = "Patient"))
comm +
  geom_text(aes(color = A$Meta$Exact_location, label = A$Meta$Original)) +
  guides(col = guide_legend(title = "Location"))
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
  geom_text(aes(color = A$Meta$IBD, label = A$Meta$Original)) +
  guides(col = guide_legend(title = "Patient"))

# With metadata ####
C <- matrix(
  0, ncol = length(A), nrow = length(A),
  dimnames = list(names(A), names(A))
)
model1 <- subSymm(C, "Micro", "Meta", 1)
model1 <- subSymm(model1, "RNAseq", "Meta", 1)
model1i <- subSymm(model1, 1, 1, 1)
model2 <- subSymm(model1, "RNAseq", "Micro", 1)
model2i <- subSymm(model2, 1, 1, 1)
model2b <- subSymm(model1, "RNAseq", "Meta", 0.1)
model2bi <- subSymm(model2b, 1, 1, 1)

models2 <- list(model1, model1i, model2, model2i, model2b, model2bi)
names(models2) <- c("model1", "model1i", "model2", "model2i", "model2b", "model2bi")
A2 <- A
A2$Meta <- model_RGCCA(A$Meta, c("ID", "AgeDiag", "diagTime", "Exact_location", "SEX"))
A2b <- lapply(A2, function(x) scale2(x, bias = TRUE)/sqrt(NCOL(x)))
shrinkage[seq(3, length(shrinkage))] <- 1
out <- lapply(models2, use, A = A2b, c1 = shrinkage)
saveRDS(out, "models2.RDS")

# Complex models ####
Localization <- model_RGCCA(A$Meta, c("Exact_location")) # With SESCD local it increase the AVE_inner
Time <- model_RGCCA(A$Meta, c("AgeDiag", "Age"))
Demographics <- model_RGCCA(A$Meta, c("ID","SEX"))

A3 <- A[1:2]
A3$Demographics <- Demographics
A3$Localization <- Localization
A3$Time <- Time
A3b <- lapply(A3, function(x) scale2(x, bias = TRUE)/sqrt(NCOL(x)))

C <- matrix(
  0, ncol = length(A3), nrow = length(A3),
  dimnames = list(names(A3), names(A3))
)
model3 <- subSymm(C, "Micro", "RNAseq", 1)
model3 <- subSymm(model3, "RNAseq", "Demographics", 1)
model3 <- subSymm(model3, "RNAseq", "Localization", 1)
model3 <- subSymm(model3, "RNAseq", "Time", 1)
model3 <- subSymm(model3, "Micro", "Demographics", 1)
model3 <- subSymm(model3, "Micro", "Localization", 1)
model3 <- subSymm(model3, "Micro", "Time", 1)
model3i <- subSymm(model3, 1, 1, 1)

model3b <- subSymm(C, "RNAseq", "Localization", 1)
model3b <- subSymm(model3b, "Demographics", "Micro", 1)
model3b <- subSymm(model3b, "Localization", "Micro", 0.5)
model3b <- subSymm(model3b, "Demographics", "Time", 1)
model3bi <- subSymm(model3b, 1, 1, 1)

models3 <- list(model3, model3i, model3b, model3bi)
shrinkage3 <- rep(1, length(A3b))
shrinkage3[1:2] <- shrinkage[1:2]
names(models3) <- c("model3", "model3i", "model3b", "model3bi")
out <- lapply(models3, use, A = A3b, c1 = shrinkage3)
saveRDS(out, "models3.RDS")


models <- list.files(pattern = "models[0-9].RDS")
models <- lapply(models, readRDS)
models <- do.call(c, models)
out <- lapply(names(models), function(x) {
  cbind.data.frame("RNAseq" = models[[x]]$Y[[1]][, 1],
                   "Micro" = models[[x]]$Y[[2]][, 1],
                   model = x,
                   AVE_inner = models[[x]]$AVE$AVE_inner[[1]],
                   AVE_outer = models[[x]]$AVE$AVE_outer[1],
                   Sample = rownames(models[[x]]$Y[[1]]))
})

out2 <- Reduce(rbind, out)
out3 <- merge(out2, A$Meta, by.x = "Sample", by.y = "Original", all.x = TRUE,
              sort = TRUE)
out3$Interaction <- ifelse(grepl("i$", out3$model), 1, 0)
out3$model <- gsub("i$", "", out3$model)
out3$model <- gsub("model", "", out3$model)
out3$model <- gsub("b$", " best", out3$model)
theme_update(strip.background = element_blank())
comm <- ggplot(out3[out3$Interaction != 1, ], aes(RNAseq, Micro)) +
  facet_wrap(~model, scale = "free")
comm +
  geom_point(aes(color = IBD))
comm +
  geom_point(aes(color = RNA_seq_batch))
comm +
  geom_point(aes(color = Exact_location))
comm +
  geom_point(aes(color = AgeDiag))
comm +
  geom_point(aes(color = diagTime))
comm +
  geom_point(aes(color = SEX))


genes <- sapply(names(models), function(x) {
  as.numeric(models[[x]]$a[[1]][, 1] != 0)
})

UpSetR::upset(as.data.frame(genes)[, grep("i$", colnames(genes), invert = TRUE)],
              nsets = 10, order.by = "freq", keep.order = TRUE)
micro <- sapply(names(models), function(x) {
  as.numeric(models[[x]]$a[[2]][, 1] != 0)
})

UpSetR::upset(as.data.frame(micro)[, grep("i$", colnames(micro), invert = TRUE)],
              nsets = 10, order.by = "freq", keep.order = TRUE)

sapply(models, function(x){x$AVE$AVE_inner[1]})
