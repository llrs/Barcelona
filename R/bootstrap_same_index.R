# See the modelling_overfitted.R file
library("inteRmodel")
library("integration")
library("RGCCA")
library("ggplot2")

# Index ####
A <- readRDS("data/RGCCA_data.RDS")
boots <- 1000
set.seed(9876156)
index <- vector("list", length = boots)
for (i in seq_len(boots)) {
  index[[i]] <- sample(nrow(A[[1]]), replace = TRUE)
}
saveRDS(index, file = "data_out/index_boot.RDS")
index <- readRDS("data_out/index_boot.RDS")

# * Model 0 ####
shrinkage <- c(0.285693348851943, 0) #Calculated from the server for the data derived from original data
shrinkage[2] <- tau.estimate(A[[2]])
Ab <- lapply(A[1:2], function(x) scale2(x, bias = TRUE)/sqrt(NCOL(x)))
names(Ab) <- names(A[1:2])
ab <- clean_unvariable(Ab)
b0 <- boot_index_sgcca(index, A = ab, c1 = shrinkage, scheme = "centroid",
                         scale = FALSE, verbose = FALSE, bias = TRUE)
saveRDS(b0, data_out/"boot_0.RDS")

# * Model 1.2 ####
model2_best <- readRDS("data_out/model2b_sgcca.RDS")
C1.2 <- model2_best$C

A$Meta <- model_RGCCA(A$Meta, c("ID", "AgeDiag", "diagTime", "Exact_location", "SEX"))
Ab <- lapply(A, function(x) scale2(x, bias = TRUE)/sqrt(NCOL(x)))
Ab <- clean_unvariable(Ab)
shrinkage <- rep(1, 3)
shrinkage[1:2] <- c(0.11503779803812, 0.318145965316924)


b1.2 <- boot_index_sgcca(index, A = Ab, c1 = shrinkage, scheme = "centroid",
                          scale = FALSE, verbose = FALSE, bias = TRUE, C = C1.2)
saveRDS(b1.2, "data_out/boot_1.2.RDS")

# * Model 2.2 ####
model3_best <- readRDS("data_out/model3_best_treatment.RDS")
A <- readRDS("data_out/model3_BCN.RDS")
C2.2 <- model3_best$C

shrinkage <- rep(1, 5)
names(shrinkage) <- names(A)
shrinkage[1:2] <- c(0.11503779803812, 0.318145965316924)
shrinkage[3:5] <- 1
names(shrinkage) <- names(A)
Ab <- lapply(A, function(x) scale2(x, bias = TRUE)/sqrt(NCOL(x)))
ab <- clean_unvariable(Ab)

b2.2 <- boot_index_sgcca(index, A = ab, C = C2.2, c1 = shrinkage, scheme = "centroid")
saveRDS(b2.2, file = "boot_2.2.RDS")


# AVE plot ####
b0 <- readRDS("data_out/boot_0.RDS")
b1.2 <- readRDS("data_out/boot_1.2.RDS")
b2.2 <- readRDS("data_out/boot_2.2.RDS")

theme_set(theme_bw())
theme_update(strip.background = element_blank(),
             panel.grid.minor = element_blank())

AVE0 <- b0$AVE
AVE1.2 <- b1.2$AVE
AVE2.2 <- b2.2$AVE
b <- rbind.data.frame(cbind.data.frame(AVE0, model = "0"),
                      cbind.data.frame(AVE1.2, model = "1.2"),
                      cbind.data.frame(AVE2.2, model = "2.2"))
b$index <- rep(seq_len(1000), 3)

models0 <- readRDS("data_out/models0.RDS")
model0 <- models0[[1]]

# This index doesn't converge
b <- b[!(b$inner == 0 & b$outer == 0), ]
ggplot(b) +
  geom_density(aes(inner, group = model, fill = model), alpha = 0.5)
ggplot(b) +
  geom_density(aes(outer, group = model, fill = model), alpha = 0.5)

AVE_names <- c("AVE_inner", "AVE_outer")
b$model <- as.character(b$model)
# * plot ####
p <- ggplot(b) +
  geom_point(aes(inner, outer, col = model, shape = model), alpha = 0.5) +
  geom_point(aes(AVE_inner, AVE_outer),
             data = as.data.frame(model3_best$AVE[AVE_names])[1, , drop = FALSE],
             fill = "blue", col = "black", shape = 21) +
  geom_point(aes(AVE_inner, AVE_outer),
             data = as.data.frame(model2_best$AVE[AVE_names])[1, , drop = FALSE],
             fill = "green", col = "black", shape = 21) +
  geom_point(aes(AVE_inner, AVE_outer),
             data = as.data.frame(model0$AVE[AVE_names])[1, , drop = FALSE],
             fill = "red", col = "black", shape = 21) +
  stat_ellipse(aes(inner, outer, col = model)) +
  labs(title = "AVE in bootstraps", x = "Inner AVE", y = "Outer AVE") +
  theme(legend.position = "bottom")
