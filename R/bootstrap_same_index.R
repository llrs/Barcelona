# See the modelling_overfitted.R file
library("inteRmodel")
library("integration")
library("RGCCA")
library("ggplot2")
library("BiocParallel")

# Index ####
A <- readRDS("data/RGCCA_data_wo_out.RDS")
boots <- 10000
set.seed(9876156)
# index <- inteRmodel::boot_index(nrow(A[[1]]), boots)
# saveRDS(index, file = "data_out/index_boot2.RDS")
index <- readRDS("data_out/index_boot2.RDS")

# * Model 0 ####
shrinkage <- c(0.322297910454825, 0.866155549496009) #Calculated from the server for the data derived from original data
shrinkage[2] <- tau.estimate(A[[2]])
Ab <- lapply(A[1:2], function(x) scale2(x, bias = TRUE)/sqrt(NCOL(x)))
names(Ab) <- names(A[1:2])
ab <- clean_unvariable(Ab)
C0 <- matrix(c(0, 1, 1, 0), nrow = 2)
mcp <- MulticoreParam(workers = 8, progressbar = TRUE)
b0 <- boot_index_sgcca(index, A = ab, c1 = shrinkage, scheme = "centroid",
                         scale = FALSE, verbose = FALSE, bias = TRUE, C = C0,
                       BPPARAM = mcp)
saveRDS(b0, "data_out/boot_0_b.RDS")

# * Model 1.2 ####
model2_best <- readRDS("data_out/model2b2_sgcca_b.RDS")
C1.2 <- model2_best$C

A$Meta <- model_RGCCA(A$Meta, c("ID", "AgeDiag", "diagTime", "Exact_location", "SEX"))
Ab <- lapply(A, function(x) scale2(x, bias = TRUE)/sqrt(NCOL(x)))
Ab <- clean_unvariable(Ab)
shrinkage <- rep(1, 3)
shrinkage[1:2] <- c(0.322020648273615, 0.866155549496009)


b1.2 <- boot_index_sgcca(index, A = Ab, c1 = shrinkage, scheme = "centroid",
                          scale = FALSE, verbose = FALSE, bias = TRUE, C = C1.2,
                         BPPARAM = mcp)
saveRDS(b1.2, "data_out/boot_1.2_b.RDS")

# * Model 2.2 ####
model3_best <- readRDS("data_out/model3_best_treatment_b.RDS")
# A <- readRDS("data_out/model3_BCN_b.RDS")
C2.2 <- model3_best$C
meta <- A$Meta

Localization <- model_RGCCA(meta, c("Exact_location")) # With SESCD local it increase the AVE_inner
Time <- model_RGCCA(meta, c("AgeDiag", "Age"))
Demographics <- model_RGCCA(meta, c("ID","SEX"))
Time$AgeDiag[is.na(Time$AgeDiag)] <- 0 # Time has NA values
Time$aTNF <- ifelse(meta$Time == "0" | is.na(meta$Time), 0, 1)
A2 <- A[1:2]
A2$Demographics <- Demographics
A2$Location <- Localization
A2$Time <- Time

shrinkage <- rep(1, 5)
names(shrinkage) <- names(A2)
shrinkage[1:2] <- c(0.322020648273615, 0.866155549496009)
shrinkage[3:5] <- 1
names(shrinkage) <- names(A2)
Ab <- lapply(A2, function(x) scale2(x, bias = TRUE)/sqrt(NCOL(x)))
ab <- clean_unvariable(Ab)

b2.2 <- boot_index_sgcca(index, A = ab, C = C2.2, c1 = shrinkage,
                         scheme = "centroid", BPPARAM = mcp)
saveRDS(b2.2, file = "data_out/boot_2.2_b.RDS")


# AVE plot ####
b0 <- readRDS("data_out/boot_0_b.RDS")
b1.2 <- readRDS("data_out/boot_1.2_b.RDS")
b2.2 <- readRDS("data_out/boot_2.2_b.RDS")

theme_set(theme_bw())
theme_update(strip.background = element_blank(),
             panel.grid.minor = element_blank())

AVE0 <- b0$AVE
AVE1.2 <- b1.2$AVE
AVE2.2 <- b2.2$AVE
b <- rbind.data.frame(cbind.data.frame(AVE0, model = "0"),
                      cbind.data.frame(AVE1.2, model = "1.2"),
                      cbind.data.frame(AVE2.2, model = "2.2"))
b$index <- rep(seq_len(10000), 3)

models0 <- readRDS("data_out/models0_b.RDS")
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
ggsave(plot = p, filename = "Figures/bootstrap_same_index_wo.png", dpi = 300)
