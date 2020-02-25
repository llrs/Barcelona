library("integration")
library("RGCCA")
library("dplyr")
library("ggplot2")
library("UpSetR")
library("grid")
library("pROC")

theme_set(theme_bw())
theme_update(strip.background = element_blank())

models0 <- readRDS("models0_colon.RDS")
model0 <- models0[[1]]
model0i <- models0[[2]]
models2 <- readRDS("models2_colon.RDS")
model1 <- models2[[1]]
model1i <- models2[[2]]
model2 <- models2[[3]]
model2i <- models2[[4]]
# bests according to TRIM
# model2_best <- models2[[5]]
# model2_besti <- models2[[6]]
model2_best <- readRDS("model2b_colon_sgcca.RDS") # Best according to antiTNF
models3 <- readRDS("models3_colon.RDS")
model3 <- models3[[1]]
model3i <- models3[[2]]
# bests according to TRIM
# model3_best <- models3[[3]]
# model3_besti <- models3[[4]]

model3_best <- readRDS("model3_colon_best.RDS") # Best according to antiTNF

A <- readRDS("data/RGCCA_data.RDS")
meta <- A$Meta %>%
  select(Original:UC_endoscopic_remission) %>% # Remove processing columns
  mutate(Ileum = case_when(Exact_location == "ileum" ~ "Ileum",
                           !is.na(Exact_location) ~ "Colon"),
         IBD = case_when(is.na(IBD) ~ "CONTROL",
                         IBD == "UC" ~ "UC",
                         IBD == "CD" ~ "CD"),
         inv = case_when(Exact_location == tolower(gsub(" COLON$", "", Aftected_area)) ~ "involved",
                         TRUE ~ "not_involved")) %>%
  select_if(function(x)any(!is.na(x))) %>% # Select just those with information
  filter(Ileum == "Colon")

# Look the samples ####
a0GE <- cbind.data.frame(tidyer(model0$Y[[1]], "0", "GE"), meta)
a0M <- cbind.data.frame(tidyer(model0$Y[[2]], "0", "M"), meta)
a0iGE <- cbind.data.frame(tidyer(model0i$Y[[1]], "0 i", "GE"), meta)
a0iM <- cbind.data.frame(tidyer(model0i$Y[[2]], "0 i", "M"), meta)
a1GE <- cbind.data.frame(tidyer(model1$Y[[1]], "1", "GE"), meta)
a1M <- cbind.data.frame(tidyer(model1$Y[[2]], "1", "M"), meta)
a1iGE <- cbind.data.frame(tidyer(model1i$Y[[1]], "1 i", "GE"), meta)
a1iM <- cbind.data.frame(tidyer(model1i$Y[[2]], "1 i", "M"), meta)
a2GE <- cbind.data.frame(tidyer(model2$Y[[1]], "1.1", "GE"), meta)
a2M <- cbind.data.frame(tidyer(model2$Y[[2]], "1.1", "M"), meta)
a2bGE <- cbind.data.frame(tidyer(model2_best$Y[[1]], "1.2", "GE"), meta)
a2bM <- cbind.data.frame(tidyer(model2_best$Y[[2]], "1.2", "M"), meta)
a3GE <- cbind.data.frame(tidyer(model3$Y[[1]], "2.1", "GE"), meta)
a3M <- cbind.data.frame(tidyer(model3$Y[[2]], "2.1", "M"), meta)
a3bGE <- cbind.data.frame(tidyer(model3_best$Y[[1]], "2.2", "GE"), meta)
a3bM <- cbind.data.frame(tidyer(model3_best$Y[[2]], "2.2", "M"), meta)
#
# a1GE <- cbind("Model" = "1", a1GE)
# a1iGE <- cbind("Model" = "1 i", a1iGE)
# a1M <- cbind("Model" = "1", a1M)
# a1iM <- cbind("Model" = "1 i", a1iM)

inter <- intersect(colnames(a0GE), colnames(a0M))
# inter <- grep("Rownames", inter, invert = TRUE, value = TRUE)


df <- rbind(
  merge(a0M, a0GE, all = TRUE, by = inter),
  merge(a0iM, a0iGE, all = TRUE, by = inter),
  merge(a1M, a1GE, all = TRUE, by = inter),
  merge(a1iM, a1iGE, all = TRUE, by = inter),
  merge(a2M, a2GE, all = TRUE, by = inter),
  merge(a2bM, a2bGE, all = TRUE, by = inter),
  merge(a3M, a3GE, all = TRUE, by = inter),
  merge(a3bM, a3bGE, all = TRUE, by = inter)
)

# Plots ####
# * AUC: Ileum/Colon ####
rock0 <- roc(meta$Ileum, model0$Y$RNAseq[, 1])
rock1.2 <- roc(meta$Ileum, model2_best$Y$RNAseq[, 1])
rock2.2 <- roc(meta$Ileum, model3_best$Y$RNAseq[, 1])
par(pty = "s")
plot(rocks0, asp = 1, xlim = c(0, 1), ylim = c(0, 1), main = "AUC", col = "red")
lines(rocks1.2, asp = 1, xlim = c(0, 1), ylim = c(0, 1), col = "green")
lines(rocks2.2, asp = 1, xlim = c(0, 1), ylim = c(0, 1), col = "blue")
text(x = 0.5, y = 0.5, labels = paste("AUC =", round(rock0$auc, 3)), col = "red")
text(x = 0.5, y = 0.6, labels = paste("AUC =", round(rock1.2$auc, 3)), col = "green")
text(x = 0.5, y = 0.7, labels = paste("AUC =", round(rock2.2$auc, 3)), col = "blue")
legend("bottomleft",
       fill = c("red", "green", "blue"),
       legend = c("Model 0", "Model 1.2", "Model 2.2"))

# * Aim bullet ####

rownames(model3_best$a$Location) <- paste0("Loc", seq_len(nrow(model3_best$a$Location)))
v <- variables(model3_best)
plot_variables(v)

# * Model 2.2 plots ####
df3b <- data.frame(R = model3_best$Y[[1]][, 1],
                   M = model3_best$Y[[2]][, 1],
                   D = model3_best$Y[[3]][, 1],
                   L = model3_best$Y[[4]][, 1],
                   T = model3_best$Y[[5]][, 1],
                   Rownames = rownames(model3_best$Y[[2]]))

df3b <- cbind(df3b, meta)
ggplot(df3b) +
  geom_point(aes(D, M, color = ID)) +
  scale_color_viridis_d() +
  labs(x = "Demographics", y = "Microbiome", title = "Model 2.2")
ggplot(df3b) +
  geom_point(aes(D, R, color = ID)) +
  scale_color_viridis_d() +
  labs(x = "Demographics", y = "Transcriptome", title = "Model 2.2")
df3b %>%
  ggplot() +
  geom_point(aes(D, R, color = Ileum)) +
  scale_color_viridis_d() +
  labs(x = "Demographics", y = "Transcriptome", title = "Model 2.2")
df3b %>%
  ggplot() +
  geom_point(aes(T, R, color = Ileum)) +
  scale_color_viridis_d() +
  labs(x = "Time", y = "Transcriptome", title = "Model 2.2")

df3b %>%
  ggplot() +
  geom_point(aes(T, D, color = ID)) +
  scale_color_viridis_d() +
  labs(x = "Time", y = "Demographics", title = "Model 2.2")

df3b %>%
  ggplot() +
  geom_point(aes(L, R, color = Ileum)) +
  scale_color_viridis_d() +
  labs(x = "Location", y = "Transcriptome", title = "Model 2.2",
       color = "Location")
df3b %>%
  ggplot() +
  geom_point(aes(L, M, color = Ileum)) +
  scale_color_viridis_d() +
  labs(x = "Location", y = "Microbiome", title = "Model 2.2",
       color = "Location")

df3b %>%
  ggplot() +
  geom_point(aes(L, R, color = Exact_location)) +
  scale_color_viridis_d() +
  labs(x = "Location", y = "Microbiome", title = "Model 2.2",
       color = "Location")

ggplot(df3b) +
  geom_point(aes(D, M, color = IBD)) +
  labs(x = "Demographics", y = "Microbiome", title = "Model 2.2")
df3b %>%
  ggplot() +
  geom_point(aes(L, R, color = Ileum)) +
  scale_color_viridis_d() +
  labs(x = "Location", y = "Transcriptome", title = "Model 2.2",
       color = "Location")
ggplot(df3b) +
  geom_point(aes(L, R, color = IBD))
ggplot(df3b) +
  geom_point(aes(T, D, color = AgeDiag))
ggplot(df3b) +
  geom_point(aes(M, D, color = AgeDiag))

# * All models plots ####
df %>%
  filter(!grepl(" i", Model),
         Component == "comp1") %>%
  ggplot() +
  geom_point(aes(GE, M, color = Ileum, shape = Ileum)) +
  labs(color = "Location", shape = "Location") +
  facet_wrap(~Model, scales = "free")
loc <- last_plot()
df %>%
  filter(!grepl(" i", Model),
         Component == "comp1") %>%
  ggplot() +
  geom_point(aes(GE, M, color = IBD, shape = IBD)) +
  labs(color = "Disease", shape = "Disease") +
  facet_wrap(~Model, scales = "free")
dis <- last_plot()

# Look the weights ####
a0GE <- tidyer(model0$a[[1]], "0", "GE")
a0M <- tidyer(model0$a[[2]], "0", "M")
a0iGE <- tidyer(model0i$a[[1]], "0 i", "GE")
a0iM <- tidyer(model0i$a[[2]], "0 i", "M")
a1GE <- tidyer(model1$a[[1]], "1", "GE")
a1M <- tidyer(model1$a[[2]], "1", "M")
a1iGE <- tidyer(model1i$a[[1]], "1 i", "GE")
a1iM <- tidyer(model1i$a[[2]], "1 i", "M")
a2GE <- tidyer(model2$a[[1]], "2", "GE")
a2M <- tidyer(model2$a[[2]], "2", "M")
a2bGE <- tidyer(model2_best$a[[1]], "2 best", "GE")
a2bM <- tidyer(model2_best$a[[2]], "2 best", "M")
a3GE <- tidyer(model3$a[[1]], "3", "GE")
a3M <- tidyer(model3$a[[2]], "3", "M")
a3bGE <- tidyer(model3_best$a[[1]], "3 best", "GE")
a3bM <- tidyer(model3_best$a[[2]], "3 best", "M")

a0GE <- cbind("Model" = "0", a0GE)
a0M <- cbind("Model" = "0", a0M)
a0iGE <- cbind("Model" = "0 i", a0iGE)
a0iM <- cbind("Model" = "0 i", a0iM)
a1GE <- cbind("Model" = "1", a1GE)
a1M <- cbind("Model" = "1", a1M)
a1iGE <- cbind("Model" = "1 i", a1iGE)
a1iM <- cbind("Model" = "1 i", a1iM)
a3GE <- cbind("Model" = "3", a3GE)
a3M <- cbind("Model" = "3", a3M)
a2GE <- cbind("Model" = "2", a2GE)
a2M <- cbind("Model" = "2", a2M)


dfGE <- rbind(a0GE, a0iGE, a1GE, a1iGE, a2GE, a2bGE, a3GE, a3bGE)
dfM <- rbind(a0M, a0iM, a1M, a1iM, a2M, a2bM, a3M, a3bM)
keepGE <- dfGE %>%
  # filter(!grepl(" i", Model)) %>%
  filter(Component == "V1" & GE != 0) %>%
  mutate(Presence = if_else(GE != 0, 1, 0)) %>%
  dplyr::select(-Component, -GE, Rownames) %>%
  group_by(Rownames) %>%
  tidyr::spread(Model, Presence) %>%
  as.data.frame()
rownames(keepGE) <- keepGE$Rownames
keepGE <- keepGE[, -grep("Rownames", colnames(keepGE))]
keepGE[is.na(keepGE)] <- 0

keepM <- dfM %>%
  # filter(!grepl(" i", Model)) %>%
  filter(Component == "V1" & M != 0) %>%
  mutate(Presence = if_else(M != 0, 1, 0)) %>%
  dplyr::select(-Component, -M, Rownames) %>%
  group_by(Rownames) %>%
  tidyr::spread(Model, Presence) %>%
  ungroup() %>%
  as.data.frame()
rownames(keepM) <- keepM$Rownames
keepM <- keepM[, -grep("Rownames", colnames(keepM))]
keepM[is.na(keepM)] <- 0


text_sizes <- c(1.3, 1.3, 1, 1, 1.5, 1.5)
dfGE %>%
  filter(Component == "V1" & GE != 0) %>%
  ggplot() +
  geom_density(aes(GE)) +
  facet_wrap("Model") +
  labs(title = "Distribution of the weights", xlab = "weights",
       subtitle = "Gene expression")
dfM %>%
  filter(Component == "V1" & M != 0) %>%
  ggplot() +
  geom_density(aes(M)) +
  facet_wrap("Model") +
  labs(title = "Distribution of the weights", xlab = "weights",
       subtitle = "Microbiome")

# * Upset plots ####
upset(keepGE, order.by = "freq", nsets = 6,
      sets = rev(colnames(keepGE)), keep.order = TRUE,
      line.size = NA, text.scale = text_sizes, scale.sets = "identity")
grid.text("Genes shared in models", x = 0.65, y = 0.95, gp = gpar(fontsize = 20))
upset(keepGE, order.by = "freq", keep.order = TRUE, sets = rev(c("0", "1" ,"2", "2 best", "3", "3 best")),
      line.size = NA, text.scale = text_sizes, scale.sets = "identity")
grid.text("Genes shared in models", x = 0.65, y = 0.95, gp = gpar(fontsize = 20))
upset(keepM, order.by = "freq", nsets = 6,
      sets = rev(colnames(keepM)), keep.order = TRUE,
      line.size = NA, text.scale = text_sizes)
grid.text("OTUs shared in models", x = 0.65, y = 0.95, gp = gpar(fontsize = 20))
upset(keepM, order.by = "freq", keep.order = TRUE, sets = rev(c("0", "1" ,"2", "2 best", "3", "3 best")),
      line.size = NA, text.scale = text_sizes, scale.sets = "identity")
grid.text("OTUs shared in models", x = 0.65, y = 0.95, gp = gpar(fontsize = 20))

