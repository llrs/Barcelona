library("integration")
library("RGCCA2")
library("dplyr")

A <- readRDS("data/RGCCA_data.RDS")
meta <- A$Meta

testing <- function(x, ...) {
  tryCatch({
  result.sgcca <- RGCCA2::sgcca(C = x,
                                scheme = "centroid",
                                verbose = FALSE,
                                scale = FALSE,
                                ...)
  analyze(result.sgcca)}
  , error = function(x){NA})
}

Localization <- model_RGCCA(A$Meta, c("Exact_location")) # With SESCD local it increase the AVE_inner
Time <- model_RGCCA(A$Meta, c("AgeDiag", "Age"))
Demographics <- model_RGCCA(A$Meta, c("ID","SEX"))
Time$AgeDiag[is.na(Time$AgeDiag)] <- 0 # Time has NA values
A2 <- A[1:2]
A2$Demographics <- Demographics
A2$Location <- Localization
A2$Time <- Time

A2 <- clean_unvariable(A2)
saveRDS(A2, "model3_BCN.RDS")

shrinkage <- sapply(A2[1:2], tau.estimate) # 0.11503779803812 0.318145965316924
shrinkage[3:5] <- 1
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
set.seed(46726279)
s <- sample(designs, size = min(length(designs)*.1, 1000))

out <- sapply(s, testing, A = Ab, c1 = shrinkage, USE.NAMES = FALSE)
out2 <- out[lengths(out) == 24]
out2 <- simplify2array(out2)
out2 <- as.data.frame(t(out2))
saveRDS(out2, "sample_model3_boot.RDS")

out %>%
  top_n(5, AVE_inner) %>%
  select(AVE_inner, AVE_outer, var12, var13, var23,
         var14, var24, var34, var15, var25, var35, var45) %>%
  arrange(desc(AVE_inner))
stop("Visual inspection of the top 5")

s2 <- sample(designs, size = min(length(designs)*.1, 1000))

out <- sapply(s2, testing, A = Ab, c1 = shrinkage, USE.NAMES = FALSE)
out2 <- out[lengths(out) == 24]
out2 <- simplify2array(out2)
out2 <- as.data.frame(t(out2))
saveRDS(out2, "sample2_model3_boot.RDS")
stop("Visual inspection of the top 5")

out1 <- readRDS("sample_model3_boot.RDS")
out0 <- rbind(out1, out2)
out <- out0[!duplicated(out0), ]
ggplot(out, aes(AVE_inner, AVE_outer)) +
  geom_point()

keep_best <- vapply(designs, function(x){
  x[1, 2] == 0 & x[2, 3] == 0
}, logical(1L))

out <- sapply(designs[keep_best], testing, A = Ab, c1 = shrinkage, USE.NAMES = FALSE)
out2 <- out[lengths(out) == 24]
out2 <- simplify2array(out2)
out2 <- as.data.frame(t(out2))
saveRDS(out2, "sample_def_model3_boot.RDS")

out0 <- readRDS("sample_model3_boot.RDS")
out1 <- readRDS("sample2_model3_boot.RDS")
out2 <- readRDS("sample_def_model3_boot.RDS")
out <- rbind(out0, out1, out2)
out <- out[!duplicated(out), ]
ggplot(out, aes(AVE_inner, AVE_outer, color = cc1)) +
  geom_point()
best3 <- out[out$AVE_inner == max(out$AVE_inner), grep("var", colnames(out))]
best3 <- symm(designs[[1]], best3)
colnames(best3) <- names(Ab)
rownames(best3) <- names(Ab)

model3_best <- sgcca(A = Ab, c1 = shrinkage, C = best3, ncomp = rep(2, 5), scheme = "centroid")
saveRDS(model3_best, "model3_best.RDS")
