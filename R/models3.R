library("integration")
library("RGCCA")
library("dplyr")
library("inteRmodel")
library("ggplot2")
library("BiocParallel")

mcp <- MulticoreParam(workers = 8, progressbar = TRUE)
A <- readRDS("data/RGCCA_data_wo_out.RDS")
meta <- A$Meta

stopifnot(all(A$Meta$Original == rownames(A$RNAseq)))
testing <- function(x, ...) {
  tryCatch({
    result.sgcca <- RGCCA::sgcca(C = x,
                                 scheme = "centroid",
                                 verbose = FALSE,
                                 scale = FALSE,
                                 ...)
    analyze(result.sgcca)}
    , error = function(x){NA})
}

meta$`Phenotype CD` <- tolower(meta$`Phenotype CD`)
meta$`Location` <- tolower(meta$`Location`)

Localization <- model_RGCCA(meta, c("Exact_location")) # With SESCD local it increase the AVE_inner
Time <- model_RGCCA(meta, c("AgeDiag", "Age"))
Demographics <- model_RGCCA(meta, c("ID","SEX"))
Time$AgeDiag[is.na(Time$AgeDiag)] <- 0 # Time has NA values
Time$aTNF <- ifelse(meta$Time == "0" | is.na(meta$Time), 0, 1)
A2 <- A[1:2]
A2$Demographics <- Demographics
A2$Location <- Localization
A2$Time <- Time

A2 <- clean_unvariable(A2)
# saveRDS(A2, "data_out/model3_BCN_treatment.RDS")
#
shrinkage <- rep(1, 5)
names(shrinkage) <- names(A2)
Ab <- lapply(A2, function(x) scale2(x, bias = TRUE)/sqrt(NCOL(x)))
# shrinkage[1:2] <- vapply(A2[1:2], tau.estimate, numeric(1L)) # 0.322297910454825 0.959997006295785
# dput(shrinkage)
shrinkage[1:2] <- c(0.322020648273615, 0.866155549496009)

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
s <- sample(designs, size = min(length(designs)*.1, 10000))
out <- bplapply(s, testing, A = Ab, c1 = shrinkage, USE.NAMES = FALSE,
                BPPARAM = mcp)
out2 <- out[lengths(out) == 24]
out2 <- simplify2array(out2)
out2 <- as.data.frame(t(out2))
saveRDS(out2, "data_out/sample_model3_boot_treatment_b.RDS")
out1 <- readRDS("data_out/sample_model3_boot_treatment_b.RDS")
out %>%
  top_n(5, AVE_inner) %>%
  select(AVE_inner, AVE_outer, var12, var13, var23,
         var14, var24, var34, var15, var25, var35, var45) %>%
  arrange(desc(AVE_inner))
stop("Visual inspection of the top 5")
#
# s2 <- sample(designs, size = 10000)
#
# out <- bplapply(s2, testing, A = Ab, c1 = shrinkage, USE.NAMES = FALSE,
# BPPARAM = mcp)
# out2 <- out[lengths(out) == 24]
# out2 <- simplify2array(out2)
# out2 <- as.data.frame(t(out2))
# saveRDS(out2, "data_out/sample2_model3_boot_b.RDS")
# out2 <- readRDS("data_out/sample2_model3_boot_b.RDS")
#
# # out1 <- readRDS("data_out/sample_model3_boot_b.RDS")
# out0 <- rbind(out1, out2)
# out <- out0[!duplicated(out0), ]
# # ggplot(out, aes(AVE_inner, AVE_outer)) +
# #   geom_point()
# out %>%
#   top_n(5, AVE_inner) %>%
#   select(AVE_inner, AVE_outer, var12, var13, var23,
#          var14, var24, var34, var15, var25, var35, var45) %>%
#   arrange(desc(AVE_inner))
# # stop("Visual inspection of the top 5")
#
# keep_best <- vapply(designs, function(x){
#   x[2, 3] == 1 & x[1, 5] == 0 & x[2, 5] == 0 & x[3, 5] == 0 & x[4, 5] != 0 & x[1, 3] != 0
# }, logical(1L))
#
# out <- bplapply(designs[keep_best], testing, A = Ab, c1 = shrinkage, USE.NAMES = FALSE,
# BPPARAM = mcp)
# # out2 <- out[lengths(out) == 24]
# # out2 <- simplify2array(out2)
# # out2 <- as.data.frame(t(out2))
# # out <- as.data.frame(t(out))
# saveRDS(out, "data_out/sample_def_model3_boot_b.RDS")
# out <- readRDS("data_out/sample_def_model3_boot_b.RDS")
out0 <- readRDS("data_out/sample_model3_boot_b.RDS")
out1 <- readRDS("data_out/sample2_model3_boot_b.RDS")
out2 <- readRDS("data_out/sample_def_model3_boot_b.RDS")
out2 <- as.data.frame(t(out2))
out <- rbind(out0, out1, out2)
out <- out[!duplicated(out), ]
# ggplot(out, aes(AVE_inner, AVE_outer, color = cc1)) +
#   geom_point()
best3 <- out[out$AVE_inner == max(out$AVE_inner), grep("var", colnames(out))]
best3 <- symm(C, best3)
# colnames(best3) <- names(Ab)
# rownames(best3) <- names(Ab)
d <- weight_design(weights = 11, size = length(Ab), which(lower.tri(best3) & best3 != 0))
check_model <- function(x){
  x[1, 2] != 0 & x[1, 3] != 0 & x[2, 3] == 1 & x[1, 4] == 0 & x[2, 4] == 1 & x[3, 4] == 1 &
    x[1, 5] == 0 & x[2, 5] == 0 & x[3, 5] == 0 & x[4, 5] == 0.1
}
keep_best <- vapply(d, check_model, logical(1L))
d <- d[keep_best]
out <- bplapply(d, testing, A = Ab, c1 = shrinkage, USE.NAMES = FALSE,
                 BPPARAM = mcp)
saveRDS(out, "data_out/refined_model3_2.RDS")
out <- readRDS("data_out/refined_model3_2.RDS")
if (is.list(out)) {
  out <- out[lengths(out) == 24]
  out <- simplify2array(out)
  }
o <- as.data.frame(t(out))

best3 <- o[o$AVE_inner == max(o$AVE_inner), grep("var", colnames(o))]
best3 <- symm(C, best3)
model3_best <- sgcca(A = Ab, c1 = shrinkage, C = best3, ncomp = rep(2, 5), scheme = "centroid",
                     scale = FALSE)
model3_best <- improve.sgcca(model3_best, names(Ab))
beepr::beep()
stopifnot(model3_best$AVE$AVE_inner[1] != max(o$AVE_inner))
saveRDS(model3_best, "data_out/model3_best_treatment.RDS")
