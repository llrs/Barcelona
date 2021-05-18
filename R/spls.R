# Run the spls method on the data 19/04/2021
library("spls")
A <- readRDS("data/RGCCA_data_wo_out.RDS")
meta <- A$Meta
cv <- cv.spls(A$RNAseq, A$Micro, eta = seq(0.1,0.9,0.1), K = c(5:10))
saveRDS(cv, "output/cv_spls.RDS")
f <- spls(yeast$x, yeast$y, eta = cv$eta.opt, K = cv$K.opt )
saveRDS(f, "output/f_spls.RDS")
spls(yeast$x, yeast$y, eta=cv$eta.opt, K=cv$K.opt)


library("mixOmics")

ncomp <- 2
result.spls <- mixOmics::spls(A$RNAseq, A$Micro, ncomp = ncomp,
                    # keepX = rep(4000, ncomp),
                    # keepY = rep(100, ncomp),
                    max.iter = 1000)
plotIndiv(result.spls, group = meta$Exact_location,
          rep.space = "XY-variate", legend = TRUE,
          legend.title = 'IBD',
          ind.names = meta$Name,
          title = 'aTNF: sPLS')

# tune.spls <- perf(result.spls, validation = 'Mfold', folds = 10,
                  # criterion = 'all', progressBar = FALSE)
