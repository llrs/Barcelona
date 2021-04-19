# Run the spls method on the data 19/04/2021
library("spls")
A <- readRDS("data/RGCCA_data_wo_out.RDS")
meta <- A$Meta
cv <- cv.spls(A$RNAseq, A$Micro, eta = seq(0.1,0.9,0.1), K = c(5:10))
f <- spls(yeast$x, yeast$y, eta = cv$eta.opt, K = cv$K.opt )
