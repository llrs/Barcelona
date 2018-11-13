#' ---
#' title: "Experimental design of Michigan sequencing"
#' author: "Lluís Revilla"
#' date: "r `date() `"
#' output: html_document
#' ---

pheno <- readRDS("samples2sequence.RDS")
#' This should match with 51 CONTROLS (They are used for both studies),
#' 115 CD of HSCT, 205 (CD+UC),
#'
#' # Classify by plates
#'
#' The aimed plate is also 96, so we can use the
#' [OSAT](https://bioconductor.org/packages/OSAT/) package. However we need to
#' remove teh numeric values in order to use it:
library("OSAT")
#' error=FALSE
gs <- setup.sample(pheno, optimal=c("Sample Name_RNA", "IBD", "SEX", "Study",
                                    "ID", "Exact_location"))
gc <- setup.container(IlluminaBeadChip96Plate, 6, batch='plates')
gSetup <- create.optimized.setup(sample=gs, container=gc, nSim=1000)
QC(gSetup)
