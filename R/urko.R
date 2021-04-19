# Compare the samples on this study with the ones send to Urko
library("dplyr")
urko <- xlsx::read.xlsx("~/Documents/users/urko/database_mostres_Urko_COMPLETE_rev2.xlsx", sheetIndex = 1, stringsAsFactors = FALSE)

A <- readRDS("data/RGCCA_data_wo_out.RDS")
meta <- A$Meta
urko_bcn <- urko %>%
  filter(Study == "BCN",
         region != "blood") %>%
  mutate(Sample = paste0(Patient, "-", week))

setdiff(meta$Original, urko_bcn$Sample)
setdiff(urko_bcn$Sample, meta$Original)
length(intersect(meta$Original, urko_bcn$Sample))
