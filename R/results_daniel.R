library("openxlsx")

# Read
colon_CD <- read_xls("output/results_060420_antiTNF/colon.CD.xls")
colon_UC <- read_xls("output/results_060420_antiTNF/colon.UC.xls")
ileum_CD <- read_xls("output/results_060420_antiTNF/ileum.CD.xls")

# Clean names
clean_names <- function(x) {
  x <- unlist(x)
  names(x) <- NULL
  x <- gsub(x, pattern = "-w", replacement = "-w00")
  x <- gsub(x, pattern = "-w00(1|4)", replacement = "-w0\\1")
  x <- gsub(x, pattern = "^([0-9])", replacement = "00\\1")
  x <- gsub(x, pattern = "^00?([0-9]{3})", replacement = "\\1")
  gsub(x, pattern = "-T-TR-", replacement = "-T-DM-")
}

s_colon_CD <- clean_names(colon_CD[5, 2:48])
s_colon_UC <- clean_names(colon_UC[4, 2:51])
s_ileum_CD <- clean_names(ileum_CD[5, 2:40])


meta <- readRDS("output/refined_meta.RDS")
meta$ileum <- ifelse(meta$Exact_location == "ileum", "ileum", "colon")

scUC <- meta$Original[meta$IBD %in% c("CONTROL", "UC") & meta$ileum == "colon"]
all(scUC %in% s_colon_UC) # All sample in
any(!s_colon_UC %in% scUC) # No extra sample
scCD <- meta$Original[meta$IBD %in% c("CONTROL", "CD") & meta$ileum == "colon"]
all(scCD %in% s_colon_CD) # All sample in
any(!s_colon_CD %in% scCD) # No extra sample
siCD <- meta$Original[meta$ileum == "ileum"]
all(siCD %in% s_ileum_CD) # All sample in
any(!s_ileum_CD %in% siCD) # No extra sample

# Conclusion, the samples are correct on the right file
# I assume they are correctly assigned on the contrasts
