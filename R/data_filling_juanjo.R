library("readxl")
library("stringr")
library("purrr")
library("dplyr")


nom <- read_xls("~/Documents/users/juanjo_lozano/aTNF/noms.xls")
pheno <- read.csv("data_out/pheno_data.csv", row.names = 1)
nam <- gsub(" reseq", "", nom$nom)
bcn <- grep("-[wW]", nam)
nam[bcn] <- tolower(nam[bcn])


correct_bcn <- function(x) {
  if (length(x) > 1) {
    a <- str_pad(x[1], width = 3, pad = "0")
    b <- str_pad(x[2], width = 3, pad = "0")
    x <- paste(a, b, sep = "-w")
  }
  x
}

nam <- nam[bcn] %>%
  str_split("-w") %>% # Ready for BCN
  map(correct_bcn) %>%
  unlist()
m <- merge(data.frame(nam), pheno, by.x = "nam", by.y = "Original",
           sort = FALSE, all.x = TRUE, all.y = FALSE)
m %>%
  group_by(nam) %>%
  map( ~ function(x){
    y <- unique(x$Exact_location)
    y[!is.na(y)]
    z <- unique(x$Segmento)
    z[!is.na(z)]

  })
write.csv(m[, c("nam", "Exact_location", "Segmento", "IBD")], "data_out/nams.csv", row.names = FALSE)
