library("clusterProfiler")
library("org.Hs.eg.db")
library("magrittr")
library("ReactomePA")
d <- read.csv("data_out/weights_2.1.csv")
d <- d[, -1]
entrez <- mapIds(org.Hs.eg.db, keys = as.character(d$r), column = "ENTREZID", keytype = "ENSEMBL")
symbol <- mapIds(org.Hs.eg.db, keys = as.character(d$r), column = "SYMBOL", keytype = "ENSEMBL")
d2 <- cbind(d, name = symbol)
write.csv(d2, file = "data_out/weight_named_2.1.csv",
          row.names = FALSE, na = "")

# Data downloaded from wikipathways on 07/04/2020
wp2gene <- read.gmt("data/wikipathways-20200310-gmt-Homo_sapiens.gmt") %>%
  tidyr::separate(ont, c("name","version","wpid","org"), "%")

e <- clusterProfiler::enricher(entrez,
                               TERM2GENE = wp2gene[, c("wpid", "gene")],
                               TERM2NAME = wp2gene[, c("wpid", "name")])
df <- as.data.frame(e)
ek <- ReactomePA::enrichPathway(entrez)
ek <- as.data.frame(ek)
ek[, 1:5]

def <- rbind(ek[, 1:5], e[, 1:5])
def[order(def$pvalue), ]

write.csv(def, file = "data_out/pathways_model_1.2.csv", row.names = FALSE)
