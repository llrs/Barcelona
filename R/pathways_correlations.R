library("integration")
library("org.Hs.eg.db")
library("psych")
library("dplyr")
library("tidyr")
library("clusterProfiler")
library("msigdbr")
library("ReactomePA")
library("ComplexUpset")
library("ggplot2")

m_t2g <- msigdbr(species = "Homo sapiens", category = "C7") %>%
  dplyr::select(gs_name, entrez_gene)
msig <- msigdbr(species = "Homo sapiens", category = "C5")
m_t2g_c5 <-  msig %>%
  dplyr::select(gs_name, entrez_gene)

genes_interest <- c("IL4", "IL2", "IL8", "IL12", "S100A8",
                    "ADCYAP1", "GMPR", "HDC", "TPSAB1", "APQ8", "AQP7", "CHI3L1",
                    "COL3A1", "DERL3", "GZMH", "HTR3E", "IFN GAMMA", "IL17A", "OSM",
                    "PDGFD", "PTGDR2", "RTNLB",  "SOX6", "TBX21", "THY1")


pathways <- msig %>%
  filter(gene_symbol %in% genes_interest) %>%
  distinct(gs_name) %>%
  pull(gs_name)
genes_pathways <- msig %>%
  filter(gs_name %in% pathways)
GOI2 <- mapIds(org.Hs.eg.db, keys = as.character(genes_pathways$entrez_gene),
              keytype = "ENTREZID", column = "ENSEMBL")
A <- readRDS("data/RGCCA_data_wo_out.RDS")

RNAseq <- A$RNAseq
micro <- A$micro

RNAseq2 <- RNAseq[, trimVer(colnames(RNAseq)) %in% GOI2]
corr2 <- readRDS("output/correlations_genes_pvalue_ASV.RDS")
# p <- corr.test(RNAseq2, micro)

corr3 <- lapply(corr2, function(x, y) {x[y, ]}, y = colnames(RNAseq2))
r3 <- corr3$r %>%
  as.data.frame() %>%
  tibble::rownames_to_column() %>%
  pivot_longer(cols = -rowname) %>%
  rename(Microorganism = name, r = value)
p3 <- corr3$p %>%
  as.data.frame() %>%
  tibble::rownames_to_column() %>%
  pivot_longer(cols = -rowname) %>%
  rename(Microorganism = name, adj.p.value = value)
tidy_corr <- full_join(r3, p3)
