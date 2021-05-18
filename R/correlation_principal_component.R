library("integration")
library("org.Hs.eg.db")
library("psych")
library("dplyr")
library("tidyr")
library("clusterProfiler")
library("msigdbr")
library("ReactomePA")
library("ComplexUpset")
library("ComplexHeatmap")
library("ggplot2")
library("inteRmodel")
library("tidyverse")
library("xlsx")

# Load data ####
A <- readRDS("data/RGCCA_data_wo_out.RDS")
models0 <- readRDS("output/models0_b.RDS")
model0 <- models0[[1]]
model2_best <- readRDS("output/model2b2_sgcca_b.RDS")
model3_best <- readRDS("output/model3_best_treatment_b.RDS")
microorganism <-  readRDS("output/taxonomy_ASV.RDS")$tax
microorganism <- cbind(microorganism, "ASV" = rownames(microorganism))
microorganism <- as.data.frame(microorganism, stringsAsFactors = FALSE)
microorganism$rowname <- as.character(seq_len(nrow(microorganism)))
rownames(microorganism) <- as.character(seq_len(nrow(microorganism)))
mit <- rownames(microorganism)[microorganism[, "Family"] == "Mitochondria"]

# Variables with their component ####
# Calculates correlation with their principal component.
cor_PC <- function(x, y, z) {
  l <- vector("list", length(x))
  consty <- 0.69
  constx <- 3.58
  for (i in seq_along(x)) {
    if (is.null(z)) {
      xx <- x[[i]]
    }
    # Offset to get positive numbers
    # consulta Juanjo (idea inicial sumar 0.01)
    xx <- x[[i]][ ,z[[i]][, 1] != 0] + constx
    l[[i]] <- corr.test(log10(xx), log10(y[[i]][, 1] + consty), method = "spearman", ci = FALSE)
  }
  names(l) <- names(x)
  l
}

cPC_m0 <- cor_PC(A[1:2], model0$Y, model0$a)
cPC_m1.2 <- cor_PC(A[1:2], model2_best$Y, model2_best$a)
cPC_m2.2 <- cor_PC(A[1:2], model3_best$Y, model3_best$a)
l <- list(m0 = cPC_m0, m1.2 = cPC_m1.2, m2.2 = cPC_m2.2)
saveRDS(l, "output/cPC2.RDS")
l <- readRDS("output/cPC2.RDS")

# *Select genes and microorganisms ####
# Parameter to control the code below #
p <- 20
{
prop_analysis <- p/100

df_cP <- function(x) {
  data.frame(r = x$r[, 1], p = x$p[, 1])
}

l2 <- lapply(l, function(x) {
  l0 <- lapply(x, function(y, prop) {
    z <- df_cP(y) %>%
      filter(p < 0.05 )
    message(nrow(z))
    z %>%
      arrange(p, -abs(r)) %>%
      slice_max(p, prop = prop) %>%
      tibble::rownames_to_column() %>%
      mutate(rowname = trimVer(rowname))
  }, prop = prop_analysis)
  do.call(rbind, l0)
})

for (name in names(l2)) {
  x <- l2[[name]]
  gene <- startsWith(x$rowname, "ENSG")
  x$Type <- ifelse(gene, "Gene", "Microorganism")
  x[, "Model"] <- name
  rownames(x) <- NULL
  l2[[name]] <- x
}

corr_models <- do.call(rbind, l2)
rownames(corr_models) <- NULL
genes <- unique(corr_models$rowname[corr_models$Type == "Gene"])
genes2 <- mapIds(org.Hs.eg.db, keys = genes, column = "SYMBOL",
                 keytype = "ENSEMBL", multiVals = "first")

corr_models$name[corr_models$Type == "Gene"] <- genes2[corr_models$rowname[corr_models$Type == "Gene"]]

microorganism[, "Genus"] <- ifelse(is.na(microorganism[, "Genus"]), "", microorganism[, "Genus"])
micro <- apply(microorganism[, c("Family", "Genus")], 1, paste, collapse = " ")
names(micro) <- rownames(microorganism)
corr_models$name[corr_models$Type != "Gene"] <- micro[corr_models$rowname[corr_models$Type != "Gene"]]
corr_models$name[corr_models$name == "NA "] <- NA
saveRDS(corr_models, paste0("output/correlations_", p, "perc.RDS"))
}
{
corr_models %>%
  filter(!is.na(name)) %>%
  xlsx::write.xlsx(paste0("output/correlations_", p, "perc.xlsx"), sheetName = "Sheet1",
                 col.names = TRUE, row.names = FALSE, append = FALSE)

entrez <- mapIds(org.Hs.eg.db, keys = corr_models$rowname, keytype = "ENSEMBL", column = "ENTREZID")
g <- split(entrez, corr_models$Model)
g2 <- lapply(g, function(x){x[!is.na(x)]})
universe <- mapIds(org.Hs.eg.db, keys = trimVer(rownames(model3_best$a$RNAseq)), keytype = "ENSEMBL", column = "ENTREZID")
genes_list <- entrez[!is.na(entrez)]
universe <- universe[!is.na(universe)]
m_t2g <- msigdbr(species = "Homo sapiens", category = "C7") %>%
  dplyr::select(gs_name, entrez_gene)
m_c7 <- lapply(g2, enricher, universe = universe, TERM2GENE = m_t2g)
sapply(m_c7, function(x){nrow(as.data.frame(x))})
m_t2g_c5 <- msigdbr(species = "Homo sapiens", category = "C5") %>%
  dplyr::select(gs_name, entrez_gene)
m_c5 <- lapply(g2, enricher, universe = universe, TERM2GENE = m_t2g_c5)
sapply(m_c5, function(x){nrow(as.data.frame(x))})
m_path <- lapply(g2, enrichPathway, universe = universe)
sapply(m_path, function(x){nrow(as.data.frame(x))})
m_dose <- lapply(g2, DOSE::enrichDO,
                 ont           = "DO",
                 pvalueCutoff  = 0.05,
                 pAdjustMethod = "BH",
                 universe      = universe,
                 minGSSize     = 5,
                 maxGSSize     = 500,
                 qvalueCutoff  = 0.05,
                 readable      = FALSE)
sapply(m_dose, function(x){nrow(as.data.frame(x))})
m_MeSH <- lapply(g2, meshes::enrichMeSH,
                 MeSHDb = "MeSH.Hsa.eg.db",
                 database = 'gendoo',
                 category = 'C')
sapply(m_MeSH, function(x){nrow(as.data.frame(x))})


# Visualize the differences
upset_data <- corr_models %>%
  filter(!is.na(name)) %>%
  select(-rowname, -r, -p) %>%
  mutate(n = 1) %>%
  distinct() %>%
  pivot_wider(names_from = Model,
              values_from = n,
              values_fill = 0)
upset(upset_data, intersect = c("m0", "m1.2", "m2.2"),
      base_annotations = list(
        'Size' = (
          intersection_size(
            mode = "exclusive_intersection",
            mapping = aes(fill = Type),
            size = 0,
            text = list(check_overlap = FALSE)
          ))
      ),
      set_sizes = (
        upset_set_size(
          geom = geom_bar(
            aes(fill = Type),
          ),
        )
      ),
      # moves legends over the set sizes
      guides = 'over'
)
ggsave(paste0("Figures/upset_corr_PC_", p, "percent.png"))

}
# * Plots correlations ####
plot_cor_PC <- function(x, y, name, filter_x) {
  constx <- 3.58
  consty <- 0.69
  pdf(name)
  for (i in seq_along(x)) {
    for (j  in seq_len(ncol(x[[i]]))) {
      name <- colnames(x[[i]])[j]
      if (!trimVer(name) %in% filter_x) {
        next
      }
      ex <- x[[i]][, j] + constx
      pc <- y[[i]][, 1] + consty
      cor <- cor.test(log10(ex), log10(pc), method = "spearman")
      if (cor$p.value > 0.05) {
        next
      }
      plot(log10(ex), log10(pc), xlab = name, ylab = names(y)[i],
           main = paste0("r=", as.character(round(cor$estimate, 3)),
                         "\np.value=", as.character(signif(cor$p.value, 4))))
    }
  }
  dev.off()
}

plot_cor_PC(A[1:2], model0$Y, "output/correlations_PC_m0_log10.pdf", corr_models$rowname[corr_models$Model == "m0"])
plot_cor_PC(A[1:2], model2_best$Y, "output/correlations_PC_m1.2_log10.pdf", corr_models$rowname[corr_models$Model == "m1.2"])
plot_cor_PC(A[1:2], model3_best$Y, "output/correlations_PC_m2.2_log10.pdf", corr_models$rowname[corr_models$Model == "m2.2"])

pc01 <- cor_PC(model0$Y, model2_best$Y)
pc02 <- cor_PC(model0$Y, model3_best$Y)
pc12 <- cor_PC(model2_best$Y[1:2], model3_best$Y)

# Correlation tests ####
# corr <- corr.test(A$RNAseq, A$Micro, method = "spearman", ci = FALSE)
# saveRDS(corr, "output/correlations_genes_ASV.RDS")
# corr2 <- corr[c("r", "p")]
# saveRDS(corr2, "output/correlations_genes_pvalue_ASV.RDS")
corr2 <- readRDS("output/correlations_genes_pvalue_ASV.RDS")
r2 <- corr2$r %>%
  as.data.frame() %>%
  tibble::rownames_to_column() %>%
  pivot_longer(cols = -rowname) %>%
  rename(Microorganism = name, r = value)
p2 <- corr2$p %>%
  as.data.frame() %>%
  tibble::rownames_to_column() %>%
  pivot_longer(cols = -rowname) %>%
  rename(Microorganism = name, adj.p.value = value)


# Variables contribution ####

variables_relations <- function(rel, comp = c(1, 2)) {
  stopifnot(length(comp) == 2, is.numeric(comp), all(comp <= length(rel)))
  x <- t(rel[[comp[1]]])
  y <- t(rel[[comp[2]]])

    # Add offset and make logarithms
  xx <- log10(x[, apply(x, 2, filter), drop = FALSE] + min(x) + 0.1)
  yy <- log10(y[, apply(y, 2, filter), drop = FALSE] + min(y) + 0.1)
  cor(xx, yy)
}

# Filter empty values and those who have many values concentrated
filter <- function(z, limit = 0.5) {
  r <- base::xtfrm(z)
  any(z != 0) && length(r) >= 3 && entropy(r) > limit
}


entropy <- function(xy) {
  pt <- prop.table(table(xy))
  sum(-pt*log(pt, max(2, length(pt))))
}


filter_astar <- function(y, percent = 0.05) {
  y <- y[y != 0]
  names(sort(abs(y)))[1:(ceiling(length(y)*percent))]
}

nam <- mapIds(org.Hs.eg.db, keys = trimVer(rownames(model3_best$astar[[1]])),
       keytype = "ENSEMBL", column = "SYMBOL")
# * Model 0 ####
vc0 <- variables_contribution(A[1:2], model0)
w0 <- variables_relations(vc0)
x0 <- filter_astar(model0$astar[[1]][, 1])
y0 <- filter_astar(model0$astar[[2]][, 1])
w0b <- w0[rownames(w0) %in% x0, !colnames(w0) %in% mit]
rownames(w0b) <- nam[trimVer(rownames(w0b))]
w0b <- w0b[!is.na(rownames(w0b)),]
Heatmap(w0b, name = "Model 0", show_row_names = FALSE)
# * Model 1.2 ####
vc1.2 <- variables_contribution(A[1:2], model2_best)
w1.2 <- variables_relations(vc1.2)
x1.2 <- filter_astar(model2_best$astar[[1]][, 1])
y1.2 <- filter_astar(model2_best$astar[[2]][, 1])
w1.2b <- w1.2[rownames(w1.2) %in% x1.2, !colnames(w1.2) %in% mit]
rownames(w1.2b) <- nam[trimVer(rownames(w1.2b))]
w1.2b <- w1.2b[!is.na(rownames(w1.2b)),]
Heatmap(w1.2b, name = "Model 1.2", show_row_names = FALSE)
# * Model 2.2 ####
vc2.2 <- variables_contribution(A[1:2], model3_best)
w2.2 <- variables_relations(vc2.2)
x2.2 <- filter_astar(model3_best$astar[[1]][, 1])
y2.2 <- filter_astar(model3_best$astar[[2]][, 1])
w2.2b <- w2.2[rownames(w2.2) %in% x2.2, !colnames(w2.2) %in% mit]
rownames(w2.2b) <- nam[trimVer(rownames(w2.2b))]
w2.2b <- w2.2b[!is.na(rownames(w2.2b)),]
Heatmap(w2.2b, name = "Model 2.2", show_row_names = FALSE)

# * Transform for sharing ####

to_long <- function(x, yx) {
  x %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    pivot_longer(cols = -rowname, names_to = "Microorganism",
                 values_to = "r") %>%
    left_join(y = yx, by = c("Microorganism" = "rowname")) %>%
    select(-ASV) %>%
    dplyr::filter(Genus != "Mithocondria", abs(r) > cor_sign(126)) %>%
    rename("Gene" = "rowname") %>%
    arrange(Gene, -abs(r)) %>%
    as.data.frame()
}

to_long(w0b, microorganism) %>%
  write.xlsx(file = "output/m0_correlations_influence_log10.xlsx", sheetName = "Sheet1",
             col.names = TRUE, row.names = FALSE, append = FALSE)
to_long(w1.2b, microorganism) %>%
  write.xlsx(file = "output/m1.2_correlations_influence_log10.xlsx", sheetName = "Sheet1",
             col.names = TRUE, row.names = FALSE, append = FALSE)
to_long(w2.2b, microorganism) %>%
  write.xlsx(file = "output/m2.2_correlations_influence_log10.xlsx", sheetName = "Sheet1",
             col.names = TRUE, row.names = FALSE, append = FALSE)

variables_plots <- function(vc, w, gene_name, micro_ref) {
  w2 <- w %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    pivot_longer(cols = -rowname, names_to = "Microorganism",
                 values_to = "r") %>%
    dplyr::filter(abs(r) > cor_sign(126))
  constx <- min(vc[[1]]) + 0.1
  consty <- min(vc[[2]]) + 0.1
  for (i in seq_len(nrow(w2))) {
    micro <- w2[i, "Microorganism", drop = TRUE]
    gene <- w2[i, "rowname", drop = TRUE]
    ref <- micro_ref[micro, c("Family", "Genus")]
    ref <- ref[!is.na(ref)]
    if (length(ref) == 0) {
      ref <- micro
    } else {
      ref <- paste0(ref)
    }
    x <- log10(vc[[1]][gene, ] + constx)
    y <- log10(vc[[2]][micro, ] + consty)
    r <- cor.test(x, y)
    plot(x, y, xlab = gene_name[trimVer(gene)],
         ylab = ref,
         main = paste("r", signif(r$estimate, 3), "p.value", signif(r$p.value, 3)))
  }
}

pdf("Figures/m0_correlations_plots_log10.pdf")
variables_plots(vc0, w0[rownames(w0) %in% x0, !colnames(w0) %in% mit],
                gene_name = nam, micro_ref = microorganism)
dev.off()
pdf("Figures/m1.2_correlations_plots_log10.pdf")
variables_plots(vc1.2, w1.2[rownames(w1.2) %in% x1.2, !colnames(w1.2) %in% mit],
                gene_name = nam, micro_ref = microorganism)
dev.off()
pdf("Figures/m2.2_correlations_plots_log10.pdf")
variables_plots(vc2.2, w2.2[rownames(w2.2) %in% x2.2, !colnames(w2.2) %in% mit],
                gene_name = nam, micro_ref = microorganism)
dev.off()
