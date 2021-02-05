# Explore the microbiome and transcripme as a result of finding two outlier
# samples from the same patient 52
library("phyloseq")
library("dada2")
library("ggplot2")
library("metacoder")
library("taxa")

# Microbiome ####
seqtab.nochim <- readRDS("data/ASV.RDS")

sample_names <- gsub("_.*", "", rownames(seqtab.nochim))
rownames(seqtab.nochim) <- sample_names
taxonomy <- readRDS("data_out/taxonomy_ASV.RDS")
taxa <- taxonomy$tax
A <- readRDS("data/RGCCA_data.RDS")
samdf <- A$Meta
rownames(samdf) <- samdf$Original
seqtab.nochim <- seqtab.nochim[rownames(samdf), ]
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
               sample_data(samdf),
               tax_table(taxa))

ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method = "NMDS", distance = "bray")
plot_ordination(ps.prop, ord.nmds.bray, color = "Pacient_id", title = "Bray NMDS")

top20 <- names(sort(taxa_sums(ps), decreasing = TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU, na.rm = TRUE))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x = "IBD", fill = "Family") +
  facet_wrap(~Time, scales = "free_x")


plot_heatmap(ps.top20, sample.label = "Pacient_id")
plot_heatmap(ps)
(p <- plot_heatmap(ps, "NMDS", "bray", "IBD", "Family"))
plot_heatmap(ps, "NMDS", "jaccard")
plot_richness(ps, measures = "Shannon")
plot_richness(ps, "IBD", measures = "Shannon", col = "Pacient_id")
dm <- estimate_richness(ps, measures = "Shannon")
pst <- merge_samples(ps, "IBD")
plot_richness(pst, measures = "Shannon")

# https://grunwaldlab.github.io/analysis_of_microbiome_community_data_in_r/07--diversity_stats.html
# ps %>%
#   tax_table() %>%
#   taxa::filter_taxa(leaf > 50) %>% # taxa:: needed because of phyloseq::filter_taxa
#   heat_tree(node_label = taxon_names,
#             node_size = leaf,
#             node_color = leaf,
#             layout = "da", initial_layout = "re",
#             title = "Taxa in leafs")


# Transcriptome ####
A <- readRDS("data/RGCCA_data_server.RDS")
rna <- A$RNAseq
d <- dist(rna)
p <- cmdscale(d)
plot(p)
points(p[startsWith(rownames(p), "052"), ], pch = 15)

## Boxplot ####
png("Figures/boxplot_rnaseq.png")
boxplot(t(rna), main = "Boxplot of normalized values",
        border = ifelse(startsWith(rownames(p), "052"), "red", "black"))
dev.off()
## Hierarchical cluster ####
png("Figures/hierarchical_rnaseq.png")
h <- hclust(d)
dendo <- as.dendrogram(h)
col <- ifelse(startsWith(labels(dendo), "052"), "red", "black")
names(col) <- labels(dendo)
dendextend::labels_colors(dendo) <- col
plot(dendo)
dev.off()
