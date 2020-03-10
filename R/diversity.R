# Calculate the alpha diversity
library("tidyr")
library("ggplot2")
library("vegan")

A <- readRDS("data/RGCCA_data.RDS")
otus <- t(A$Micro)
genus <- read.csv("data/genus.csv", row.names = 1)

# Alpha diversity ####
alpha <- prop.table(otus, 2)*100
a <- as.data.frame(alpha)
a$otus <- genus[, 1]
a <- pivot_longer(a, colnames(alpha))

ggplot(a) +
  geom_col(aes(name, value, fill = otus), col = "transparent") +
  guides(fill=FALSE) +
  theme_minimal()

# Beta diversity ####
beta <- vegdist(otus, method = "jaccard")
cmd <- cmdscale(d = beta)
plot(cmd)
