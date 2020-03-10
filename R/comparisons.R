# Comparisons a grosso modo
library("limma")
library("dplyr")
library("ggplot2")

# Start ####
A <- readRDS("data/RGCCA_data.RDS")
rna <- t(A$RNAseq)
meta <- A$Meta
meta$Ileum <- ifelse(meta$Exact_location == "ileum", "ileum", "colon")
levels(meta$IBD) <- c("CD", "UC", "C")
meta$IBD[is.na(meta$IBD)] <- "C"
treatment <- ifelse(meta$Time == "0" | is.na(meta$Time), "No", "Yes")

meta %>%
  group_by(IBD) %>%
  summarise(mean = mean(Age), median = median(Age), sd = sd(Age),
            min = min(Age), max = max(Age))

ggplot(meta) +
  geom_jitter(aes(Age, IBD, col = IBD, shape = IBD)) +
  theme_minimal() +
  labs(title = "Age of the samples", size = "Samples",
       col = "Status",
       shape = "Status", y = element_blank()) +
  scale_size(range = c(1, 3))
ggsave("Figures/Age_status.png")
meta %>%
  mutate(Ileum = case_when(Exact_location == "ileum" ~ "Ileum",
                           TRUE ~ "Colon")) %>%
  ggplot() +
  geom_jitter(aes(Age, Ileum, col = Ileum, shape = Ileum)) +
  theme_minimal() +
  labs(title = "Age of the samples", size = "Samples",
       col = "Status",
       shape = "Status", y = element_blank()) +
  scale_size(range = c(1, 3))
ggsave("Figures/Age_location.png")

meta %>%
  filter(IBD != "C") %>%
  ggplot() +
  geom_jitter(aes(diagTime, IBD, col = IBD, shape = IBD)) +
  theme_minimal() +
  labs(title = "Years since diagnostic", size = "Samples",
       col = "Status", x = element_blank(),
       shape = "Status", y = element_blank()) +
  scale_size(range = c(1, 3))
ggsave("Figures/years_diagnostic.png")

meta %>%
  filter(IBD != "C") %>%
  mutate(Treatment = case_when(Time == "0" ~ "No",
                               TRUE ~ "Yes")) %>%
  count(Treatment)

# model ####
grouping <- as.factor(paste(meta$IBD, meta$Ileum, sep = "_"))

design <- model.matrix(~0 + grouping, data = grouping)
stopifnot(all(colSums(design) > 1))
colnames(design) <- levels(grouping)
design <- cbind(Intercept = 1, design)

# v <- voom(rna, design = design, plot = TRUE)
# The data I am using has already been normalized by voom
# The data had some batch effect detected by looking at the PUR
# Percentage unique reads, which was corrected with everything together...

# comparisons ####
# Location ileum vs colon
# Disease C vs CD
# Treatment No vs Yes
contr <- makeContrasts(ileum_vs_colon = CD_ileum + C_ileum - C_colon - CD_colon - UC_colon,
                       CD_vs_C = + CD_ileum + CD_colon - C_ileum - C_colon,
                       UC_vs_C = UC_colon - C_colon - C_ileum,
                       UC_vs_CD = UC_colon  - CD_colon - CD_ileum,
                       levels = colnames(design))
contr
pick <- design %*% contr
pick[pick == 0] <- NA
pick[pick == -1] <- 2
stopifnot(all(colSums(pick != 0) > 3))

# Fits ####
l <- vector("list", ncol(pick))
for (i in seq_len(ncol(pick))) {
  fit <- lmFit(rna, pick[, i, drop = FALSE])
  fit2 <- eBayes(fit)
  l[[i]] <- topTreat(fit2, coef = 1, number = Inf)
}

fit2 <- contrasts.fit(fit, contr)
results <- decideTests(fit2, lfc = log2(1.5), adjust.method = "BH", p.value = 0.05)
summary(results)
vennDiagram(results)
res <- results
res[res == -1] <- 1
UpSetR::upset(as.data.frame(res), order.by = "freq", nintersects = NA)

v <- vector("list", ncol(contr))
for (i in seq_len(ncol(contr))) {
  tt <- topTreat(fit2, coef = i, number = Inf)
  v[[i]] <- tt
}
names(v) <- colnames(contr)

unlog <- function(x){
  sign(x)*2^abs(x)
}
