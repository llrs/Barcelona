library("forcats")
library("ggplot2")
library("dplyr")
library("tidyr")
library("scales")

family <- read.delim("data/Partek_Michigan3_Kraken_Classified_family.tsv", check.names = FALSE)


family.tidy <- gather(family, Sample, Count, -'Sample name') %>%
  filter(Count != 0) %>%
  rename(Microorganism = 'Sample name') %>%
  group_by(Sample) %>%
  mutate(ratio = Count/sum(Count)) %>%
  mutate(
    Study = case_when(
      grepl("-w0", Sample) ~ "BCN",
      grepl("^500", Sample) ~ "NC",
      grepl("^C", Sample) ~ "Controls",
      grepl("-T", Sample) ~ "TRIM"
    )
  ) %>%
  ungroup %>%
  mutate(Sample = gsub("_S.+", "", Sample)) %>%
  separate(Sample, c("Original", "Rep"), sep = "_p", convert = TRUE,
           fill = "right", remove = FALSE) %>%
  mutate(Rep = gsub("(^[0-9]).+", "\\1", Rep)) %>%
  mutate_if(is.character, as.factor)

saveRDS(family.tidy, "Samples_counts.RDS")

family.tidy %>%
  filter(!is.na(Rep)) %>%
  ggplot() +
  geom_col(aes(Sample, ratio, fill = Microorganism)) +
  guides(fill = FALSE) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()) +
  facet_grid(~Rep)

family.tidy %>%
  filter(!is.na(Rep)) %>%
  ggplot(aes(group = Sample, Count, fill = Microorganism)) +
  geom_bar()


BCN <- grep("-w0", colnames(family), value = TRUE)
NC <- grep("^500", colnames(family), value = TRUE)
TRIM <- grep("-T", colnames(family), value = TRUE)
C <- grep("^C", colnames(family), value = TRUE)
TRIM <- TRIM[!TRIM %in% C]
order_names <- c(BCN, NC, TRIM, C)


family.tidy$Sample <-  fct_relevel(family.tidy$Sample, order_names)
family.tidy$Study <-  fct_relevel(family.tidy$Study, c("BCN", "NC", "TRIM", "Controls"))
a <- theme_bw() +
  theme(axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      panel.grid.major.x = element_blank())
theme_set(a)

family.tidy %>%
  group_by(Sample) %>%
  summarise(counts = sum(Count),
            study = unique(Study)) %>%
  ungroup() %>%
  droplevels() %>%
  ggplot() +
  geom_col(aes(lvls_reorder(Sample, order(counts)), log10(counts), fill = study)) +
  labs(title = "BCN", x = "Samples", fill = "Study") +
  geom_hline(yintercept = c(log10(1000), log10(20000)), col = c("red", "green"))

family.tidy %>%
  filter(Study == "BCN") %>%
  ggplot() +
  geom_col(aes(Sample, log10(Count), fill = Microorganism)) +
  guides(fill = FALSE) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()) +
  labs(title = "BCN")
bcn <- family.tidy %>%
  filter(Study == "BCN") %>%
  ggplot() +
  geom_col(aes(Sample, ratio, fill = Microorganism)) +
  guides(fill = FALSE) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()) +
  labs(title = "BCN")
nc <- family.tidy %>%
  filter(Study == "NC") %>%
  ggplot() +
  geom_col(aes(Sample, ratio, fill = Microorganism)) +
  guides(fill = FALSE) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()) +
  labs(title = "Water")
trim <- family.tidy %>%
  filter(Study == "TRIM") %>%
  ggplot() +
  geom_col(aes(Sample, ratio, fill = Microorganism)) +
  guides(fill = FALSE) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()) +
  labs(title = "TRIM")
Controls <- family.tidy %>%
  filter(Study == "Controls") %>%
  ggplot() +
  geom_col(aes(Sample, ratio, fill = Microorganism)) +
  guides(fill = FALSE) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()) +
  labs(title = "Controls", caption = "Michigan")

library("patchwork")
(bcn + nc + trim)/Controls
(bcn + trim)/(nc + Controls)
