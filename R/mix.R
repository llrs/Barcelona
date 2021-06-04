# Mix the counts and the metadata about the plate
# Compare if the different replicates how do they evolve on each plate
library("dplyr")
library("tidyr")
library("ggplot2")

theme_set(theme_bw())

location <- readRDS("output/Samples_concentration_distribution.RDS")
location$Name <- gsub("-TTR-", "-T-DM-", location$Name)
counts <- readRDS("Samples_counts.RDS")
counts$Sample <- gsub("-TTR-", "-T-DM-", counts$Sample)
counts$Original <- gsub("-TTR-", "-T-DM-", counts$Original)

info <- merge(counts, location, by.y = "Name", by.x = "Sample", all.x = TRUE)
info2 <- info %>%
  group_by(Sample) %>%
  summarise(Count = sum(Count),
            Concen = unique(Concentr),
            Row = unique(Row),
            Column = unique(Column),
            Plate = unique(Plate),
            Original = unique(Original)
            ) %>%
  ungroup() %>%
  filter(!is.na(Plate))

info2 %>%
  ggplot() +
  geom_point(aes(Count, Concen))
info2 %>%
  ggplot() +
  geom_point(aes(Count, Concen)) +
  scale_x_log10() +
  scale_y_log10()
# No relationship between higher concentration and Counts

info2 %>%
  ggplot() +
  geom_point(aes(Count, Concen)) +
  facet_wrap(~Plate) +
  scale_x_log10() +
  scale_y_log10()
# Without differences in by plate
info2 %>%
  ggplot() +
  geom_point(aes(Count, Concen, col = as.factor(Plate))) +
  facet_grid(Column~Row) +
  scale_x_log10() +
  scale_y_log10()

info2 %>%
  ggplot() +
  geom_point(aes(Count, Concen)) +
  facet_grid(Column~Plate) +
  scale_x_log10() +
  scale_y_log10()
info2 %>%
  ggplot() +
  geom_point(aes(Count, Concen)) +
  facet_grid(Column~Plate) +
  scale_x_log10() +
  scale_y_log10()

info2 %>%
  filter(Plate == 1) %>%
  ggplot() +
  geom_point(aes(Count, Concen)) +
  facet_grid(Column~Row) +
  scale_x_log10() +
  scale_y_log10()

info2 %>%
  filter(Plate == 1) %>%
  ggplot() +
  geom_point(aes(Count, Concen)) +
  facet_grid(Column~Row)

info2 %>%
  filter(as.character(Sample) != as.character(Original)) %>%
  ggplot() +
  geom_point(aes(Count, Concen, col = Original))
# Equal concentration, but different counts!!

df <- info2 %>%
  filter(as.character(Sample) != as.character(Original))
keepOriginals <- as.character(unique(df$Original))
info2 %>%
  filter(Original %in% keepOriginals) %>%
  mutate(Plate = as.factor(Plate)) %>%
  filter(!is.na(Plate)) %>%
  ggplot() +
  geom_point(aes(Plate, Count, group = Original, col = Original)) +
  facet_wrap(~Original) +
  scale_y_log10()
# Sample  010-w000 in the 4rth plaque is wrong
# But sample C6-TTR-CA is increasing on the later plates
# Also sample 116-w014 and sample 092-46 has some weird pattern
# Some of them are bettern seen without scale transformation


info2 %>%
  filter(as.character(Sample) != as.character(Original)) %>%
  group_by(Original) %>%
  summarise(Count_sd = sd(Count),
            Count_IQR = IQR(Count),
            Count_mad = mad(Count),
            Count_max = max(Count),
            Count_min = min(Count),
            ) %>%
  mutate(Range = Count_max-Count_min) %>%
  ggplot() +
  geom_point(aes(Original, Range))

info2 %>%
  group_by(Plate) %>%
  summarise(Count_mean = mean(Count),
            Count_range = max(Count)-min(Count),
            Count_sd = sd(Count),
            Count_IQR = IQR(Count),
            Count_mad = mad(Count),
            Count_max = max(Count),
            Count_min = min(Count)) %>%
  ggplot(aes(x = Plate)) +
  geom_col(aes(y = Count_IQR, )) +
  geom_point(aes(y = Count_min)) +
  # geom_point(aes(y = Count_max))
  geom_point(aes(y = Count_mean))

# Without waters
info %>%
  mutate(Plate = as.factor(Plate)) %>%
  filter(!is.na(Plate)) %>%  # Something is wrong with the merge
  filter(!grepl("^500", Original)) %>%
  ggplot() +
  geom_boxplot(aes(Plate, Count)) +
  scale_y_log10()

# With waters
info %>%
  mutate(Plate = as.factor(Plate)) %>%
  filter(!is.na(Plate)) %>%  # Something is wrong with the merge
  ggplot() +
  geom_boxplot(aes(Plate, Count)) +
  scale_y_log10()


info %>%
  mutate(Plate = as.factor(Plate)) %>%
  filter(!is.na(Plate)) %>%
  group_by(Plate) %>%
  summarize(mean(Count), median(Count))
# The plate 4 has lower counts
info %>%
  mutate(Plate = as.factor(Plate)) %>%
  filter(!is.na(Plate)) %>%
  filter(!grepl("^500", Original)) %>%
  group_by(Plate) %>%
  summarize(mean(Count), median(Count))
# Excluding water increases the mean number of counts
info %>%
  mutate(Plate = as.factor(Plate)) %>%
  filter(!is.na(Plate)) %>%
  group_by(Plate) %>%
  summarize(n_distinct(Microorganism))
# Similar number of microorganisms per plate ~250 from a total of 265 families different

