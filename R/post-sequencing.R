library("readxl")
library("tidyverse")

# Dealing with the names ####
Names <- "data/dataDelivered/amplicon_names.xlsx"
names1 <- read_xlsx(Names, range = anchored("A4", dim = c(9, 13)))
names2 <- read_xlsx(Names, range = anchored("A21", dim = c(9, 13)))
names3 <- read_xlsx(Names, range = anchored("A38", dim = c(9, 13)))
names4 <- read_xlsx(Names, range = anchored("A55", dim = c(9, 13)))

n1 <- names1 %>%
  gather("Column", "Name", 2:13, convert = TRUE) %>%
  rename(Row = `X__1`) %>%
  mutate(Plate = 1)
n2 <- names2 %>%
  gather("Column", "Name", 2:13, convert = TRUE) %>%
  rename(Row = `X__1`) %>%
  mutate(Plate = 2)
n3 <- names3 %>%
  gather("Column", "Name", 2:13, convert = TRUE) %>%
  rename(Row = `X__1`) %>%
  mutate(Plate = 3)
n4 <- names4 %>%
  gather("Column", "Name", 2:13, convert = TRUE) %>%
  rename(Row = `X__1`) %>%
  mutate(Plate = 4)

nam <- rbind(n1, n2, n3, n4)

# Dealing with the Concentrations ####
Concentrations <- "data/dataDelivered/concentr_FinalProject.xlsx"
plate1 <- read_xlsx(Concentrations, range = anchored("B4", dim = c(9, 13)),
                    skip = 4, n_max = 8)
p1 <- plate1 %>%
  gather("Column", "Concentr", 2:13, convert = TRUE) %>%
  rename(Row = `X__1`) %>%
  mutate(Plate = 1, Concentr = as.numeric(Concentr))

plate2 <- read_xlsx(Concentrations, range = anchored("B22", dim = c(9, 13)),
                    skip = 4, n_max = 8)
p2 <- plate2 %>%
  gather("Column", "Concentr", 2:13, convert = TRUE) %>%
  rename(Row = `X__1`) %>%
  mutate(Plate = 2, Concentr = as.numeric(Concentr))

plate3 <- read_xlsx(Concentrations, range = anchored("B40", dim = c(9, 13)),
                    skip = 4, n_max = 8)
p3 <- plate3 %>%
  gather("Column", "Concentr", 2:13, convert = TRUE) %>%
  rename(Row = `X__1`) %>%
  mutate(Plate = 3, Concentr = as.numeric(Concentr))

plate4 <- read_xlsx(Concentrations, range = anchored("B58", dim = c(9, 13)),
                    skip = 4, n_max = 8)
p4 <- plate4 %>%
  gather("Column", "Concentr", 2:13, convert = TRUE) %>%
  rename(Row = `X__1`) %>%
  mutate(Plate = as.character(4), Concentr = as.numeric(Concentr))

plates <- rbind(p1, p2, p3, p4)

theme_set(theme_bw())
plates <- plates %>%
  mutate(Plate = as.factor(Plate), Row = as.factor(Row),
         Column = as.factor(Column))

ggplot(plates) +
  geom_violin(aes(Plate, Concentr))

ggplot(plates) +
  geom_boxplot(aes(Plate, Concentr)) +
  facet_wrap(~Row)

ggplot(plates) +
  geom_boxplot(aes(Plate, Concentr)) +
  facet_wrap(~Column)

library("experDesign")
i <- split(seq_along(plates$Plate), plates$Plate)
res <- evaluate_index(i, plates[, "Concentr", drop = FALSE])
res[, "Concentr", ]


# Merging with previously used data ####
out <- as.tibble(merge(nam, plates))
spli <- strsplit(out$Name, "_")
out$root <- sapply(spli, function(x){x[1]})
out$duplicated <- sapply(spli, function(x){ifelse(length(x) > 1, x[2], NA)})
out$Plate <- as.factor(out$Plate)
out$Column <- as.factor(out$Column)
out$Row <- as.factor(out$Row)

# Read each sheet
namesSamples <- lapply(1:4, function(x) {
  read_xlsx("data/dataDelivered/michigan_samples_modified.xlsx",
            sheet = x)
})

namesSamples <- do.call(rbind, namesSamples)
# drop columns
namesSamples <- namesSamples[, -c(9:11, 13:15)]
m <- merge(out, namesSamples, by.x = "root", by.y = "X__1",
           all.x = TRUE, all.y = FALSE)
# Discard the duplicated values
mD <- m[is.na(m$duplicated), ]
# The distribution seems quite ok, entropy is >= 0.65
i <- split(seq_along(mD$Plate), mD$Plate)
res <- evaluate_index(i, mD[, c("IBD", "SEX", "Study", "Exact_location"), drop = FALSE])
apply(res[c("na", "entropy"), , ], c(1,2), sum)
