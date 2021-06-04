library("readxl")
library("tidyverse")
library("experDesign")

# Dealing with the names ####
Names <- "data/dataDelivered/amplicon_names.xlsx"
names1 <- read_xlsx(Names, range = anchored("A4", dim = c(9, 13)))
names2 <- read_xlsx(Names, range = anchored("A21", dim = c(9, 13)))
names3 <- read_xlsx(Names, range = anchored("A38", dim = c(9, 13)))
names4 <- read_xlsx(Names, range = anchored("A55", dim = c(9, 13)))

plating <- function(x, plate) {
  x %>%
    gather("Column", "Name", 2:13, convert = TRUE) %>%
    rename(Row = `...1`) %>%
    mutate(Plate = !!plate)
}
n1 <- plating(names1, "1")
n2 <- plating(names2, "2")
n3 <- plating(names3, "3")
n4 <- plating(names4, "4")

nam <- rbind(n1, n2, n3, n4)

# Dealing with the Concentrations ####
Concentrations <- "data/dataDelivered/concentr_FinalProject.xlsx"
plate1 <- read_xlsx(Concentrations, range = anchored("B4", dim = c(9, 13)),
                    skip = 4, n_max = 8)
plate2 <- read_xlsx(Concentrations, range = anchored("B22", dim = c(9, 13)),
                    skip = 4, n_max = 8)
plate3 <- read_xlsx(Concentrations, range = anchored("B40", dim = c(9, 13)),
                    skip = 4, n_max = 8)
plate4 <- read_xlsx(Concentrations, range = anchored("B58", dim = c(9, 13)),
                    skip = 4, n_max = 8)

concentrate <- function(x, plate) {
  x %>%
    gather("Column", "Concentr", 2:13, convert = TRUE) %>%
    rename(Row = `...1`) %>%
    mutate(Plate = !!plate, Concentr = as.numeric(Concentr))
}
# Warnings expected
p1 <- concentrate(plate1, "1")
p2 <- concentrate(plate2, "2")
p3 <- concentrate(plate3, "3")
p4 <- concentrate(plate4, "4")

plates <- rbind(p1, p2, p3, p4)


# Merging with previously used data ####
out <- merge(nam, plates) %>%
  as.tibble() %>%
  filter(!is.na(Name)) %>%
  mutate(Plate = as.factor(Plate), Row = as.factor(Row),
         Column = as.factor(Column))
# Because of duplicated names we need to tweak the names
any(duplicated(out$Name))
index <- grep("^22-T52-T-DM", out$Name)
out$Name[index] <- paste0(out$Name[index], "_p", out$Plate[index])
saveRDS(out, "output/Samples_concentration_distribution.RDS")

theme_set(theme_bw())

ggplot(out) +
  geom_violin(aes(Plate, Concentr))

ggplot(out) +
  geom_boxplot(aes(Plate, Concentr)) +
  facet_wrap(~Row)

ggplot(out) +
  geom_boxplot(aes(Plate, Concentr)) +
  facet_wrap(~Column)


i <- use_index(out$Plate)
res <- evaluate_index(i, out[, "Concentr", drop = FALSE])
res[, "Concentr", ]

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
# Compare with the other variables
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
