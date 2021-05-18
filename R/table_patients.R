meta <- readRDS("")
m <- group_by(meta, IBD)
m %>%
  count(SEX) %>%
  pivot_wider(id_cols = 2, names_from = IBD, values_from = n) %>%
  select(c(1, 3, 4, 2))
m %>%
  count(ileum) %>%
  pivot_wider(id_cols = 2, names_from = IBD, values_from = n) %>%
  select(c(1, 3, 4, 2))
m %>%
  count(Time) %>%
  pivot_wider(id_cols = 2, names_from = IBD, values_from = n) %>%
  select(c(1, 3, 4, 2))
m %>%
  summarise(min = min(AgeDiag), max = max(AgeDiag), median = median(AgeDiag)) %>%
  filter(IBD != "CONTROL")
m %>%
  summarize(Patients = n_distinct(Pacient_id))
pivot_wider(id_cols = 2, names_from = IBD, values_from = n) %>%
  select(c(1, 3, 4, 2))

# No active
meta %>% filter(IBD == "CD") %>%
  group_by(Ileum) %>%
  count(`CDEIS_partial` <= 4)
meta %>% filter(IBD == "UC") %>% count(`partial MAYO UC` <= 1)
m2 <- mutate(meta,
       active = case_when(
         IBD == "CONTROL" ~ "INACTIVE",
         IBD == "CD" & `CDEIS_partial` <= 4 ~ "INACTIVE",
         IBD == "UC" & `partial MAYO UC` <= 1 ~ "INACTIVE",
         TRUE ~ "ACTIVE"))
