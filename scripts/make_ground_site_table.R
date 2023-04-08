library(tidyverse)
library(officer)
library(flextable)
library(here)

tt <- read_csv(here("data/field_sites.csv"))

tt2 <- tt %>% 
  mutate(lon = round(lon, digits = 4),
         lat = round(lat, digits = 5)) %>% 
  arrange(visit_date) %>% 
  mutate(n_samples  = as.character(n_samples)) %>% 
  replace_na(list(n_samples = "-"))


colnames(tt2) <-c("Site name", "S2 aquisition date", "Visit date","Longitude", "Latitude", 
                      "N. samples", "Notes")

tt2 %>% 
  flextable() %>%
  save_as_docx( path = here("figs/ground_truth/groundtruth_tbl.docx"))
