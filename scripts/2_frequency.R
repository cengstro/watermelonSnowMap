library(tidyverse)
library(here)
library(fs)

regions <- read_csv(here("data/s2_classifier_map/training_region_name_dict.csv")) %>% 
  rename(region=name) %>% 
  clean_names() %>% 
  select(region, subregion_of)

# @ 50 m resolution
freq <- dir_ls(here("data/s2_classifier_map/freq_histograms/"), regexp="csv") %>% 
  map_df(read_csv, id = "path") %>% 
  mutate(region = basename(path) %>% str_remove(".csv"),
         .keep="unused")

parse_dict <- function(dict){
  dict %>% 
    str_remove("\\{") %>% 
    str_remove("\\}") %>% 
    str_split(", ", simplify = TRUE) %>% 
    str_split("=", simplify = TRUE) %>% 
    as_tibble() %>% 
    mutate_all(as.numeric)
}

hist_data <- freq %>% 
  janitor::clean_names() %>% 
  mutate(histogram = map(n_yrs, parse_dict), .keep="unused") %>% 
  unnest(histogram) %>% 
  rename(freq = V1, count = V2) %>% 
  arrange(freq) 
hist_data

# pre region statistics
hist_data %>% 
  group_by(region) %>% 
  mutate(frac = count/sum(count)) %>% 
  ungroup() %>% 
  select(-count) %>% 
  pivot_wider(names_from = freq, values_from = frac)

# global stats
hist_data %>% 
  group_by(freq) %>% 
  summarise(sum = sum(count)) %>% 
  filter(freq!=0) %>% 
  mutate(frac = sum/sum(sum)) 
# # A tibble: 5 × 3
# freq      sum   frac
# <dbl>    <dbl>  <dbl>
#   1     0 2896615. 0.262 
# 2     1 4664226. 0.422 
# 3     2 2291354. 0.208 
# 4     3 1030566. 0.0933
# 5     4  157355. 0.0143


# where are the perennial sites?
hist_data %>% 
  filter(freq==4) %>% 
  left_join(regions) %>% 
  group_by(subregion_of) %>% 
  summarise(count = sum(count)) %>% 
  mutate(frac = count/sum(count))
  
