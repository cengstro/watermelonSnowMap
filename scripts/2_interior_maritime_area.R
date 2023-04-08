library(tidyverse)
library(here)
library(janitor)
library(kableExtra)
library(snakecase)

areas <- read_csv(here("data/s2_classifier_map/area/interiorVsMaritime.csv"))


areaclean <- areas %>% 
  clean_names() %>%
  rename(sum_rgnd = mean_rgnd,
         cloud_mean = mean, 
         cloud_sd = std_dev,
         cloud_count = count) %>% 
  select(name, algae_km2, glacier_km2, sum_rgnd, contains("cloud")) %>% 
  mutate(algae_frac_cover = algae_km2/glacier_km2, .before = sum_rgnd) %>% 
  add_column(side = c("Interior","Maritime","Interior","Maritime","Interior","Maritime"),
             .before  = 2) %>% 
  mutate(name = name %>% str_remove("Interior|Maritime|North|South")) %>% 
  mutate(name = case_when(name=="bcCoast"~"coastSouth",
                          name=="alaskaCoast"~"coastNorth",
                          name=="denali"~"denaliHayesDelta")) %>% 
  rename(region=name)
  
areaclean 
# in Alaska the south:north algae frac cover ratio is 2x

round100 <- function(x){
  case_when(x>10000 ~ round(x,-3),
            x>1000 ~ round(x,-2),
            x>100 ~ round(x,-1),
            .default = round(x, 0))
}


areaclean %>% 
  mutate(across(algae_km2:glacier_km2, ~round100(.x))) %>% 
  mutate(algae_frac_cover = round(algae_frac_cover*100, 0)) %>% 
  mutate(across(c(cloud_mean, cloud_sd, sum_rgnd), ~round(.x,digits=1))) %>% 
  mutate(region =  to_title_case(as.character(region))) %>% 
  mutate(region = if_else(region=="Denali Hayes Delta", "Denali+Hayes+Delta", region)) %>% 
  kable(col.names = c("Region", "Side","Algae $\\text{km}^2$", "Glacier $\\text{km}^2$", 
                      "Algae % cover","Mean. Ann. Biomass (mean RGND)", "Cloud Pix Count", 
                      "Cloud Freq Mean", "Cloud Freq SD"))  %>% 
  kable_classic("striped", full_width = F, html_font = "Cambria")



# anova comparison of cloud cover

areaclean %>% 

