library(tidyverse)
library(here)
library(kableExtra)

theme_set(theme_bw())

my_read <- function(df, string){
  df %>% 
    read_csv() %>% 
    janitor::clean_names() %>%
    select(name, subregion_of, aspect_masked, dem, slope) %>% 
    add_column(type = string)
}

alga <- my_read(here("data/s2_classifier_map/covariates/demHistAlgae.csv"), string = "alga")
glac <- my_read(here("data/s2_classifier_map/covariates/demHistGlacier.csv"), string="glac")
both <-  alga %>% bind_rows(glac)

# wrangle --------------------

parse_dict <- function(dict){
  dict %>% 
    str_remove_all("[{}]") %>% 
    str_split(", ", simplify = TRUE) %>% 
    str_split("=", simplify = TRUE) %>% 
    as_tibble() %>% 
    mutate_all(as.numeric)
}

my_tidy <- function(col){
  both %>% 
    select(type, name, subregion_of, {{col}}) %>% 
    mutate(tidy_data = map({{col}}, parse_dict), .keep="unused") %>% 
    unnest(cols = c(tidy_data)) %>% 
    drop_na()
}


aspect <- my_tidy(aspect_masked) %>% 
  rename(aspect_deg = V1, count = V2) 

slope <- my_tidy(slope) %>% 
  rename(slope = V1, count = V2) 

elev <- my_tidy(dem) %>% 
  rename(elev_m = V1, count = V2) 


# plots ---------------------------------------



## slope -------------------------------

slope %>%   
  # use the no sheet, otherwise we are going arbitrarily far into the ice sheet
  filter(!(name!="greenlandNoSheet" & subregion_of=="greenland")) %>% #distinct(name) %>% print(n=24)
  filter(!(subregion_of %in% c("caucasus", "newZealand","antarctica"))) %>% # no validated blooms in there regions
  filter(subregion_of!="arctic") %>% #very high error in this region
  ggplot(aes(x = slope, weight =count, color = type, fill = type)) +
  geom_density(alpha = 0.2) +
  facet_wrap(vars(subregion_of))

ggsave(here("figs/distribution/slope_density.png"))

slope2 <- slope %>% 
  filter(!str_detect(name, "NoSheet")) %>% # these overlap the other greenland and antarctica polygons
  # calculate fraction algal cover, per slope bin
  pivot_wider(names_from = type, values_from = count, values_fill = 0) %>% 
  group_by(name, slope) %>% 
  summarise(alga = sum(alga), glac = sum(glac)) %>% 
  mutate(frac = alga/glac, frac_wt = frac*glac) %>% 
  arrange(slope) %>% 
  ungroup()

slope2 %>% 
  ggplot(aes(x = slope, y = frac_wt, color = glac)) +
  geom_point()+
  facet_wrap(vars(name), scales = "free_y")
# could this just be that there is more area of low slope, so by chance there is more algae there?



# KS statistic
ks.test(slope2$alga, slope2$glac)

slope2 %>% 
  group_by(name) %>% 
  nest() %>% 
  mutate(ks = map(data, ~ks.test(.x$alga, .x$glac) %>% broom::tidy())) %>% 
  unnest(ks) %>% 
  arrange(-p.value)




## aspect ----------

aspect_global <- aspect %>% 
  filter(!(name!="greenlandNoSheet" & subregion_of=="greenland")) %>% #distinct(name) %>% print(n=24)
  filter(!(subregion_of %in% c("caucasus", "newZealand","antarctica"))) %>% # no validated blooms in there regions
  filter(subregion_of!="arctic") %>% #very high error in this region
  mutate(hemisphere = case_when(name %in% c("antarctica","andes","newZealand")~ "south",
                                .default = as.character("north"))) %>% 
  mutate(solar_exp = if_else(hemisphere=="north",
                             -1 * cos((pi*aspect_deg)/180),
                             cos((pi*aspect_deg)/180)) )
view(aspect_global)

aspect_global %>% 
  ggplot(aes(x = solar_exp, weight =count, color = type, fill = type)) +
  geom_density(alpha = 0.2) +
  facet_wrap(vars(name)) +
  geom_vline(xintercept = 0)
# where 1 is south and -1 is north in the northern hemisphere


        



    
  
aspect_global %>% 
  ggplot(aes(x = direction, y = frac)) +
  geom_col()+
  facet_wrap(vars(hemisphere), scales = "free")
  

aspect %>% 
  filter(!str_detect(name, "NoSheet")) %>% # these overlap the other greenland and antarctica polygons
  mutate(hemisphere = case_when(name %in% c("antarctica","andes","newZealand")~ "south",
                                .default = as.character("north")) ) %>% 
  # normalize per glacier area
  pivot_wider(names_from = type, values_from = count) %>% 
  mutate(frac = alga/glac) %>% 
  ggplot(aes(x = aspect_deg, y = frac)) +
    geom_point(alpha = 0.1)+
  lims(y=c(0,0.25)) +
  geom_smooth(method = "lm")




# elevation ---------------------

elev %>% 
  filter(!(name!="greenlandNoSheet" & subregion_of=="greenland")) %>% #distinct(name) %>% print(n=24)
  filter(!(subregion_of %in% c("caucasus", "newZealand","antarctica"))) %>% # no validated blooms in there regions
  filter(subregion_of!="arctic") %>% #very high error in this region
  filter(!(name %in% c("brooks", "highMtnAsia")) ) %>%
  # mutate(name = if_else(str_detect(name, "norway"), "norway", name)) %>% 
  ggplot(aes(x = elev_m, weight =count, color = type, fill = type)) +
  geom_density(alpha = 0.2) +
  facet_wrap(vars(name), scales="free", ncol = 3) +
  theme(legend.position="bottom")
ggsave(here("figs/distribution/elev_density.png"), width = 8, height = 8, units = "in")

# scratch------------------------------



aa
aa %>% 
  mutate(name = fct_relevel(name, my_order)) %>% 
  ggplot(aes(x = aspect_bin, weight = aspect_count, fill = type, color = type)) +
  geom_density(alpha = 0.3) +
  facet_grid(rows=vars(name), scales="free", switch="y") +
  labs(y = "density") +
  scale_y_continuous(position = "right") +
  theme_bw() +
  theme(
    strip.text.y.left = element_text(angle = 0),
    axis.title.x = element_text(),
    axis.text.x = element_text(),
    axis.ticks.x = element_line(),
    # axis.text.y = element_text(size = 5),
    # panel.grid.minor = element_blank(),
    legend.position = "none")

my_order <- c(
  # northAmerica
  "alaskaPeninsula", "alaskaRange", "brooks", "coastNorth", "coastSouth","interiorNorth","interiorSouth","cascades",
  # greenland
  "southwestGreenland", "southeastGreenland","westCentralGreenland", "eastCentralGreenland", "northGreenland","greenlandNoSheet", 
  # europe
  "iceland", "norwayNorth", "alps","caucasus",
  # arctic
  "ellesmere", "baffin", "svalbard","russianArctic",
  # centralAsia
  "highMtnAsia","altai", 
  # southern hemisphere
  "andes",
  "newZealand",
  "kamchatka",
  "antarctica","antarcticaNoSheet"
)

  




  
  # minpix <- 10 # must be at least this n pix in a region 
  
#   # scratch
  
#   filter(!name %in% c("northGreenland","northwestGreenland","northeastGreenland","eastCentralGreenland"))
# mutate(name = case_when(
#   name %in% c("northwestGreenland","northeastGreenland") ~ "northGreenland",
#   name %in% c("norwayNorth","norwaySouth") ~ "norway",
#   name=="eastCentralGreenland"~"eastGreenland",
#   name=="westCentralGreenland"~"westGreenland",
#   name %in% c("ellesmere", "baffin","svalbard","russianArctic")~"arctic",
#   .default = as.character(name) ) )






