#*
#* Confidence intervals?
#*
#*
#*

library(tidyverse)
library(here)
library(flextable)
library(officer)

# algaeMask = algaeMap.select('CLASS').mosaic()
algal_pix_per_glacier <- read_csv(here("data/s2_classifier_map/algae_per_glacier/algae_pix_count_per_glacier.csv"))
glacier_pix_per_glacier <- read_csv(here("data/s2_classifier_map/algae_per_glacier/glacier_pix_count.csv"))

parse_dict <- function(dict){
  dict %>% 
    str_remove_all("[{}]") %>% 
    str_split(", ", simplify = TRUE) %>% 
    str_split("=", simplify = TRUE) %>% 
    as_tibble() %>% 
    mutate_all(as.numeric)
}
  
# the number of algae pixels per glacier (at 50 m resolution)
apix <- algal_pix_per_glacier %>% 
   janitor::clean_names() %>% 
   select(name, subregion_of, histogram) %>% 
   mutate(histogram = map(histogram, parse_dict)) %>% 
   unnest(histogram) %>% 
   rename(glac_id = V1, n_alg_pix = V2) %>% 
   filter(n_alg_pix>=5)# at least 3 algae pixels to remove some noise

# the number of glacier pixels per glacier (at 50 m resolution)
gpix <- glacier_pix_per_glacier %>% 
   janitor::clean_names() %>% 
   select(name, subregion_of, histogram) %>% 
   mutate(histogram = map(histogram, parse_dict)) %>% 
   unnest(histogram) %>% 
   rename(glac_id = V1, n_glac_pix = V2)

# get number of glaciers per region (excluding <1 km2)
glac_count_subregions <- count(gpix, name) %>% filter(name!="highMtnAsia") 
glac_count_regions <- count(gpix, subregion_of) %>% rename(name = subregion_of)
glac_count_world <- count(gpix) %>% mutate(name = "world")
n_glaciers_per_region <- glac_count_subregions %>%
  bind_rows(glac_count_regions) %>%
  bind_rows(glac_count_world) %>% 
  distinct() %>% 
  rename(n_glaciers = n)

# compute algal cover per glacier
# filter out observations where algae cover <1% of glacier to remove some noise
 algae_frac_cover <- apix %>% 
   left_join(gpix) %>% 
   mutate(frac = n_alg_pix / n_glac_pix) %>% 
   filter(frac>0.01) %>% 
   arrange(-frac)
 algae_frac_cover
 
 # export list of glac_id with high algal cover in North America (collect snowmelt time series for these)
 algae_frac_cover %>% 
   filter(frac>0.1, subregion_of=="northAmerica") %>% #cn196
   distinct(glac_id, frac) %>% 
   mutate(northing = str_sub(glac_id, -5, -1),
          easting = glac_id %>% str_remove(northing),
          glac_id = paste0("G", easting, "E", northing,"N")) %>% 
   select(glac_id, frac) %>% 
   write_csv(here("data/s2_classifier_map/algae_per_glacier/gids_north_america_gt10.csv"))
# 2128 glaciers in North America with frac> 10%
 
 algae_frac_cover %>% 
   arrange(-frac)

# bin by decile
algae_frac_decile_subregions <- algae_frac_cover %>% 
   drop_na() %>% 
   mutate(decile_bin = floor(frac*10)) %>% 
   group_by(name, decile_bin, subregion_of) %>% 
   summarise(n = n()) %>%
   arrange(-decile_bin) %>% 
  # include the higher deciles in each decile bin (eg the >40% bin should include glaciers with >50%)
   group_by(name) %>% # ensure that bins are not grouped!
   mutate(n_gt = cumsum(n)) %>% # ie the number of glaciers with AT LEAST this much algal cover
   ungroup() %>% #filter(name=="coastNorth") %>% print
   arrange(name, decile_bin) %>% 
  select(subregion_of, name, decile_bin, n, n_gt) 
algae_frac_decile_subregions #check that csum is the sum of n

# lump for global stats
algae_frac_decile_world <- algae_frac_decile_subregions %>% 
  group_by(decile_bin) %>% 
  summarise(n_gt = sum(n_gt)) %>% 
  mutate(name = "world", subregion_of = "world")

# lump for regional stats
algae_frac_decile_continents <- algae_frac_decile_subregions %>%
  filter(subregion_of !=name) %>% 
  group_by(subregion_of, decile_bin) %>% 
  summarise(n_gt = sum(n_gt)) %>% 
  ungroup() %>% 
  mutate(name = subregion_of) 

algae_frac_decile <- algae_frac_decile_subregions %>% 
  select(subregion_of, name, decile_bin, n_gt) %>% 
  bind_rows(algae_frac_decile_continents) %>% 
  bind_rows(algae_frac_decile_world) %>% 
  mutate(type = if_else(name==subregion_of, "continent", "subrange")) %>% 
  select(-subregion_of) %>% 
  filter(! (type=="subrange" & str_detect(name, "Greenland")) )



# display the results in a table
ord <- c("northAmerica","alaskaPeninsula","alaskaRange","coastNorth","coastSouth","interiorNorth","interiorSouth","cascades",        
         "greenland","greenlandNoSheet", 
         "europe","iceland","norwayNorth","norwaySouth","alps",           
         "arctic","ellesmere","baffin","svalbard","russianArctic",  
         "highMtnAsia","altai","caucasus","kamchatka",      
         "andes", "antarctica","newZealand",
         "world" )

table_s4 <- algae_frac_decile %>% 
  mutate(decile_bin = paste0(">", decile_bin*10, "%"),
         # because we filtered out the 1 %, update decile bin label
         decile_bin = if_else(decile_bin == ">0%", ">1%", decile_bin)) %>% 
  pivot_wider(names_from = decile_bin, values_from = n_gt, values_fill = 0) %>% 
  left_join(n_glaciers_per_region) %>% 
  mutate(region = fct_relevel(name, ord), .keep="unused", .before=1) %>% 
  arrange(region) %>% 
  relocate(n_glaciers, .after = region) %>% 
  mutate(region =  to_title_case(as.character(region))) %>% 
  mutate(region = if_else(type=="subrange", paste("- ",region), region)) %>% 
  select(-type) 

colnames(table_s4) <- c("Region","N. glaciers", ">1%", ">10%", ">20%", ">30%", ">40%", ">50%", ">60%", ">70%")
table_s4 %>% 
  flextable() %>% 
  colformat_double(big.mark = ",", digits = 0, na_str = "-") %>% 
  save_as_docx(path = here("figs/distribution/glacier_table_s4.docx"))


 