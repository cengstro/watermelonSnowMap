# algae % glacier coverage and biomass per region

# create algal range maps for figure2 
# calculate per glacier algal cover stats for select regions


# UPSTREAM
# 06_export_low_res_map.js

# DOWNSTREAM

# IMPORTS
# world_algae_percent.kml
# north_america_10km_grid_03thresh_nab2.kml
# gris_10km_grid_nab2.kml
# algaePerGlacier_*.csv
# algae_percent_region.csv

# EXPORTS:
# global_dist_map.png (fig 2a), 
# north_america_map_2.png (fig 2b) 
# greenland_map_2.png (fig 2c)
# per_glacier_table.png
# glac_ids_gt_10.csv
# glac_ids_gt_1.csv



library(tidyverse)
library(here)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(RColorBrewer)
library(ggspatial) # scalebar
# library(kableExtra)
# library(snakecase)


# import data -----------------

algaeworld_raw <- st_read(here("data/s2_classifier_map/area/world_algae_map50.kml") ) 
northamerica_raw <- st_read(here("data/s2_classifier_map/area/north_america_algae_map10.kml") ) 
greenland_raw <- st_read(here("data/s2_classifier_map/area/greenland_algae_map10.kml") ) 
hma_raw <- st_read(here("data/s2_classifier_map/area/hmaAlgaeMap10.kml"))

# background map data
world <- ne_countries(scale = "medium", returnclass = "sf")
states <- rnaturalearth::ne_states(c("united states of america", "canada"), returnclass = "sf")
geo_lines <- sf::st_read(here("figs/distribution/3_area/ne_110m_geographic_lines/ne_110m_geographic_lines.shx"))


# tidy  -----------------------------
wag_proj_string <-"+proj=wag4 +lon_0=0 +datum=WGS84 +units=m +no_defs" # Wagner IV proj

bins <- c("0", "1–10","10–20","20–30","30–40", ">40")

# loosly based on 9 class brewer YlOrRd
# , "#ffffcc",
my_pal <- c("#FFFFFF", "#fed976", "#fd8d3c", "#fc4e2a", "#bd0026", "#800026") #"#810f7c" #  "#bd0085"

length(bins)==length(my_pal)

bin_percent <- function(sf){
  sf %>% 
    mutate(algae_km2 = as.numeric(algae_km2),
           glacier_km2 = as.numeric(glacier_km2),
           fraction = algae_km2/glacier_km2,
           percent = fraction*100,
           percent_binned = case_when(fraction < 0.01 ~ bins[1],
                                      # (fraction >= 0.001 & fraction <= 0.01) ~ bins[2], # I pre-rounded the data to nearest percent, so use <= 
                                      (fraction >= 0.01 & fraction < 0.1) ~ bins[2],
                                      (fraction >= 0.1 & fraction < 0.2) ~ bins[3],
                                      (fraction >= 0.2 & fraction < 0.3) ~ bins[4],
                                      (fraction >= 0.3 & fraction < 0.4) ~ bins[5],
                                      # (fraction >= 0.4 & fraction < 0.5) ~ bins[6],
                                      fraction >= 0.4 ~ bins[6] ), 
           percent_binned = percent_binned %>% fct_relevel(bins) )
}

algaeworld <- algaeworld_raw %>% 
  bin_percent() %>% 
  st_transform(wag_proj_string) %>% 
  st_centroid(algaeworld) %>% 
  arrange(-fraction) # so that higher fraction points are plotted on top

max(algaeworld$fraction)

# just arctic and antarctic circles
polar_circles <- geo_lines %>% 
  sf::st_transform(wag_proj_string) %>% 
  as.data.frame() %>% 
  filter(name %>% str_detect("rctic")) %>% 
  st_as_sf()

ggplot() + 
  geom_sf(data = polar_circles, color="gray30", linetype = "dashed", size = 0.5) +
  geom_sf(data = st_transform(world, wag_proj_string), fill="gray20", color=NA) +
  geom_sf(data = algaeworld, aes(color = percent_binned), size = 0.5) +#size  = glacier_km2 
  theme(panel.background = element_rect(fill = "grey50"),
        panel.grid.major = element_line(colour="gray30", size=0.1)) +
  labs(color = "% algal cover") +
  scale_color_manual(values = my_pal) +
  theme(legend.position = "none")
  # theme(legend.position = "bottom")

ggsave(here("figs/distribution/global_dist_map.png"), width = 10, height =10, units="in")

# use object to pattern, division 
# Inkscape export png crashes above 250 DPI--- use smaller bitmaps





# North america ----------------------------------------
# define an Albers equal area projection, with custom center and standard parallels 
aea_string <- function(center_on, parallels){
  lon_0 <- center_on[1]
  lat_0 <- center_on[2]
  lat_1 <- parallels[1]
  lat_2 <- parallels[2]
  paste0("+proj=aea +lat_0=",lat_0," +lon_0=",lon_0," +lat_1=",lat_1," +lat_2=",lat_2,
         " +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +type=crs")
}

america_proj <- aea_string(c(-135, 57),c(49, 63))

northamerica_algae <- northamerica_raw %>% 
  bin_percent() %>% 
  st_transform(america_proj) #%>% 
  # st_centroid(algaeworld)
max(northamerica_algae$fraction)

levels(northamerica_algae$percent_binned)

# theme just for the zoomed in maps
my_map_theme <-   theme(panel.background = element_rect(fill = "gray50", colour = NA), # make the ocean near black, 
                        panel.grid.major = element_line(colour="gray30", size=0.1)) # political boundaries in grey


ggplot() +  
  geom_sf(data = st_transform(states, america_proj), fill="gray20", color="gray40") + 
  geom_sf(data = northamerica_algae, aes(fill = percent_binned), color=NA, size=0.02) +
  coord_sf(xlim = c(-1.4e6, 1.4e6), ylim = c(-1e6, 0.9e6), expand = FALSE) + 
  scale_fill_manual(values = my_pal) +
  labs(fill = "% algae cover") +
  my_map_theme +
  theme(legend.position = "none")+
  ggspatial::annotation_scale(location = 'bl', width_hint = 0.15) +
  scale_x_continuous(breaks=seq(-150, -120, 10)) + # set longitude at 10 degree breaks
  scale_y_continuous(breaks=seq(50, 60, 5), position = "right")  # set latltude at 5 degree breaks
ggsave(here("figs/distribution/northamerica_dist_map.png"), width = 9, height = 6, units = "in")

# the same map, with out the red color
ggplot() +  
  geom_sf(data = st_transform(states, america_proj), fill="gray20", color="gray40") + 
  geom_sf(data = northamerica_algae, aes(fill = "white"), color=NA, size=0.02) +
  coord_sf(xlim = c(-1.4e6, 1.4e6), ylim = c(-1e6, 0.9e6), expand = FALSE) + 
  scale_fill_manual(values = my_pal) +
  labs(fill = "% algae cover") +
  my_map_theme +
  theme(legend.position = "none")+
  ggspatial::annotation_scale(location = 'bl', width_hint = 0.15) +
  scale_x_continuous(breaks=seq(-150, -120, 10)) + # set longitude at 10 degree breaks
  scale_y_continuous(breaks=seq(50, 60, 5), position = "right")  # set latltude at 5 degree breaks
ggsave(here("figs/distribution/northamerica_glaciers.png"), width = 9, height = 6, units = "in")



# Greenland --------------------------------------------------

gr_proj <- aea_string(c(-33, 66), c(49, 63))#c(60, 75)


greenland_algae <- greenland_raw %>% 
  bin_percent() %>% 
  st_transform(gr_proj) 

max(greenland_algae$fraction)
levels(greenland_algae$percent_binned)

ggplot() +  
  geom_sf(data = st_transform(polar_circles, gr_proj), color="gray30", linetype = "dashed", size = 0.5) +
  geom_sf(data = st_transform(world, gr_proj), fill="gray20", color="gray40") + 
  geom_sf(data = greenland_algae, aes(fill = percent_binned), color=NA, size=0.02) +
  coord_sf(xlim = c(-1.1e6, 1e6), ylim = c(-0.7e6, 0.7e6), expand = FALSE) + # for gris proj
  scale_fill_manual(values = my_pal) +
  labs(fill = "% algae cover") +
  my_map_theme +
  theme(legend.position = "bottom")+
  guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
  ggspatial::annotation_scale(location = 'bl', width_hint = 0.15) +
  scale_x_continuous(breaks=seq(-50, 0, 10)) + # set longitude at 10 degree breaks
  scale_y_continuous(breaks=seq(60,80, 5), position = "right") + # set latltude at 5 degree breaks
  theme(legend.text=element_text(size=14),legend.title=element_text(size=14))
ggsave(here("figs/distribution/greenland_dist_map.png"), width = 9, height = 6, units = "in")




# HMA ---------------------------------------
hma_proj <- aea_string(c(89, 38), c(25, 50) )

hma_algae <- hma_raw %>% 
  mutate(algae_km2 = as.numeric(algae_km2),
         glacier_km2 = as.numeric(glacier_km2),
         fraction = algae_km2/glacier_km2,
         percent = fraction*100) %>% 
  st_transform(hma_proj)

max(hma_algae$percent)

ggplot() +  
  geom_sf(data = st_transform(world, hma_proj), fill="gray20", color="gray40") + 
  geom_sf(data = hma_algae, aes(fill = algae_km2), color=NA, size=0.02) +
  coord_sf(xlim = c(-2e6, 2e6), ylim = c(-2e6, 2e6), expand = FALSE) + 
  scale_fill_gradient2(low="#ffffff", high = "#ff0000") +
  labs(fill = "Algal cover (km2)") +
  my_map_theme +
  ggspatial::annotation_scale(location = 'bl', width_hint = 0.15) +
  scale_x_continuous(breaks=seq(-50, 0, 10)) + # set longitude at 10 degree breaks
  scale_y_continuous(breaks=seq(60,80, 5), position = "right") 
ggsave(here("figs/distribution/greenland_dist_map.png"), width = 9, height = 6, units = "in")




# scratch-----------------------------------------------------
# 
# # create a zoomed in version
# ll <- 0.02
# ggplot() + 
#   geom_sf(data = polar_circles, color="gray30", linetype = "dashed", size =ll) +
#   geom_sf(data = st_transform(world, wag_proj_string), fill="gray20", color="gray50", size = ll) +
#   geom_sf(data = st_transform(states, wag_proj_string), fill=NA, color="gray50", size = ll) +
#   geom_sf(data = algaeworld, shape =15, aes(color = percent_binned, size  = glacier_km2)) +
#   theme(panel.background = element_rect(fill = "gray80"),
#         panel.grid.major = element_line(colour="gray50", size=ll)) +
#   labs(color = "% algal cover", size = "Glacier area (km^2)") +
#   scale_color_manual(values = my_pal) +
#   scale_size(range = c(0, 1)) 
# # ggspatial::annotation_scale(location = 'bl', width_hint = 0.15) # not working
# 
# ggsave(here("figs/distribution/global_dist_map10.png"), width = 30, height =30, units="in", limitsize = FALSE)
# 



    
    
# Per glacier % --------------------------

per_glac_2 <- per_glac %>% 
  mutate(range = str_split_fixed(basename(fn), pattern="_",n=2)[,2] %>% str_remove(".csv")) %>% # distinct(range)
  select(range, glac_id, algae_frac) %>% print()

# stats per region
ord <- c("<0.1%", "<1%","1-5%", "5-10%", "10-25%", "25-50%", "50-75%", "75-100%")
per_glac_3 <- per_glac_2 %>%
  mutate(algae_frac_bin = case_when(
    # algae_frac==0~"0",
    algae_frac<0.001~"<0.1%",
    algae_frac<0.01~"<1%",
    algae_frac<0.05~"1-5%",
    algae_frac<0.1~"5-10%",
    algae_frac<0.25~"10-25%",
    algae_frac<0.5~"25-50%",
    algae_frac<0.75~"50-75%",
    algae_frac<1~"75-100%"
  ) %>% fct_relevel(ord))

per_glac_3 %>% count(range, algae_frac_bin)

per_glac_2 %>% 
  ggplot(aes(algae_frac)) +
  geom_histogram()+
  facet_wrap(vars(range))

per_glac_2 %>%
  filter(algae_frac>0.05) %>% 
  ggplot(aes(algae_frac)) +
  geom_histogram()+
  facet_wrap(vars(range))

per_glac_2 %>% 
  filter(algae_frac>0.1) %>%
  # the number of glaciers in North America with >1km2 area and >10% algal cover
  count(range) %>% filter(!range %>% str_detect("Greenland|kamchatka|orway|svalbard")) %>% summarise(sum(n))
1305/8419 # percent of glaciers in NA with >10% algal cover

# count glaciers per range
total_glac <-  per_glac_2 %>% 
  count(range) %>% 
  rename(total_glac=n)

my_order <- c("alaska", "chugatch", "wrangellStElias",
              "northCoast","southCoast", 
              "mackenzie","northernRockies","canadianRockies",
              "westCentralGreenland", "southwestGreenland", "southeastGreenland", "eastCentralGreenland",
              "norway","northNorway","svalbard", "kamchatka")


per_glac_3 %>%
  count(range, algae_frac_bin) %>% 
  pivot_wider(names_from = "algae_frac_bin", values_from = "n", values_fill =0) %>% 
  left_join(total_glac) %>% 
  mutate(range = range %>% fct_relevel(my_order)) %>% 
  arrange(range) %>% 
  # relocate(range, total_glac) %>% 
  kbl() %>% 
  kable_classic(full_width = F, html_font = "Cambria")
ggsave(here("figs/distribution/per_glacier_table.png"))


# select glaciers for albedo analysis
per_glac_2 %>% 
  filter(!range %>% str_detect("Greenland|kamchatka|orway|svalbard")) %>% 
  filter(algae_frac>0.1) %>% 
  rownames_to_column(var = "gt10_rowname") %>% 
  write_csv(here("data/s2_classifier_map/algae_per_glacier/glac_ids_gt_10.csv"))
  
per_glac_2 %>% 
  filter(!range %>% str_detect("Greenland|kamchatka|orway|svalbard")) %>% 
  filter(algae_frac>0.01) %>% # nrow() 3486
  rownames_to_column(var = "gt1_rowname") %>% 
  write_csv(here("data/s2_classifier_map/algae_per_glacier/glac_ids_gt_1.csv"))









# SCRATCH------------------ ---------------
z_local_dir <- here("figs/distribution/v5/")
z_goog_paths <- drive_ls(path = "figs", pattern = "kml") %>% print()
z_local_paths <- paste0(z_local_dir,z_goog_paths %>% pull(name))

# download from google drive
# z_goog_paths %>% 
#   pull(id) %>% 
#   walk2(z_local_paths, ~drive_download(as_id(.x), path = .y, overwrite = TRUE)) #

norway <- st_read(z_local_paths[z_local_paths %>% str_detect("norge")] )
kamchatka <- st_read(z_local_paths[z_local_paths %>% str_detect("kam")] )
# # also download the ridge kml
# ridge_goog_path <- drive_ls(path = "figs", pattern = "ridge_features.kml")
# ridge_local_path <- paste0(z_local_dir,ridge_goog_path %>% pull(name))
# drive_download(ridge_goog_path, ridge_local_path)


# Norway -------------------------------------------

norway_proj <- aea_string(c(16,64), c(60, 68))

norway_tidy <- tidy_algae_map(norway, norway_proj)


ggplot() +  
  geom_sf(data = st_transform(world, norway_proj), fill="gray20", color="gray40") +
  geom_sf(data = norway_tidy, aes(fill = fraction_cut), color="black", size=0.05) + #color = NA)+#
  coord_sf(xlim = c(-0.75e6, 0.5e6), ylim = c(-1e6, 0.75e6), expand = FALSE) +
  labs(fill = "% algae cover") +
  scale_fill_manual(values = my_discrete_pal, limits = my_fraction_cut_vals) +
  my_map_theme +
  scale_x_continuous(breaks=seq(10, 20, 10)) + # set longitude at 10 degree breaks
  scale_y_continuous(breaks=seq(60, 65, 5)) # set latltude at 5 degree breaks

ggsave(here("figs/distribution/v5/norway_map.png"))



# Kamchatka -------------------------------------

kamchatka_proj <- aea_string(c(159,56), c(50,60))
kamchatka_tidy <- tidy_algae_map(kamchatka, kamchatka_proj)

ggplot() +  
  geom_sf(data = st_transform(world, kamchatka_proj), fill="gray20", color="gray40") +
  geom_sf(data = kamchatka_tidy, aes(fill = fraction_cut), color="black", size=0.05) + #color = NA)+#
  coord_sf(xlim = c(-0.75e6, 0.5e6), ylim = c(-1e6, 0.75e6), expand = FALSE) +
  labs(fill = "% algae cover") +
  scale_fill_manual(values = my_discrete_pal, limits = my_fraction_cut_vals) +
  my_map_theme +
  scale_x_continuous(breaks=seq(150, 160, 10)) + # set longitude at 10 degree breaks
  scale_y_continuous(breaks=seq(50, 60, 5)) # set latltude at 5 degree breaks

ggsave(here("figs/distribution/v5/kamchatka_map.png"))



# cover %>%
#   mutate(code_range = code_range %>% fct_relevel(my_order_coderange)) %>% 
#   ggplot(aes(y = code_range, x = algae_area_km2)) +
#   geom_col(fill = "gray60") +
#   theme_minimal() +
#   geom_errorbarh(aes(xmin = algae_area_km2 - (fpr*algae_area_km2), xmax = algae_area_km2, height = 0.1))
# 
# # biogeography is the relevant factor, since total area is just a function of the glacier area
# 
# cover %>% 
#   group_by(region) %>% 
#   summarise(algae_area_km2 = sum(algae_area_km2),
#             glacier_area_km2 = sum(glacier_area_km2),
#             percent = algae_area_km2/glacier_area_km2)





# # error data
# drive_error <- drive_ls(path = "s2_classifier_map", pattern = "s2sr_error")
# local_error_path <- paste0(local_dir, "s2sr_error.csv")
# # drive_download(as_id(drive_error), path = local_error_path, overwrite = TRUE) # update
# error <- read_csv(local_error_path)
# glimpse(error)

# Tidy the data -----------------------------------------------

# # note 2022-10-07, take mean of FPR grouped by range, version
# err_clean <- error %>% 
#   janitor::clean_names() %>% 
#   select(region, range, version, fpr_pts, fpr) %>% 
#   mutate(fpr = if_else(fpr %in% c("#DIV/0!", "#VALUE!"), NA_character_, fpr) %>% as.numeric())

# cover <- cover_raw %>% 
#   janitor::clean_names() %>% 
#   select(range, region, version, algae_area_m2, glacier_area_m2) %>% 
#   mutate(algae_area_km2 = algae_area_m2 / (1000*1000), 
#          glacier_area_km2 = glacier_area_m2 / (1000*1000),
#          glacier_minus_algae_km2 = glacier_area_km2 - algae_area_km2,
#          percent = algae_area_km2 / glacier_area_km2, # should be relabelled as fraction
#          .keep = "unused") %>% 
#   filter(glacier_area_km2>0) %>% 
#   inner_join(err_clean) %>% 
#   mutate(fpr_pct = percent*fpr,
#          fpr = if_else(is.na(fpr),0,fpr)) %>% # since we are filtering by min(fpr), need to convert NA to 0
#   group_by(range) %>% 
#   filter(fpr == min(fpr)) %>%  # choose version that minimizes error
#   filter(version == max(version)) %>% # if FPR is identical, choose the more recent version
#   ungroup() %>% 
#   filter(!(algae_area_km2==0 & fpr_pts>0))# remove entries lacking area but with TP algae found in error assessment
# glimpse(cover)
# 
# # error rates missing for:
# filter(cover, is.na(fpr))

  


# range_codes <- c(
#   paste0("N", 1:10),
#   paste0("S", 1:3),
#   paste0("E", 1:3),
#   paste0("Gr", 1:7),
#   paste0("Ar", 1:8),
#   paste0("C", 1:8),
#   paste0("Z", 1),
#   paste0("K", 1),
#   paste0("An", 1:4)
# ) 

# names(my_order) <- range_codes
# 
# cover <- cover %>% 
#   left_join(enframe(my_order, name = "code", value = "range")) %>% 
#   mutate(code_range = paste(code, range, sep = "_"))


# my_order <- my_order[my_order %in% cover$range] # avoid error with mismatched lengths
# my_order <- rev(my_order) # horizontal axis barplot will put first at bottom, want on top
# 
# my_order_coderange <- paste(names(my_order), my_order, sep="_")
# 
# 
# 
# cover %>% 
#   filter(range %in% c("norway", "kamchatka"))
# 
# 
# # save data for use with other scripts
# # write_csv(cover, here("data/s2_classifier_map/cover_final.csv"))
# 
# 
# 
# 
# 
# 
# 









# SCRATCH --------------------------

# convert fraction to categorical (easier to view)
num_to_fct <- function(sf, breaks, labs){
  labs <- enquo(labs)
  
  sf$pct_cat <- as.numeric(sf$fraction) %>% 
    cut(breaks = breaks, include.lowest = TRUE) 
  
  lvls <- sf$pct_cat %>% levels() 
  names(lvls) <- my_labs
  # return sf, with recoded levels
  sf$pct_cat <- sf$pct_cat %>% fct_recode(lvls[1],lvls[2],lvls[3],lvls[4],lvls[5])
  sf
}
# 
# 
# glacier_outlines <- st_read(here("figs/archive/glacierVectors.kml"))
# #plot(s[1])
# 
# # annotate the vectors with percent area
# algae_map_data <- glacier_outlines %>% 
#   mutate(lng = as.numeric(lng), 
#          lat = as.numeric(lat)) %>% 
#   inner_join(cover) 
# 
# 
#   # rev(kovesi.linear_kryw_5_100_c64(7))
# 
# # percent area covered map
# ggplot() + 
#   geom_map(data = background_dat, map=background_dat,
#            aes(x=long, y=lat, map_id=region),
#            fill="gray20", size=0) +
#   # geom_sf(aes(fill = percent), color = NA, data = algae_map_data) +#, color = percent
#   # scale_color_gradientn(colors = pal)  +
#   scale_fill_gradientn(colors = pal) +
#   # geom_text_repel(aes(x = lng, y = lat, label = code), data = algae_map_data) +
#   geom_tile(data = test_df, aes(x = x, y=y, fill = value)) +
#   labs(fill = "fraction")
# 
# 
# # how to get the 
# 
# 
# 
# 
# # in polar orthographic projection
# library(sf)
# library(ggplot2)
# library(mapview)
# library(lwgeom)
# library(rnaturalearth)
# 
# # world data
# world <- rnaturalearth::ne_countries(scale = 'small', returnclass = 'sf')
# 
# # Fix polygons so they don't get cut in ortho projection
# world  <- st_cast(world, 'MULTILINESTRING') %>%
#   st_cast('LINESTRING', do_split=TRUE) %>%
#   mutate(npts = npts(geometry, by_feature = TRUE)) %>%
#   st_cast('POLYGON')
# library(raster)
# library(stars) # can handle curved reference coords
# ss <- read_stars("~/Downloads/randomRaster.tif")
# rr <- raster("~/Downloads/randomRasterMercator.tif")
# rr
# 
# test_spdf <- as(rr, "SpatialPixelsDataFrame")
# test_df <- as.data.frame(test_spdf)
# colnames(test_df) <- c("value", "x", "y")
# as_tibble(test_df) #%>% pivot_wider(names_from = "x", values_from = "value")
# 
# 
# # map
# ggplot() +
#   geom_sf(data=world, fill="gray20", color = "gray20", size = 0.1) + # remove color and size args to show national borders
#   theme(panel.background = element_rect(fill = 'gray50'),
#         panel.grid.major = element_line(size=0.265)) +
#   geom_sf(aes(fill = percent), color = NA, data = algae_map_data) +#, color = percent
#   # scale_color_gradientn(colors = pal)  +
#   scale_fill_gradientn(colors = pal)+
#   coord_sf( crs= "+proj=ortho +lat_0=90 +lon_0=0")  
# 
# # how to make the oceans one color and the background white?
# # shift center slightly downwards top is cut off
# 
# # south pole
# ggplot() +
#   geom_sf(data=world, color="gray80", fill="gray20") +
#   coord_sf( crs= "+proj=ortho +lat_0=-90 +lon_0=0") +
#   theme(panel.background = element_rect(fill = 'gray50'))

# clip in a circle :
# https://gis.stackexchange.com/questions/79350/plotting-a-certain-range-in-azimuthal-projection-using-ggplot2

# more tips for polar stereographic plotting in ggplot
# https://stackoverflow.com/questions/48816773/polar-stereographic-map-in-r


# scratch ---------
# # this is a poor visualization, it draws too much attention to the glacier area and not enough to the algae
# long_dat <- cover %>% 
#   mutate(percent = paste(as.character(round(percent,3)*100), "%"),
#          range = range %>% fct_relevel(my_order)) %>% # order factors for plotting
#   rename(algae = algae_area_km2, glacier = glacier_minus_algae_km2) %>% 
#   pivot_longer(cols = c(algae,glacier), names_to = "type", values_to = "area_km2") 
# 
# long_dat %>%
#   ggplot(aes(y = range, x = area_km2, fill = type)) +
#   geom_col() +
#   theme_minimal() + 
#   facet_wrap(vars(region), scales = "free_y") +
#   # annotate with percent
#   geom_text(aes(label = percent), data = long_dat %>% filter(type=="glacier"), nudge_x = 3500)
# 
# 
