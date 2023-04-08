#* map of modis trends in north america
#* 
#* 
#* 
#* UPSTREAM
#* 09_trend.js
#* 
#* 
#* OUTPUT
#* north_america_median.png (fig 3)
#* mk_p_val_intensity.png (the p values, using the mean per grid cell, then MK test in R method)
#* mean_max_rgnd_composite_vs_year.png (the per range RGND trends)
#* Mann-Kendall trend analysis, for both grid cells, and per ranges

library(tidyverse)
library(googledrive)
library(here)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggpubr)
library(scales)
library(Kendall)
library(broom)
library(cowplot) # ggdraw, sim ggpubr::ggarrange
library(zyp) 

# read --------------------------------

# read into R
trends <- st_read(here("figs/modis/rgnd_intensity/westcoast_grid/westcoast_sens_by_elv_20k_bc_albers.kml"))

range_meta_raw <- read_csv(here("data/cells_annotated.csv"))

# the annual data, mean values per grid cell
meanmax <- dir_ls(here("data/time_series/modis/"), 
                  regexp="modis_mean_max_20k_bcaea.csv") %>%
  map_df(read_csv, id = "filename")

# regional mean/median RGND
algae_raw <- read_csv(here("data/time_series/modis/regional_yrly_rgnd_stats.csv"))
algae_raw_exp <- read_csv(here("data/time_series/modis/regional_yrly_rgnd_stats_expanded.csv")) # with the 5k downslope buffer
# same data, but just for BC coast and interior
algae_raw_2r <- read_csv(here("data/time_series/modis/regional_yrly_rgnd_stats_coastvinterior.csv"))
algae_raw_exp_2r <- read_csv(here("data/time_series/modis/regional_yrly_rgnd_stats_expanded_coastvinterior.csv"))


algae_raw_nr <- read_csv(here("data/time_series/modis/regional_yrly_rgnd_stats_newranges.csv"))

# tidy ------------------------------------

# define an Albers equal area projection, with custom center and standard parallels 
aea_string <- function(center_on, parallels){
  lon_0 <- center_on[1]
  lat_0 <- center_on[2]
  lat_1 <- parallels[1]
  lat_2 <- parallels[2]
  paste0("+proj=aea +lat_0=",lat_0," +lon_0=",lon_0," +lat_1=",lat_1," +lat_2=",lat_2,
         " +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +type=crs")
}

bc_aea <- aea_string(c(-135, 57), c(49, 63))

trends2 <-
  trends %>% 
  as_tibble() %>% 
  select(system_index, hi_sens_count:geometry) %>% 
  rename(cell_id = system_index) %>% 
  mutate_at(vars(matches("sens")), as.numeric) %>% 
  mutate_at(vars(matches("mean")), ~.x/1000) %>% # I multiplied by 1000 to simplify computation, so return to units of RGND 
  mutate(diff = hi_sens_mean - lo_sens_mean) %>% 
  st_as_sf() %>% 
  st_transform(bc_aea)

range_meta <- range_meta_raw %>% 
  janitor::clean_names() %>% 
  rename(cell_id =system_index) %>% 
  mutate(cell_id = str_sub(cell_id, 3,8))

meanmax2 <- meanmax %>%
  clean_names() %>% 
  select(-system_index, -geo) %>% 
  mutate(elev = basename(filename) %>% str_sub(1,2), 
         mean = mean/1000,std_dev=std_dev/1000,
         .keep="unused")

my_clean <- function(df){
  df %>% 
    janitor::clean_names() %>% 
    select(range, year, mean) %>% 
    mutate(rgnd = mean/1000, .keep="unused")
}

rr <- my_clean(algae_raw)
rr_exp <- my_clean(algae_raw_exp)
rr_2r <- my_clean(algae_raw_2r)
rr_exp_2r <- my_clean(algae_raw_exp_2r)

# plot ----------------

trend_pal <- c('#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac') %>% rev()

# more white in center
trend_pal2 <- c('#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7','#f7f7f7','#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac') %>% rev()
viridis_trend_pal <- c( '#21918c','#ffffff','#fde725')

# load background map data, transform to desired proj
# states <- ne_download(scale = 50, type = 'states', category = 'cultural', returnclass = "sf")
states <- rnaturalearth::ne_states(c("united states of america", "canada"), returnclass = "sf")

my_map_theme <-   theme(panel.background = element_rect(fill = 'gray10', colour = NA),
                        panel.grid.major = element_line(colour="gray30", size=0.1))



## median Sens per grid cell map ----------
ggplot() + 
  geom_sf(data = st_transform(states, bc_aea), fill="gray20") +
  geom_sf(data = trends2, aes(fill = all_sens_median), color="black", size=0.1 )+ #color = NA)+#
  coord_sf(xlim = c(-1.1e6, 1.4e6), ylim = c(-0.8e6, 0.9e6), expand = FALSE) +
  scale_fill_gradient2(high=muted("red"), low = muted("blue")) +
  # theme_void() +
  theme(panel.grid.major = element_line(colour="gray70"),
        panel.background = element_rect(fill = "gray40"),
        plot.background = element_rect(fill = "white", color=NA)) +
  labs(fill = "Median Sen's slope")
ggsave(here("figs/modis/rgnd_intensity/north_america_median.png"))

trends2 %>% 
  as_tibble() %>% 
  arrange(-all_sens_median) %>% 
  filter(all_sens_median<0.4)
# max: 58,-25
# max below 0.4 -19,20


# annotate by range
trends2 %>% 
  left_join(range_meta) %>% 
  as_tibble() %>% 
  drop_na(diff) %>% 
  group_by(range) %>% 
  summarise(mean(diff), median(diff))
# not much difference by elevation.. interesting that interior south has expanded downhill??


## mean Sens per grid cell map ----------
# caution: high discrepancy between mean and median!
my_north_map <- function(fill, title, count_thresh=5, count_var, 
                         hi=0.001, hi_break=0.001, pal = trend_pal){
  
  # adjust the highest two values in canadian rockies---indicate on legend >0.001
  trends2 <- trends2 %>% 
    mutate_at(vars(matches("mean")), ~if_else(.x>0.001,0.001,.x)) %>% 
    mutate(diff = if_else(diff>0.001,0.001,diff))
  
  if(count_thresh>5){
    trends_filt <- trends2 %>% 
      filter({{count_var}}>=count_thresh) #inclusive
  }else{
    trends_filt <- trends2 
  }
  
  ggplot() + 
    geom_sf(data = st_transform(states, bc_aea), fill="gray20") +
    geom_sf(data = trends_filt, aes(fill = {{fill}}), color="black", size=0.1) + #color = NA)+#
    coord_sf(xlim = c(-1.1e6, 1.4e6), ylim = c(-0.8e6, 0.9e6), expand = FALSE) +
    scale_fill_gradientn(colours = viridis_trend_pal, limits=c(-hi,hi), breaks = c(-hi_break, 0, hi_break)) +
    # theme_void() +
    theme(panel.grid.major = element_line(colour="gray70"),
          panel.background = element_rect(fill = "gray40"),
          plot.background = element_rect(fill = "white", color=NA)) +
    labs(fill = title)
}

my_north_map(all_sens_mean, "23-yr trend", 6, all_sens_count) 
ggsave(here("figs/modis/rgnd_intensity/fig3a.png"), width=7,height=4,units="in")

## difference between the means (evidence for upslope expansion) ----------
my_north_map(diff, "Upslope\nexpansion", 6, all_sens_count)  #hi=0.0008, hi_break =0.0005, pal =trend_pal2
ggsave(here("figs/modis/rgnd_intensity/fig3.png"), width=7,height=9,units="in")

##  hi mean Sens ------------
my_north_map(hi_sens_mean, "High elev\n 23-yr trend", 6, hi_sens_count)
ggsave(here("figs/modis/rgnd_intensity/hi_sens_na.png"), width=7, height=4)
##  lo mean Sens ------------
my_north_map(lo_sens_mean, "Low elev\n 23-yr trend", 6, lo_sens_count)
ggsave(here("figs/modis/rgnd_intensity/lo_sens_na.png"), width=7, height=4)

# for reference, 0.001 RGND corresponds to a percent area coverage increase of 0.23%, a quarter of a percent, *22 years a 5% total increase
# see algae_abundance_regressions::frac_mod_simple (I used a rgnd^2 model in my final paper)


## regional RGND trends --------------------- 
sen <- function(..., weights = NULL) {
  mblm::mblm(...)
}

# look at the regional RGND trends... (same code as in modis_trend_map_2.R)
rr %>%
  mutate(range =range %>% fct_relevel("alaska","chugatch","wrangellStElias", "northCoast", "southCoast", "mackenzie", "northernRockies","canadianRockies")) %>% 
  ggplot(aes(x = year, y = rgnd)) +
  geom_point()+
  geom_smooth(method = sen) +
  facet_wrap(vars(range)) +
  labs(y="mean(max rgnd composite)")
ggsave(here("figs/modis/rgnd_intensity/mean_max_rgnd_composite_vs_year.png"))


# same, with extended downhill slope
rr_exp %>%
  mutate(range =range %>% fct_relevel("alaska","chugatch","wrangellStElias", "northCoast", "southCoast", "mackenzie", "northernRockies","canadianRockies")) %>% 
  ggplot(aes(x = year, y = rgnd)) +
  geom_point()+
  geom_smooth(method = sen) +
  facet_wrap(vars(range)) +
  labs(y="mean(max rgnd composite)", title = "23 year trend with downhill extended patches")
ggsave(here("figs/modis/rgnd_intensity/mean_max_rgnd_composite_vs_year_extended.png"))

# same but for the new ranges
algae_raw_nr %>% 
  mutate(mean = mean/1000) %>% 
  filter(range %in% c("coast_south", "interior_south")) %>% 
  mutate(range =range %>% fct_relevel("alaska","coast_north","coast_south", "interior_north", "interior_south")) %>%
  mutate(range = fct_recode(range, "BC Coast"= "coast_south", "Southern Interior Ranges"="interior_south")) %>% 
  ggplot(aes(x = year, y = mean)) +
  geom_point() +
  geom_smooth(method = "lm" ) +
  facet_wrap(vars(range)) +
  labs(y="Algal bloom intensity (mean RGND)", x = "Year") +
  theme_minimal()
ggsave(here("figs/modis/rgnd_intensity/mean_max_rgnd_composite_vs_year_newranges.png"),
       width = 6, height = 3, units="in")

  


# model --------------------------------------

## MK mean RGND per grid cell---------------
mk_res <- meanmax2 %>% 
  filter(elev=="bo") %>% 
  group_by(cell_id, elev) %>% 
  group_modify(~pull(.x, mean) %>% MannKendall() %>% broom::tidy()) %>% 
  ungroup() %>%
  clean_names()
  
# join test results to above df
trends3 <- trends2 %>% 
  left_join(mk_res)

my_cell_id <- "-19,20"#"58,-25" #
filter_dat <- meanmax2 %>% 
  filter(elev=="bo",cell_id==my_cell_id)
ggplot(filter_dat, aes(x = year, y =mean))+
  geom_point() +
  geom_smooth(method = "lm") +
  labs(title = my_cell_id) 
lm(mean~year, data = filter_dat) %>% predict(data.frame(year = c(2000,2023)))
# a near doubling

## plot P values -------------------
p_pal <- c('#fde725','#fde725', '#5ec962','#21918c','#3b528b',rep('#440154',14))
a1_raw <- ggplot() + 
  geom_sf(data = st_transform(states, bc_aea), fill="gray20") +
  geom_sf(data = trends3, aes(fill = p_value), color="black", size=0.1) + #color = NA)+#
  coord_sf(xlim = c(-1.1e6, 1.4e6), ylim = c(-0.8e6, 0.9e6), expand = FALSE) +
  scale_fill_gradientn(colours = p_pal, breaks = c(0, 0.05)) +
  theme(panel.grid.major = element_line(colour="gray70"),
        panel.background = element_rect(fill = "gray40"),
        plot.background = element_rect(fill = "white", color=NA)) +
  labs(fill = "Mann-Kendall\np value")
a1_raw
ggsave(here("figs/modis/rgnd_intensity/mk_p_val_intensity.png"))
# extract the legend from A1
p_val_legend <- cowplot::get_legend(a1_raw) 

a1 <- a1_raw +
  theme(legend.position = "none")


## MK diff hi vs lo ---------------------
elv_dif_mk <- meanmax2 %>% 
  filter(elev !="bo") %>% 
  select(cell_id, year, elev,mean) %>% 
  distinct() %>% 
  pivot_wider(names_from = "elev",values_from = "mean") %>% 
  mutate(dif = hi-lo) %>% 
  group_by(cell_id) %>% 
  group_modify(~pull(.x, dif) %>% MannKendall() %>% broom::tidy()) %>% 
  ungroup() %>%
  clean_names()

trends4 <- trends2 %>% 
  left_join(elv_dif_mk)
## plot P values ------------
b1 <- ggplot() + 
  geom_sf(data = st_transform(states, bc_aea), fill="gray20") +
  geom_sf(data = trends4, aes(fill = p_value), color="black", size=0.1) + #color = NA)+#
  coord_sf(xlim = c(-1.1e6, 1.4e6), ylim = c(-0.8e6, 0.9e6), expand = FALSE) +
  scale_fill_gradientn(colours = p_pal, breaks = c(0, 0.05)) +
  theme(panel.grid.major = element_line(colour="gray70"),
        panel.background = element_rect(fill = "gray40"),
        plot.background = element_rect(fill = "white", color=NA)) +
  labs(fill = "MK\np value", tag = "B1") +
  theme(legend.position = "none")
b1



## Sen's slope for regional trends -----------------



algae_raw_nr %>% 
  group_by(range) %>% 
  nest() %>% 
  mutate(vec = map(data, ~pull(.x, mean)),
         mod = map(vec, MannKendall),
         mod_res = map(mod, tidy)) %>% 
  unnest(mod_res)
  


# get the intercept and Sen's slope
rr_2r %>% 
  mutate(year = year-2000) %>% # set x intercept (x=0) to year 2000
  group_by(range) %>% 
  group_modify(~zyp.sen(rgnd~year, data = .x) %>% coef() %>% enframe()) %>% 
  ungroup() %>% 
  pivot_wider() %>%
  rename(slope=year) %>% 
  mutate(change = slope*23,
         pct = (change/Intercept)*100)
# rockies, 50% increase relative to 2000 levels
#  mackenzie, 108% increase, southcoast 70%, northcoast 112%


rr_2r %>% 
  group_by(range) %>% 
  group_modify(~pull(.x, rgnd) %>% sens.slope() %>% tidy()) %>% 
  arrange(p.value)
# p<0.05: mackenzie, canadian rockies, south coast, north coast 
# p<0.1: northern rockies 

rr_exp %>% 
  group_by(range) %>% 
  group_modify(~pull(.x, rgnd) %>% sens.slope() %>% tidy()) %>% 
  arrange(p.value)
# still statistically significant for mackenzie and canadian rockies. 
# p<0.1: north coast, northern rockies, south coast. 






algae_raw_nr %>% 
  select()















# SCRATCH --------------------------------------------------

# plot raw trend vals with green square around significant cells
# edit the trend map above to highlight significant cells
a2_raw <- a +
  geom_sf(data = trends3 %>% filter(p_value<0.05),color = "purple", size=0.3, fill=NA) +
  coord_sf(xlim = c(-1.1e6, 1.4e6), ylim = c(-0.8e6, 0.9e6), expand = FALSE) +
  labs(tag="A2") +
  theme(legend.position="bottom")

slope_legend <- cowplot::get_legend(a2_raw) 

a2 <- a2_raw +
  theme(legend.position = "none")


b2 <- b +
  geom_sf(data = trends4 %>% filter(p_value<0.05),color = "purple", size=0.3, fill=NA) +
  coord_sf(xlim = c(-1.1e6, 1.4e6), ylim = c(-0.8e6, 0.9e6), expand = FALSE) +
  labs(tag = "B2")+
  theme(legend.position = "none")
b2


# assemble the final supp fig with modis p values
cowplot::ggdraw(plot_grid(plot_grid(a1,a2,b1,b2,ncol=2),
                          plot_grid(p_val_legend,slope_legend, ncol=2),
                          nrow=2,
                          rel_heights =c(1,0.1) ))
ggsave(here("figs/modis/rgnd_intensity/fig3_pval_supps.png"), height = 6, width =10, units="in")

# 
# # zoomed in maps-----------------------------------------------------------------------------------------
# 
# # data = sum rgnd / year per 10 km albers grid cell
# # load into R to compute the time series statistics
# 
# # download the data from google drive
# zz_local_dir <- here("figs/modis/rgnd_intensity/10km_grid/")
# zz_goog_paths <- drive_ls(path = "modis", pattern = "10k_grid_sum_rgnd_no_qa.kml") %>% print()
# zz_local_paths <- paste0(zz_local_dir,zz_goog_paths %>% pull(name))
# # download the files
# zz_goog_paths %>%
#   pull(id) %>%
#   walk2(zz_local_paths, ~drive_download(as_id(.x), path = .y))#, overwrite = TRUE
# 
# # read into R
# read_plus <- function(path){
#   region <- path %>% 
#     basename() %>% 
#     str_sub(1,5)
#   
#   st_read(path) %>% 
#     add_column(region = region) %>% 
#     mutate(grid_index = str_split_fixed(system_index,"_",2)[,2]) %>% 
#     mutate(across(.cols = year:sum, as.numeric))
# } 
# raw_10k <- zz_local_paths %>% map_df(read_plus) # the per year RGND data
# 
# 
# library(Kendall)
# library(broom)
# 
# mk_res <- raw_10k %>%  
#   arrange(year) %>% 
#   group_by(grid_index, region) %>%
#   nest() %>% 
#   ungroup() %>% 
#   mutate(mk = map(data, ~pull(.x, sum) %>% MannKendall() %>% broom::tidy()) ) %>% 
#   unnest(mk) %>% 
#   janitor::clean_names() %>% print()
# 
# mk_res %>% 
#   filter(p_value<0.1) %>%
#   ggplot(aes(statistic)) +
#   geom_histogram() +
#   geom_vline(xintercept = 0)
# 
# 
# 
# # BELOW HERE edit for brevity ------------------------------
# 
# 
# 
# 
# 
# 
# ggplot() +  
#   geom_sf(data = states, fill="gray20", color="gray40") +
#   geom_sf(data = mk_sf_rprj, aes(fill = statistic), color="black", size=0.05) + #color = NA)+#
#   coord_sf(xlim = c(-1.3e6, 1.4e6), ylim = c(-1e6, 1e6), expand = FALSE) +
#   labs(fill = "MK estimate") +
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
#   my_map_theme +
#   theme(legend.position="bottom") 
# 
# ggsave(here("figs/modis/rgnd_intensity/10km_grid/map_all_p.png"))
# 
# 
# ggplot() +  
#   geom_sf(data = states, fill="gray20", color="gray40") +
#   geom_sf(data = mk_sf_rprj %>% 
#             filter(p_value<0.05), aes(fill = statistic), color="black", size=0.05) + #color = NA)+#
#   coord_sf(xlim = c(-1.3e6, 1.4e6), ylim = c(-1e6, 1e6), expand = FALSE) +
#   labs(fill = "MK estimate") +
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
#   my_map_theme +
#   theme(legend.position="bottom") 
# 
# ggsave(here("figs/modis/rgnd_intensity/10km_grid/map_signif.png"))
# 
# 
# 
# # same process for GREENLAND ------------------------------------
# 
# # read into R
# ggrids <- st_read(z_local_paths[z_local_paths %>% str_detect("gris") ])
# 
# ggrid_data <- ggrids %>% st_drop_geometry() %>% 
#   as.data.frame() %>% 
#   mutate(grid_index = str_split_fixed(system_index,"_",2)[,2]) %>% 
#   select(grid_index, year, count, sum) %>% 
#   as_tibble() %>% 
#   mutate(across(.cols = year:sum, as.numeric))
# 
# # lookup table
# ggrid_lookup <- ggrids %>% 
#   as.data.frame() %>% 
#   mutate(grid_index = str_split_fixed(system_index,"_",2)[,2]) %>% 
#   distinct(grid_index, geometry)
# 
# 
# gmk_res <- ggrid_data %>%  
#   arrange(year) %>% 
#   group_by(grid_index) %>%
#   nest() %>% 
#   ungroup() %>% 
#   mutate(mk = map(data, ~pull(.x, sum) %>% MannKendall() %>% broom::tidy()) ) %>% 
#   unnest(mk) %>% 
#   janitor::clean_names() %>% print()
# 
# gmk_res %>% 
#   filter(p_value<0.1) %>%
#   ggplot(aes(statistic)) +
#   geom_histogram() +
#   geom_vline(xintercept = 0)
# 
# # transform into an sf object
# gmk_sf <- gmk_res %>% 
#   left_join(ggrid_lookup) %>% 
#   select(grid_index, statistic, p_value, geometry) %>% 
#   st_as_sf()
# 
# gris_proj <- aea_string(c(-45, 65), c(60, 70))
# 
# # reproject
# gmk_sf_rprj <- gmk_sf %>% st_transform(gris_proj)
# 
# 
# world <- ne_countries(scale = "medium", returnclass = "sf")
# world_aea <- sf::st_transform(world, gris_proj)
# 
# 
# 
# 
# ggplot() +  
#   geom_sf(data = world_aea, fill="gray20", color="gray40") +
#   geom_sf(data = gmk_sf_rprj, aes(fill = statistic), color="black", size=0.05) + #color = NA)+#
#   coord_sf(xlim = c(-0.6e6, 0.9e6), ylim = c(-0.7e6, 1e6), expand = FALSE) +
#   labs(fill = "MK estimate") +
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
#   my_map_theme +
#   theme(legend.position="bottom") 
# 
# ggsave(here("figs/modis/rgnd_intensity/10km_grid/gris_map_all_p.png"))
# 
# 
# ggplot() +  
#   geom_sf(data = world_aea, fill="gray20", color="gray40") +
#   geom_sf(data = gmk_sf_rprj %>% 
#             filter(p_value<0.05), aes(fill = statistic), color="black", size=0.05) + #color = NA)+#
#   coord_sf(xlim = c(-0.6e6, 0.9e6), ylim = c(-0.7e6, 1e6), expand = FALSE) +
#   labs(fill = "MK estimate") +
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
#   my_map_theme +
#   theme(legend.position="bottom") 
# 
# ggsave(here("figs/modis/rgnd_intensity/10km_grid/gris_map_signif.png"))
# 
# 
# 
# 
# ## same for NORWAY --------------------------------
# 
# # read into R
# ngrids <- st_read(z_local_paths[z_local_paths %>% str_detect("norge5") ])
# 
# nmk_res <- mk_kml(ngrids)
# 
# nmk_res %>% 
#   filter(p_value<0.1) %>%
#   ggplot(aes(statistic)) +
#   geom_histogram() +
#   geom_vline(xintercept = 0)
# 
# # transform into an sf object
# nmk_sf <- nmk_res %>% 
#   left_join(ngrid_lookup) %>% 
#   select(grid_index, statistic, p_value, geometry) %>% 
#   st_as_sf()
# 
# norge_proj <- aea_string(c(13, 65), c(60, 70))
# 
# 
# 
# 
# ggplot() +  
#   geom_sf(data = sf::st_transform(world, norge_proj), fill="gray20", color="gray40") +
#   geom_sf(data = st_transform(nmk_sf, norge_proj), aes(fill = statistic), color="black", size=0.05) + #color = NA)+#
#   coord_sf(xlim = c(-0.6e6, 0.9e6), ylim = c(-1e6, 1e6), expand = FALSE) +
#   labs(fill = "MK estimate") +
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
#   my_map_theme +
#   theme(legend.position="bottom") 
# 
# ggsave(here("figs/modis/rgnd_intensity/10km_grid/gris_map_all_p.png"), width = 6, height = 6, units = "in")
# 
# 
# ggplot() +  
#   geom_sf(data = sf::st_transform(world, norge_proj), fill="gray20", color="gray40") +
#   geom_sf(data = st_transform(nmk_sf, norge_proj) %>% 
#             filter(p_value<0.05), aes(fill = statistic), color="black", size=0.05) + #color = NA)+#
#   coord_sf(xlim = c(-0.6e6, 0.9e6), ylim = c(-0.7e6, 1e6), expand = FALSE) +
#   labs(fill = "MK estimate") +
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
#   my_map_theme +
#   theme(legend.position="bottom") 
# 
# ggsave(here("figs/modis/rgnd_intensity/10km_grid/gris_map_signif.png"), width = 6, height = 6, units = "in")
# 
# 
# 
# # just svalbard
# ggplot() +  
#   geom_sf(data = sf::st_transform(world, norge_proj), fill="gray20", color="gray40") +
#   geom_sf(data = st_transform(nmk_sf, norge_proj), aes(fill = statistic), color="black", size=0.05) + #color = NA)+#
#   coord_sf(xlim = c(-0.2e6, 0.4e6), ylim = c(1.2e6, 1.8e6), expand = FALSE) +
#   labs(fill = "MK estimate") +
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
#   my_map_theme +
#   theme(legend.position="bottom") 
# 
# ggsave(here("figs/modis/rgnd_intensity/10km_grid/sval_map_all_p.png"), width = 6, height = 6, units = "in")
# 
# ggplot() +  
#   geom_sf(data = sf::st_transform(world, norge_proj), fill="gray20", color="gray40") +
#   geom_sf(data = st_transform(nmk_sf, norge_proj) %>% 
#             filter(p_value<0.05), aes(fill = statistic), color="black", size=0.05) + #color = NA)+#
#   coord_sf(xlim = c(-0.2e6, 0.4e6), ylim = c(1.2e6, 1.8e6), expand = FALSE) +
#   labs(fill = "MK estimate") +
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
#   my_map_theme +
#   theme(legend.position="bottom") 
# 
# ggsave(here("figs/modis/rgnd_intensity/10km_grid/sval_map_sig.png"), width = 6, height = 6, units = "in")
# 
# 
# 
# ## same for KAMCHATKA --------------------------------
# 
# 
# # read into R
# kgrids <- st_read(z_local_paths[z_local_paths %>% str_detect("kam5") ])
# kam_proj <- aea_string(c(150, 56), c(50,60))
# 
# 
# 
# mk_kml <- function(grid_ts, proj){
#   grid_data <- grid_ts %>% 
#     st_drop_geometry() %>% 
#     as.data.frame() %>% 
#     mutate(grid_index = str_split_fixed(system_index,"_",2)[,2]) %>% 
#     select(grid_index, year, count, sum) %>% 
#     as_tibble() %>% 
#     mutate(across(.cols = year:sum, as.numeric))
#   
#   # lookup table
#   grid_lookup <- grid_ts %>% 
#     as.data.frame() %>% 
#     mutate(grid_index = str_split_fixed(system_index,"_",2)[,2]) %>% 
#     distinct(grid_index, geometry)
#   
#   
#   mk_res <- grid_data %>%  
#     arrange(year) %>% 
#     group_by(grid_index) %>%
#     nest() %>% 
#     ungroup() %>% 
#     mutate(mk = map(data, ~pull(.x, sum) %>% MannKendall() %>% broom::tidy()) ) %>% 
#     unnest(mk) %>% 
#     janitor::clean_names() 
#   
#   # transform into an sf object
#   mk_sf <- mk_res %>% 
#     left_join(grid_lookup) %>% 
#     select(grid_index, statistic, p_value, geometry) %>% 
#     st_as_sf() %>% 
#     st_transform(proj)
#   return(mk_sf)
# }
# kam_mk_res <- mk_kml(kgrids, kam_proj)
# 
# # plot a histogram
# kam_mk_res %>% 
#   # filter(p_value<0.1) %>%
#   ggplot(aes(statistic)) +
#   geom_histogram() +
#   geom_vline(xintercept = 0)
# 
# 
# ggplot() +  
#   geom_sf(data = st_transform(world, kam_proj), fill="gray20", color="gray40") +
#   geom_sf(data = kam_mk_res, aes(fill = statistic), color="black", size=0.05) + #color = NA)+#
#   coord_sf(xlim = c(-0.1e6, 0.9e6), ylim = c(-1e6, 1e6), expand = FALSE) +
#   labs(fill = "MK estimate") +
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
#   my_map_theme +
#   theme(legend.position="bottom") 
# 
# ggsave(here("figs/modis/rgnd_intensity/10km_grid/kam_map_all_p.png"), width = 6, height = 6, units = "in")
# 
# 
# ggplot() +  
#   geom_sf(data = st_transform(world, kam_proj), fill="gray20", color="gray40") +
#   geom_sf(data = kam_mk_res %>% filter(p_value<0.05), aes(fill = statistic), color="black", size=0.05) + #color = NA)+#
#   coord_sf(xlim = c(-0.1e6, 0.9e6), ylim = c(-1e6, 1e6), expand = FALSE) +
#   labs(fill = "MK estimate") +
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
#   my_map_theme +
#   theme(legend.position="bottom") 
# 
# ggsave(here("figs/modis/rgnd_intensity/10km_grid/kam_map_sig.png"), width = 6, height = 6, units = "in")
