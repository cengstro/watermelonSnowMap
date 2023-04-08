library(tidyverse)
library(here)
library(janitor)
library(sf)
library(lme4)
library(fs)

dir <- here("data/s2_temporal/")

my_read <- function(path, remove_index = TRUE) {
  df <- read_csv(path) 
  df_tidy <- df %>% 
    janitor::clean_names() %>% 
    select(-geo, -subregion_of) %>% 
    rename(region=name)
  if(remove_index == TRUE){
    df_tidy %>% select(-system_index)
  }else{
    df_tidy
  }
}

freq_areas_raw <- my_read(paste0(dir, "freqAreas.csv"))
freq_map <- st_read(here("data/s2_temporal/freqMapPixCountGt100.kml")) %>% mutate(mean = as.numeric(mean))

algal_areas_raw <- my_read(paste0(dir, "algalAreas.csv"))

# mean_rgnd <- read_csv(here("data/s2_temporal/annMaxRgndGlaciermask.csv")) #within the glacier mask
mean_rgnd_algaemask <- read_csv(here("data/s2_temporal/annMaxRgndAlgaemask.csv"))
glacier_area_raw <- dir_ls(here("data/s2_classifier_map/area/per_region_F1"), regexp = "T.kml$") %>%
  map_df(st_read)

doy_max_rgnd <- my_read(paste0(dir, "meanMaxRgndDoyStats.csv"))

# 500 random pts per region per year of peak bloom date, with corresponding elevation, latitude, longitude 
predictor_raw <- read_csv(here("data/s2_temporal/maxRgndDoyPredictorsStrat.csv"))
  

biomass_coefs <- read_csv(here("data/engstrom_etal_2022_model_coefs/abundance_coefs.csv")) # also want the dry biomass coefs
states <- rnaturalearth::ne_states(c("united states of america", "canada"), returnclass = "sf")





# # read in data for snowmelt stats
# melt2020 <- my_read(paste0(dir,"2020algalMeltTotals.csv"))
# gl_melt <-  my_read(paste0(dir,"2020algaMeltPerGlac.csv"))
# mb <- read_csv(paste0(dir,"summer_mass_balance_data.csv")) # summer mass balance data for context
# 
# # algal frac cover per glacier (from 2_algae_per_glacier)
# algal_cover_glacier <- read_csv(here("data/s2_classifier_map/algae_per_glacier/gids_north_america_gt10.csv"))
# 
# 




# tidy ------------------------

## frequency data -------------------------
freq_areas <- freq_areas_raw %>% 
  pivot_longer(cols = starts_with("constant"), 
               names_prefix = "constant_", 
               names_to = "n_yrs", 
               values_to = "area") %>% 
  mutate(area_km2 = area/(1000*1000), 
         n_yrs = as.numeric(n_yrs)+1, 
         .keep = "unused") %>%
  group_by(region) %>% 
  mutate(frac = area_km2 / sum(area_km2)) %>% 
  ungroup()

region_areas <- freq_areas %>% 
  group_by(region) %>% 
  summarise(area_km2 = sum(area_km2))

cont_freq_area <- freq_areas %>% 
  group_by(n_yrs) %>% 
  summarise(area_km2 = sum(area_km2)) %>% 
  mutate(frac = area_km2 / sum(area_km2))


##biomass data ---------------------
glacier_area <- glacier_area_raw %>%
  distinct(region, .keep_all=TRUE) %>% 
  as_tibble() %>% 
  mutate(glacier_km2 = as.numeric(glacier), algae_km2 = as.numeric(algae)) %>%
  select(region, algae_km2, glacier_km2) %>% 
  filter(region %in% c("alaskaRange", "coastNorth", "coastSouth", "interiorNorth", "interiorSouth"))

algal_areas <- algal_areas_raw %>% 
  select(-algal_area_above_snowline_km2) %>% 
  left_join(glacier_area) %>% 
  mutate(algae_glac_ratio = algal_area_km2/glacier_km2)

# import the coefficients from my algal biomass-rgnd linear model (engstrom et al 2022)
toc_slope <- pull(biomass_coefs[4,"estimate"]) # in g/m2
toc_intercept <- pull(biomass_coefs[3,"estimate"]) 
toc_ci <- 1.96*(pull(biomass_coefs[4,"std.error"])) 
tn_slope <- (pull(biomass_coefs[6,"estimate"]))/1000 # convert mg/m2 -> g/m2
tn_intercept <- (pull(biomass_coefs[5,"estimate"]))/1000 
tn_ci <- (1.96*(pull(biomass_coefs[6,"std.error"])))/1000
# for plotting two y-axis
toc_to_tn_ratio <- tn_slope/toc_slope # TN = (Kn/Kc) * TOC

# scale sat data by coefficients to estimate algal biomass per year
# theoretically within the glacier mask should yield same results, but not working
# biomass <- mean_rgnd %>% 
#   clean_names() %>% 
#   rename(region=name) %>% 
#   left_join(glacier_area) %>% 
#   # convert RGND within the glacier mask to units of biomass
#   mutate(tn_g_m2 = tn_slope*mean + tn_intercept ,
#          toc_g_m2 = toc_slope*mean + toc_intercept,
#          toc_m2_ci = mean*toc_ci,
#          tn_m2_ci = mean*tn_ci) %>% 
#   mutate(total_tn_g = tn_g_m2 *glacier_km2*1000*1000,
#          total_toc_g = toc_g_m2 *glacier_km2*1000*1000,
#          total_tn_ci = tn_m2_ci *glacier_km2*1000*1000,
#          total_toc_ci = toc_m2_ci *glacier_km2*1000*1000)

# check that I get the same answer using the algaemask data
biomass <- mean_rgnd_algaemask %>% 
  clean_names() %>% 
  rename(region=name) %>% 
  left_join(glacier_area) %>% 
  mutate(tn_g_m2 = tn_slope*mean + tn_intercept ,
         toc_g_m2 = toc_slope*mean + toc_intercept,
         toc_m2_ci = mean*toc_ci,
         tn_m2_ci = mean*tn_ci) %>% 
  mutate(total_tn_g = tn_g_m2 *algae_km2*1000*1000,
         total_toc_g = toc_g_m2 *algae_km2*1000*1000,
         total_tn_ci = tn_m2_ci *algae_km2*1000*1000,
         total_toc_ci = toc_m2_ci *algae_km2*1000*1000)
  

# # ratio of biomass above and below snowline, per region per year
# biomass_snowline_ratio <- mean_rgnd %>% 
#   select(region, year,algal_rgnd_above_snowline_sum, algal_rgnd_sum) %>% 
#   mutate(frac_above = algal_rgnd_above_snowline_sum/algal_rgnd_sum, .keep="unused")


## DOY covariates data ---------------------

predictors <- predictor_raw %>% 
  mutate(region = case_when(region_id==1~"alaskaRange",
                            region_id==2~"coastNorth",
                            region_id==3~"coastSouth",
                            region_id==4~"interiorNorth",
                            region_id==5~"interiorSouth")) %>% 
  select(region, year, max_rgnd_doy, dem, latitude, longitude)




##  mapping setup ---------
# function for defining Albers equal area projection, with custom center and standard parallels 
make_aea_proj <- function(center_on, parallels){
  lon_0 <- center_on[1]
  lat_0 <- center_on[2]
  lat_1 <- parallels[1]
  lat_2 <- parallels[2]
  paste0("+proj=aea +lat_0=",lat_0," +lon_0=",lon_0," +lat_1=",lat_1," +lat_2=",lat_2,
         " +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +type=crs")
}
north_america_aea_prj <- make_aea_proj(c(-137, 52 ),c(50, 60) )


region_ord <- c("alaskaRange","coastNorth","interiorNorth","coastSouth","interiorSouth")



# ANALYSIS:
# doy of max rgnd ------------------------

doy_max_rgnd %>% 
  arrange(mean)

doy_max_rgnd %>% 
  # weight by algae area
  left_join(region_areas) %>% 
  summarise(mean = weighted.mean(x = mean, w = area_km2), ci = 1.96*weighted.mean(x=std_dev, w = area_km2)) %>% 
  mutate(date = as.Date(mean, '1900-01-01'))
# July 25, +/- 31 d

doy_max_rgnd %>% 
  mutate(ci = std_dev*1.96, date_dummy = as.Date(mean,"1900-01-01"),
         region = fct_relevel(region, rev(region_ord))) %>% 
  ggplot(aes(x = region, y = date_dummy)) +
  geom_point() +
  geom_errorbar(aes(ymin = date_dummy-ci, ymax = date_dummy+ci), width=0.1) +
  facet_grid(rows = vars(year), scales = "free") +
  labs(y = "Date of peak bloom", x = "Region") +
  coord_flip() 

# all regions are significantly different, the sim sample size is so large. 
# # ANOVA with simulated data
# # https://stackoverflow.com/questions/29260139/r-function-to-perform-anova-and-tukeyhsd-from-sample-mean-sd-and-n
# gen_data <- function(means, sds, samplesizes){
#   n.grp <- length(means)
#   grps <- factor(rep(1:n.grp, samplesizes))
#   dat <- lapply(1:n.grp, function(i) {scale(rnorm(samplesizes[i]))*sds[i] + means[i]})
#   y <- do.call(rbind, dat)
#   out <- data.frame(group = grps, y = y)
#   out
# }
# sim_data <- gen_data(means = doy_max_rgnd$mean, sds = doy_max_rgnd$std_dev, samplesizes = doy_max_rgnd$count)
# glimpse(sim_data)
# lm(y~group, data = sim_data) %>% summary()
# # just alaska vs the others
# sim2 <- sim_data %>% 
#   mutate(group = if_else(group=="1", "1", "2")) %>% 
#   mutate(group = fct_relevel(group, "2", "1"))
# lm(y~group, data = sim2) %>% summary()
# # alaska range is 7 days earlier than the average
# 
# # sanity check that there is no difference between the others
# lm(y~group, data = sim_data %>% mutate(group = fct_relevel(group, "2"))) %>% summary()
# # hmmmmmm


## predictors ------------------------

nrow(predictors)

# elevation distribution per region (all four years)
predictors %>% 
  ggplot(aes(dem))+
  geom_histogram()+
  facet_grid(rows =vars(region)) +
  labs(x = "Elevation (m)", y = "Pixel count")
ggsave(here("figs/s2_temporal/north_america_slope_hist.png"), width = 4, height = 5, units = "in")

predictors %>% 
  pivot_longer(cols = c(latitude, longitude) ) %>% 
  ggplot(aes(x = value, y = max_rgnd_doy)) +
  geom_point(alpha = 0.1)+
  facet_grid(cols = vars(name), rows = vars(year), scales = "free") +
  geom_smooth(method = "lm")


# maybe big bloom years have a elevational gradient?
predictors_plus <- predictors %>% 
  left_join(mean_rgnd_algaemask, by = c("region"="name", "year")) 

ggplot(predictors_plus, aes(x = dem, y = max_rgnd_doy, color = mean)) +
  geom_point(alpha = 0.1)+
  geom_smooth(method = "lm") +
  scale_color_viridis_c()

mod1 <- lm(max_rgnd_doy ~ dem + mean, data = predictors_plus)
summary(mod1) # .1 day per 100 m

##elevation vs doy plot ------------------

predictors_plus %>% 
  mutate(date_dummy = as.Date(max_rgnd_doy, "1900-01-01")) %>% 
  ggplot(aes(x = dem, y = date_dummy, color = mean)) +
  geom_point(alpha = 0.1)+
  facet_grid(rows = vars(year), cols = vars(region), scales = "free") +
  geom_smooth(method = "lm", color = "blue") +
  scale_color_distiller(palette = "Reds", direction = 1) +
  labs(x = "Elevation (m)", y = "Date of Peak Bloom", color = "Bloom intensity") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
ggsave(here("figs/s2_temporal/elevation.png"))

mod2 <- lmer(max_rgnd_doy ~ dem + mean + (1|region) + (1|year), data = predictors_plus)
summary(mod2) # 1 day per 100 m, 3.2 days per 0.01 increase in max intensity 
# higher = later
# norther = earlier !!
# easter = later (this is likely a secondary effect of the latitude difference)

mod1 <- lm(max_rgnd_doy ~ latitude + longitude + dem, data = predictors)
mod2 <- lm(max_rgnd_doy ~ dem, data = predictors)
mod3 <- lm(max_rgnd_doy ~ latitude + dem, data = predictors)
mod4 <- lm(max_rgnd_doy ~ longitude + dem, data = predictors)
AIC(mod1, mod2, mod3, mod4) %>% 
  arrange(AIC)
summary(mod1) # report this

# for every 10 deg further north (roughly the difference between Vancouver and Anchorage), peak bloom occurs 5 days earlier 
# for every 10 deg further east (roughly the distance from Vancouver to Calgary), peak bloom occurs 3.6 days earlier
# for every 100 m higher, peak bloom occurs 1 days later (+/-0.15)

# mod6 <- lmer(max_rgnd_doy ~ latitude + longitude + dem + (1|region), data = predictors) # allow intercept to vary by region
# summary(mod6)
# 
# mod7 <- lmer(max_rgnd_doy ~dem + (1|region), data = predictors) 
# summary(mod7)
# mod8 <- lmer(max_rgnd_doy ~dem + (dem|region), data = predictors) 
# AIC(mod6, mod7, mod8) %>% 
#   arrange(AIC)
# 



# frequency stats -------------------
cont_freq_area

# map of frequency
# mean per pixel over 4 years
ggplot() + 
  geom_sf(data = st_transform(states, north_america_aea_prj), fill="gray70") +
  geom_sf(data = st_transform(freq_map, north_america_aea_prj), aes(fill = mean), size = 0.1, color = NA) +
  coord_sf(xlim = c(-1e6, 1.6e6), ylim = c(-0.3e6, 1.4e6), expand = FALSE) +
  theme(panel.background = element_rect(fill = "gray20"),
        panel.grid.major = element_line(colour="gray50", size=0.1)) +
  labs(fill = "n years\nwith bloom") +
  scale_fill_distiller(palette = "RdPu", direction = 1) +
  theme(legend.position="bottom")

ggsave(here("figs/s2_temporal/freqmap.png"), width = 8, height = 8, units="in")

# does this visually overlap with the regions of high % cover?


# biomass plots and stats -------------------------------


pa <- biomass_v2 %>% 
  select(region, year,total_toc_g, total_toc_ci) %>% 
  mutate(total_toc_kg = total_toc_g/1000,
         total_toc_ci = total_toc_ci/1000,
         .keep="unused") %>% 
  ggplot(aes(x = year, color = region, y = total_toc_kg)) +
  geom_point() +
  geom_errorbar(aes(ymin = total_toc_kg-total_toc_ci, ymax = total_toc_kg+total_toc_ci), width = 0.1, alpha = 0.3)+
  geom_line(alpha = 0.5) +
  scale_y_continuous(
    "total kg TOC", 
    sec.axis = sec_axis(~ . * toc_to_tn_ratio, name = "total kg TN")
  ) +
  labs(tag = "A")
pa

# per unit area
pb <- biomass_v2 %>% 
  select(region, year, toc_g_m2, toc_m2_ci) %>% 
  ggplot(aes(x = year, color = region, y = toc_g_m2)) +
  geom_point() +
  geom_errorbar(aes(ymin = toc_g_m2-toc_m2_ci, ymax = toc_g_m2+toc_m2_ci), width = 0.1, alpha = 0.3)+
  geom_line(alpha = 0.5) +
  scale_y_continuous(
    bquote('g TOC'~m^-2), 
    sec.axis = sec_axis(~ . * toc_to_tn_ratio, name = bquote('g TN'~m^-2))
  ) +
  labs(tag = "B")
pb
# combine in a 2 panel plot
ggpubr::ggarrange(pa, pb, nrow = 2)
ggsave(here("figs/s2_temporal/biomass_yr.png"), units = "in", height = 6, width = 6)


# stats 

# totals
biomass %>% 
  select(region, year, contains("total")) %>% 
  mutate(total_toc_kg = total_toc_g/1000,
         total_toc_ci = total_toc_ci/1000,
         total_tn_kg = total_tn_g/1000,
         total_tn_ci = total_tn_ci/1000,
         .keep="unused") %>% 
  group_by(year) %>% 
  summarise(total_toc_kg = sum(total_toc_kg),
            total_toc_ci = sum(total_toc_ci),
            total_tn_kg = sum(total_tn_kg),
            total_tn_ci = sum(total_tn_ci)) %>% 
  arrange(-year)
# year total_toc_kg total_toc_ci total_tn_kg total_tn_ci
#  2022      607344.      110939.      28465.       2444.
# 6.1 ± 1.1 × 10e5 kg TOC
# 2.9 ± 0.2 × 10e4 kg TN


# totals per region
biomass %>% 
  select(region, year, contains("total")) %>% 
  mutate(total_toc_kg = total_toc_g/1000,
         total_toc_ci = total_toc_ci/1000,
         .keep="unused") %>% 
  arrange(-total_toc_kg)


# per area
biomass %>% 
  select(region, year, contains("m2")) %>% 
  arrange(-toc_g_m2)
#   region         year  toc_g_m2 toc_m2_ci tn_g_m2 tn_m2_ci
# interiorSouth  2020     0.677    0.134    0.0289  0.00295 
# 677 ± 134 mg TOC /m2 and 29 ± 3 mg TN /m2






# SCRATCH--------------------------------

##  ELA snowline biomass ratio----------------

biomass_snowline_ratio %>% 
  left_join(algal_areas %>% select(region, year, algal_area_km2)) %>% 
  filter(year !=2022) %>% 
  group_by(year) %>% 
  summarise(frac_above = weighted.mean(frac_above, algal_area_km2)) %>% 
  summarise(mean = mean(frac_above))
#45% above, 55% below

# what is going on in 2022? must be very high snowline relative to the algal bloom
# 

# 
# # Snowmelt analysis ----------------------------
# melt2020 %>% 
#   summarise(sum(algal_error_mmswe), sum(algal_melt_mmswe)) # units actually in L
# # 8.9 +/- 1.7 e9 L snowmelt in 2020 throughout North America
# 
# melt2020 %>% 
#   arrange(-algal_melt_mmswe)
# # 3.3 +/- 0.6 e9 in Interior South in 2020
# 
# gl_select <- gl_melt %>% 
#   select(glac_id, glac_name, algal_melt_mmswe_mean, algal_error_mmswe_mean, db_area)
# 
# gl_select %>% 
#   arrange(-algal_melt_mmswe_mean)
# # 5.8 +/- 1.1 mm swe
# gl_select %>% 
#   filter(algal_melt_mmswe_mean>2) %>% 
#   distinct(glac_id)
# 
# 
# 
# 
# 
# # mass balance compariton ---------------
# 
# mb %>% distinct(glac_name, reference)
# 
# mb %>% 
#   group_by(glac_name) %>% 
#   summarise(mean = mean(summer_mass_balance_mwe)) %>% 
#   inner_join(gl_select) %>% 
#   mutate(frac = (-algal_melt_mmswe_mean /1000) / mean,
#          pct = frac*100) %>% 
#   select(-glac_id) %>% 
#   arrange(-pct) %>% 
#   summarise(mean(pct), 1.96*sd(pct), mean(mean), 1.96*sd(mean),
#             mean(algal_melt_mmswe_mean), 1.96*sd(algal_melt_mmswe_mean))
# 
# # get the algal % cover for these glaciers
# gl_select %>% 
#   filter(glac_name %in% mb$glac_name) %>% 
#   left_join(algal_cover_glacier) %>% 
#   summarise(mean(frac), 1.96*sd(frac))
# # 18 +/- 7 %