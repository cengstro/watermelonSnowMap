#*
#*

library(tidyverse)
library(here)
library(janitor)
library(lubridate)
library(sf)
library(snakecase)
library(flextable)
library(officer)

theme_set(theme_bw())

albedo_coefs <- read_csv(here("data/engstrom_etal_2022_model_coefs/albedo_coefs.csv"))

# per region
algae_region_raw <- read_csv(here("data/s2_temporal/rgndTimeSeriesRegions.csv"))  # contains mean RGND within both algal and glacier masks
insol_region_raw <- read_csv(here("data/s2_temporal/insolTimeSeriesRegions.csv")) 
region_meta_raw <- read_csv(here("data/s2_classifier_map/area/thrsRF1Area.csv"))

# per glacier
algae_glacier_raw <- read_csv(here("data/s2_temporal/rgndTimeSeriesGlaciers.csv"))   # contains mean RGND within both algal and glacier masks
insol_glacier_raw <- read_csv(here("data/s2_temporal/insolTimeSeriesGlaciers.csv"))# mean per glacier (units: Joules/m^2)
glacier_frac <- read_csv(here('data/s2_classifier_map/algae_per_glacier/gids_north_america_gt10.csv'))
glacier_sf <- st_read(here('data/glimsMeta.kml')) 

# per grid cell
algae_grid_raw <- read_csv(here("data/s2_temporal/rgndTimeSeriesGrid.csv"))  
insol_grid_raw <- read_csv(here("data/s2_temporal/insolTimeSeriesGrid.csv"))
grid_sf <- st_read(here("data/s2_temporal/northAmerica20kGrid.kml"))
states <- rnaturalearth::ne_states(c("united states of america", "canada"), returnclass = "sf")

massbal_raw <- read_csv(here("data/s2_temporal/summer_mass_balance_data.csv")) 

# map data
regions <- st_read(here("data/s2_classifier_map/glacierRegionsV15.kml")) 




# tidy --------------

albedo_coef <- albedo_coefs %>% 
  filter(term == "mean_rgnd") %>% 
  pull(estimate)

albedo_ci <- albedo_coefs %>% 
  filter(term=="mean_rgnd") %>% 
  mutate(ci = std.error*1.96) %>% 
  pull(ci)


my_clean_date <- function(df){
  df %>% 
    janitor::clean_names() %>% 
    drop_na() %>% 
    mutate(year = year(date),
           yday = yday(date),
           date_dummy = as.Date(yday,"1900-01-01"))
}


## region --------------------
# the algaemask data is used for removing low pixel obs, and the glaciermask is used for the snowmelt imputation with splines
algae_region <- algae_region_raw %>% 
  my_clean_date() %>% 
  rename_with(~str_remove(.x, "algal_rgnd_")) %>% 
  rename(glacier_mean = mean,
         glacier_std_dev = std_dev,
         glacier_count = count,
         region = name)

insol_region <- insol_region_raw %>%  
  my_clean_date() %>% 
  rename(region = name)

region_meta <- region_meta_raw %>% 
  rename(region = name, algae_km2 = algae, glacier_km2 = glacier, redmask_km2 = redmask) %>% 
  mutate(algal_frac = algae_km2/glacier_km2) %>% 
  filter(region %in% c("alaskaRange", "coastNorth", "coastSouth", "interiorNorth", "interiorSouth"))

 ## glacier ----------------------
algae_glacier <- algae_glacier_raw %>%
  my_clean_date() %>% 
  select(-system_index, -geo) %>% 
  rename_with(~str_remove(.x, "algal_rgnd_")) %>% 
  rename_with(~str_remove(.x, "in_")) %>% 
  rename(glacier_mean = mean,
         glacier_std_dev = std_dev,
         glacier_count = count) %>% 
  relocate(glac_id, date)

insol_glacier <- insol_glacier_raw %>%
  my_clean_date() %>% 
  relocate(glac_id, date) 

glacier_meta <- glacier_sf %>% 
  st_centroid() %>% 
  as_tibble() %>% 
  select(glac_id, db_area, glac_name, geometry) %>% 
  left_join(glacier_frac) %>% 
  mutate(glac_name = if_else(glac_name=="None", NA_character_, glac_name)) %>% 
  rename(glacier_km2 = db_area, algal_frac = frac)


## grid ----------------------------

algae_grid <- algae_grid_raw %>%
  my_clean_date() %>% 
  mutate(system_index = str_remove(system_index, ".*_")) %>% 
  rename_with(~str_remove(.x, "algal_rgnd_")) %>% 
  rename(glacier_mean = mean,
         glacier_std_dev = std_dev,
         glacier_count = count)

insol_grid <- insol_grid_raw %>% 
  my_clean_date() %>% 
  mutate(system_index = str_remove(system_index, ".*_")) 



## massbal ------
massbal <- massbal_raw %>% 
  group_by(glac_name, glac_id, region, range, reference) %>% 
  summarise(minyr = min(year), maxyr = max(year), nyr = n(),
            mean_bal_mwe = mean(summer_mass_balance_mwe),
            sd_bal_mwe = sd(summer_mass_balance_mwe)) %>% 
  ungroup() %>% 
  mutate(minyr = if_else(glac_name=="Peyto Glacier", 2011, minyr),
         maxyr = if_else(glac_name=="Peyto Glacier", 2013, maxyr),
         nyr = if_else(glac_name=="Peyto Glacier", 3, nyr)) %>% 
  relocate(-region, -range, -reference) %>% 
  inner_join(glacier_meta %>% select(-glac_name), by = "glac_id")


# interpolate algae timeseries with splines --------------------------

# like augment, but includes predictions for each day
my_augment <- function(mod){
  data <- mod$data
  ydays <- min(mod$x):max(mod$x) # region specific range of days
  predlist <- predict(mod, ydays)
  tibble(yday = predlist$x, .fitted = predlist$y) %>% 
    left_join(tibble(glacier_mean = data$y, yday = data$x, w = data$w), by = "yday") %>%
    filter(.fitted>=0) # prediction should not be below 0
}

# wrapper for my_augment
tidy_splines <- function(data){
  data %>% 
    mutate(splpreds = map(spl, my_augment) ) %>% 
    unnest(cols = c(splpreds)) %>% 
    ungroup() %>% 
    mutate(date_dummy = as.Date(yday,"1900-01-01")) %>% 
    select(-data, -spl)
}

## region --------------------------------
algae_region_splines <- algae_region %>% 
  group_by(region, year) %>% 
  mutate(nobs =n(), df = 0.35*nobs) %>% 
  nest() %>% 
  mutate(spl = map(data, ~smooth.spline(.x$yday, .x$glacier_mean, .x$algaemask_count, df =.x$df))) %>% 
  tidy_splines()

algae_region_splines %>% 
  group_by(region, year) %>% 
  mutate(pct_clear = (w/max(w, na.rm=TRUE))*100 ) %>% 
  ungroup() %>% 
  ggplot(aes(x = date_dummy, y = glacier_mean, color = pct_clear)) +
  geom_point()+
  facet_grid(cols = vars(year), rows = vars(region), scales = "free") +
  geom_line(aes(y = .fitted), alpha = 0.5, color = "red")+
  scale_color_distiller(palette = "PuBu", direction = 1)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(y = "Snow redness index", x = "Date", color = "% clear")
ggsave(here("figs/s2_temporal/splines.png"), height = 6, width = 8, units = "in")

## glacier ----------------------------------

algae_glacier_splines <- algae_glacier %>% 
  # since the goal is to find maximums, filter to glaciers with at least one high mean RGND observations
  group_by(glac_id, year) %>%
  mutate(yrmax = max(glacier_mean), n = n(), maxdate = yday[which.max(glacier_mean)]) %>% 
  filter(yrmax>0.01, n>=6, maxdate<222) %>%
  nest() %>%
  mutate(spl = map(data, ~smooth.spline(.x$yday, .x$glacier_mean, .x$algaemask_count))) %>% 
  tidy_splines()

set.seed(123)
algae_glacier_splines_sample <- algae_glacier_splines %>% 
  group_by(glac_id, year) %>% 
  nest() %>% 
  ungroup() %>% 
  slice_sample(n=9) %>% 
  unnest()

algae_glacier_splines_sample %>% 
  ggplot(aes(x = date_dummy, y = glacier_mean, color = w)) +
  geom_point()+
  facet_wrap(vars(glac_id, year), scales = "free") +
  geom_line(aes(y = .fitted), alpha = 0.5, color = "red")+
  scale_color_distiller(palette = "PuBu", direction = 1)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(y = "Snow redness index", x = "Date", color = "Weight \n(% clear\nskies)")

### select glaciers ---------------------------------------
algae_massbal_splines <- algae_glacier %>%
  filter(glac_id %in% massbal$glac_id) %>%
  # since the goal is to find maximums, filter to glaciers with at least one high mean RGND observations
  group_by(glac_id, year) %>%
  mutate(yrmax = max(glacier_mean), n = n(), maxdate = yday[which.max(glacier_mean)], df = 0.55*n) %>% 
  filter(n>3) %>%
  nest() %>%
  mutate(spl = map(data, ~smooth.spline(.x$yday, .x$glacier_mean, .x$algaemask_count, df =.x$df ))) %>% 
  tidy_splines()


algae_massbal_splines %>% 
  left_join(glacier_meta) %>% 
  ggplot(aes(x = date_dummy, y = glacier_mean, color = w)) +
  geom_point()+
  facet_grid(rows = vars(glac_name), cols = vars(year)) +
  geom_line(aes(y = .fitted), alpha = 0.5, color = "red")+
  scale_color_distiller(palette = "PuBu", direction = 1)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(y = "Snow redness index", x = "Date", color = "Weight \n(% clear\nskies)")


## grid -----------------

algae_grid_splines <- algae_grid %>% 
  group_by(system_index, year) %>%
  mutate(yrmax = max(glacier_mean), n = n(), maxdate = yday[which.max(glacier_mean)]) %>% 
  filter(yrmax>0.002, n>=6, maxdate<222) %>%
  nest() %>%
  mutate(spl = map(data, ~smooth.spline(.x$yday, .x$glacier_mean, .x$algaemask_count))) %>% 
  tidy_splines()
algae_grid_splines %>% count(system_index, year) %>% nrow()

set.seed(123)
algae_grid_splines_sample <- algae_grid_splines %>% 
  group_by(system_index, year) %>% 
  nest() %>% 
  ungroup() %>% 
  slice_sample(n=9) %>% 
  unnest()

algae_grid_splines_sample %>% 
  ggplot(aes(x = date_dummy, y = glacier_mean, color = w)) +
  geom_point()+
  facet_wrap(vars(system_index, year), scales = "free") +
  geom_line(aes(y = .fitted), alpha = 0.5, color = "red")+
  scale_color_distiller(palette = "PuBu", direction = 1)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(y = "Snow redness index", x = "Date", color = "Weight \n(% clear\nskies)")



# snowmelt -----------------
latent_heat_snow <- 3.3*10^5 #J/kg

estimate_melt <- function(algae_spline, insol){
  algae_spline %>%
    inner_join(insol) %>% # no insol data for 2022
    mutate( delta_albedo = .fitted*albedo_coef, 
            delta_albedo_ci = .fitted*albedo_ci,
            algal_energy = -1*delta_albedo * insol, # units: J/m2
            algal_energy_ci = delta_albedo_ci * insol,
            algal_drf = -1 * delta_albedo * srad, # units: w/m2
            algal_drf_ci = delta_albedo_ci *srad,
            # convert from energy absorbed by snowpack due to algal albedo
            # to mm SWE by dividing by latent heat of snowmelt
            # (J/m2) / (J/kg) = kg/m2, 
            # (kg/m2) * (1m3/1000kg water) * (1000 mm / 1 m) = mm SWE
            swe_mm = algal_energy / latent_heat_snow, # see above
            swe_mm_ci = algal_energy_ci / latent_heat_snow)
}


## region -------------------------

melt_region <- algae_region_splines %>% 
  filter(year<2022) %>% # 2022 DAYMET data not avail on GEE
  estimate_melt(insol_region)

# date of peak runoff
melt_region %>% 
  group_by(year, region) %>% 
  summarise(peak_melt_date = date_dummy[.fitted == max(.fitted)]) %>% 
  ungroup() %>% 
  summarise(mean = mean(peak_melt_date))
# mean date of peak melt on July 19

# plot daily radiative forcing
melt_region %>% 
  ggplot(aes(x = date_dummy, y = algal_drf)) +
  geom_line()+
  facet_grid(cols = vars(year), rows = vars(region)) +
  geom_ribbon(aes(ymin = algal_drf - algal_drf_ci, ymax = algal_drf + algal_drf_ci), alpha = 0.3) +
  labs(x = "Date") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_y_continuous(
    bquote('Daily radiative forcing over glacier area'~(Wm^-2)), 
    sec.axis = sec_axis(~ . / latent_heat_snow, name = bquote('mm SWE'~d^-1))
  )
ggsave(here("figs/s2_temporal/region_drf.png"), width = 8, height = 8, units = "in")

# plot melt per day in L per region
melt_region %>% 
  left_join(region_meta %>% select(region, glacier_km2)) %>% 
  # kg/m2 = mm SWE = L H2O / m2
  mutate(l_total = swe_mm * glacier_km2 * 1000*1000,
         l_ci =swe_mm_ci * glacier_km2 * 1000*1000 ) %>% 
  ggplot(aes(x = date_dummy, y = l_total)) +
  geom_line()+
  facet_grid(cols = vars(year), rows = vars(region)) +
  geom_ribbon(aes(ymin = l_total - l_ci, ymax = l_total + l_ci), alpha = 0.3) +
  labs(y = "Algal-derived snowmelt (L)", x  = "Date") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(here("figs/s2_temporal/daily_melt_L.png"), width = 8, height = 8, units = "in")

# summary table per year & region
melt_region_litres <- melt_region %>%
  group_by(region, year) %>% 
  summarise(swe_mm  = sum(swe_mm, na.rm=TRUE), # missing insol data for a couple of days in early July in some years, due to variation in the 5 day mosaics
            swe_mm_ci = sum(swe_mm_ci, na.rm=TRUE)) %>% 
  ungroup() %>% 
  left_join(region_meta %>% select(region, glacier_km2)) %>% 
  # kg/m2 = mm SWE = L H2O / m2
  mutate(l_total = swe_mm * glacier_km2 * 1000*1000,
         l_ci =swe_mm_ci * glacier_km2 * 1000*1000 ) %>% 
  arrange(-l_total)

melt_region_litres %>% 
  group_by(year) %>% 
  summarise(sum = sum(l_total), ci =sum(l_ci)) %>% 
  filter(sum==max(sum))
# 3.3e11 ± 1.3e11 L in North America in 2020

# export a formatted table
e_notation <- function(num) formatC(num, format = "e", digits = 0) %>% str_remove("\\+")

snowmelt_table <- melt_region_litres %>% 
  mutate(swe_mm = round(swe_mm, digits = 1), swe_mm_ci = round(swe_mm_ci, digits = 1),
         glacier_km2 = round(glacier_km2, digits = 0),
         l_total = e_notation(l_total), l_ci = e_notation(l_ci),
         region  = snakecase::to_title_case(region)) %>%
  mutate(swe_mm = paste(swe_mm, "±", swe_mm_ci), 
         l_total = paste(l_total, "±", l_ci),
         .keep = "unused")
colnames(snowmelt_table) <- c("Region","Year", "Algal melt (mm SWE ± 95% CI)", "Glacier area (km2)", "Algal melt (L ± 95% CI)")

snowmelt_table %>% 
  mutate(Year = as.character(Year)) %>% 
  flextable() %>%
  save_as_docx( path = here("figs/s2_temporal/snowmelt.docx"),
                sect_properties = prop_section(
                  page_size(orient = "landscape",
                    width = 8.3, height = 11.7
                  )
                ))



## glacier -----

melt_glacier <- estimate_melt(algae_glacier_splines, insol_glacier)

melt_glacier_yrly <- melt_glacier %>% 
  group_by(glac_id, year) %>% 
  summarise(swe_mm  = sum(swe_mm, na.rm=TRUE), # missing insol data for a couple of days in early July in some years, due to variation in the 5 day mosaics
            swe_mm_ci = sum(swe_mm_ci, na.rm=TRUE)) %>% 
  ungroup() %>% 
  arrange(-swe_mm) %>% 
  left_join(glacier_meta)
melt_glacier_yrly
# 62 mm SWE across the surface of a 1.15 km2 glacier in Interior South, with 24% algal cover

glacier_meta %>% arrange(-algal_frac)

# check that the annual spline curve is reasonable in years of high melt
yr_maxmelt <- deframe(melt_glacier_yrly[1,1:2])

algae_glacier_splines %>%
  filter(glac_id == names(yr_maxmelt), year == yr_maxmelt) %>% 
  ggplot(aes(x = date_dummy, y = glacier_mean, color = w)) +
  geom_point()+
  facet_wrap(vars(glac_id, year), scales = "free") +
  geom_line(aes(y = .fitted), alpha = 0.5, color = "red")+
  scale_color_distiller(palette = "PuBu", direction = 1)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(y = "Snow redness index", x = "Date", color = "Weight \n(% clear\nskies)")


# plot the locations of the top 10 watermelon snowmelt glaciers
top10 <- melt_glacier_yrly %>% slice_head(n=10)

aea_string <- function(center_on, parallels){
  lon_0 <- center_on[1]
  lat_0 <- center_on[2]
  lat_1 <- parallels[1]
  lat_2 <- parallels[2]
  paste0("+proj=aea +lat_0=",lat_0," +lon_0=",lon_0," +lat_1=",lat_1," +lat_2=",lat_2,
         " +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +type=crs")
}
america_proj <- aea_string(c(-135, 57),c(49, 63))

my_map_theme <-   theme(panel.background = element_rect(fill = "gray100", colour = NA), # make the ocean near black, 
                        panel.grid.major = element_blank(),
                        axis.ticks.y=element_blank(),
                        axis.text.y=element_blank(),
                        axis.text.x=element_blank(),
                        axis.ticks.x=element_blank() )

# where are these top 25 glacier located
ggplot() +  
  geom_sf(data = st_transform(states, america_proj), fill="gray90", color="gray40") + 
  geom_sf(data = glacier_sf %>% 
            select(-glac_name) %>% 
            inner_join(top10, by = "glac_id") %>% 
            st_transform(america_proj) %>% 
            st_centroid() %>% 
            arrange(swe_mm),
          aes(color = swe_mm, size = as.numeric(glacier_km2)) )+
  coord_sf(xlim = c(-1.3e6, 1.5e6), ylim = c(-1e6, 0.9e6), expand = FALSE) + 
  scale_color_distiller(palette = "YlGnBu", direction = 1, limits = c(0,65)) +
  labs(color = "Snowmelt\n(mm SWE)", size = "Glacier area\n(km2)") +
  my_map_theme
ggsave(here("figs/s2_temporal/top10_glacier_snowmelt.png"))

# which years
ggplot() +  
  geom_sf(data = st_transform(states, america_proj), fill="gray90", color="gray40") + 
  geom_sf(data = glacier_sf %>% 
            select(-glac_name) %>% 
            inner_join(top10, by = "glac_id") %>% 
            st_transform(america_proj) %>% 
            st_centroid(),
          aes(color = as.factor(year)) )+
  coord_sf(xlim = c(-1.3e6, 1.5e6), ylim = c(-1e6, 0.9e6), expand = FALSE) + 
  labs(color = "Year") +
  my_map_theme

## vowell catamount compare w Engstrom et al 2022 ---------
melt_glacier_yrly %>% 
  filter(glac_name %in% c("Catamount Glacier", "Vowell Glacier"), year==2020) 
# 4 cm SWE Vowell Glacier    
# 3 cm SWE Catamount Glacier   

# Engstrom et al estimates: 
# Vowell: 4.2 cm SWE
# Catamount: 4.5 cm SWE

# current study is -5% and -33%


algae_glacier_splines %>%
  left_join(glacier_meta) %>% 
  filter(glac_name %in% c("Catamount Glacier", "Vowell Glacier"), year==2020) %>% 
  ggplot(aes(x = date_dummy, y = glacier_mean, color = w)) +
  geom_point()+
  facet_grid(rows = vars(glac_name), cols = vars(year)) +
  geom_line(aes(y = .fitted), alpha = 0.5, color = "red")+
  scale_color_distiller(palette = "PuBu", direction = 1)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(y = "Snow redness index", x = "Date", color = "Weight \n(% clear\nskies)")


## grid ---------------------------
melt_grid <- estimate_melt(algae_grid_splines, insol_grid)

grid_glac_area <- grid_sf %>% 
  as_tibble() %>% 
  select(system_index, algae_km2, glacier_km2) %>% 
  mutate(across(contains("km2"), as.numeric)) %>% 
  expand_grid(year = 2019:2021) # add year, so that zero years still show up on the map

# annotate with glacier area
melt_grid_yrly <- melt_grid %>% 
  full_join(grid_glac_area) %>% 
  group_by(system_index, year, glacier_km2) %>% 
  # total SWE in mm over the glaciated area, per year
  summarise(swe_mm  = sum(swe_mm, na.rm=TRUE),
            swe_mm_ci = sum(swe_mm_ci, na.rm=TRUE)) %>% 
  ungroup() %>% 
  # kg/m2 = mm SWE = L H2O / m2
  mutate(l_total = swe_mm * glacier_km2 * 1000*1000,
         l_ci =swe_mm_ci * glacier_km2 * 1000*1000 ) %>% 
  replace_na(list(swe_mm = 0)) # so zero years still are plotted



# convert to SF
melt_grid_yrly_sf <- grid_sf %>% 
  select(-glacier_km2) %>% 
  inner_join(melt_grid_yrly, multiple = "all") %>% 
  st_transform(america_proj)

melt_grid_yrly_sf %>% arrange(-swe_mm) # one value =75mm, the secondhighest is 47

# theme just for the zoomed in maps
my_map_theme <-   theme(panel.background = element_rect(fill = "gray100", colour = NA), # make the ocean near black, 
                        panel.grid.major = element_blank(),
                        axis.ticks.y=element_blank(),
                        axis.text.y=element_blank(),
                        axis.text.x=element_blank(),
                        axis.ticks.x=element_blank() )

# # Total snowmelt
# ggplot() +  
#   geom_sf(data = st_transform(states, america_proj), fill="gray40", color="gray70") + 
#   geom_sf(data = melt_grid_yrly_sf,
#           # %>%  filter(year==2020),
#           aes(fill = l_total), color=NA, size=0.02) +
#   coord_sf(xlim = c(-1.3e6, 1.5e6), ylim = c(-1e6, 0.9e6), expand = FALSE) + 
#   scale_fill_distiller(palette = "YlGnBu", direction = 1) +
#   labs(fill = "Snowmelt (L)") +
#   my_map_theme +
#   facet_wrap(vars(year)) 
# 

# melt per glacier area
ggplot() +  
  geom_sf(data = st_transform(states, america_proj), fill="gray100", color="gray40") + 
  geom_sf(data = melt_grid_yrly_sf %>% 
            # filter(year==2020) %>% 
            st_centroid() %>% 
            mutate(swe_mm = if_else(swe_mm>50,50,swe_mm)) %>% # set plotmax-- the single value >75 will be lumped
            arrange(swe_mm), # plot higher melt on top
          aes(color = swe_mm, size = glacier_km2, alpha = swe_mm)) +
  coord_sf(xlim = c(-1.3e6, 1.5e6), ylim = c(-1e6, 0.9e6), expand = FALSE) + 
  scale_color_distiller(palette = "PuBu", direction = 1) +
  labs(color = "Mean snowmelt\n(mm SWE)", size = "Glacier area\n(km2)") +
  my_map_theme +
  facet_grid(rows = vars(year)) 
  # theme(legend.position = "bottom")
ggsave(here("figs/s2_temporal/total_melt.pdf"), width = 7, height = 11, units = "in")



# summer mass balance comparison --------------------------
melt_massbal <- estimate_melt(algae_massbal_splines, insol_glacier)

melt_massbal_yrly <- melt_massbal %>% 
  group_by(glac_id, year) %>% 
  summarise(swe_mm  = sum(swe_mm, na.rm=TRUE), # missing insol data for a couple of days in early July in some years, due to variation in the 5 day mosaics
            swe_mm_ci = sum(swe_mm_ci, na.rm=TRUE)) %>% 
  ungroup() %>% 
  arrange(-swe_mm) %>% 
  left_join(glacier_meta)
melt_massbal_yrly


melt_table_select_glaciers <- melt_massbal_yrly %>% 
  filter(year==2020) %>% 
  left_join(massbal %>% select(-glac_name, -algal_frac, -glacier_km2, -geometry), by = c("glac_id")) %>% 
  mutate(mean_bal_cmwe = mean_bal_mwe*100,
         sd_bal_cmwe = sd_bal_mwe*100,
         swe_cm = swe_mm/10, 
         swe_cm_ci = swe_mm_ci/10,
         lon = st_coordinates(geometry)[,1] %>% round(digits = 3)%>% as.character(), 
         lat = st_coordinates(geometry)[,2] %>% round(digits = 3) %>% as.character(),
         glac_name = str_remove(glac_name," Glacier"),
         glacier_km2 = as.numeric(glacier_km2),
         algal_frac = algal_frac*100,
         percent_of_total= (swe_cm/mean_bal_cmwe * -1 )*100) %>% 
  select(glac_name, glac_id, region, lon, lat, glacier_km2, minyr, maxyr, mean_bal_cmwe, sd_bal_cmwe, algal_frac, swe_cm, swe_cm_ci, percent_of_total, reference) %>% 
  arrange(-swe_cm) %>% 
  mutate(across(where(is.numeric), ~round(.x, digits =1))) %>% 
  mutate(algal_frac = round(algal_frac),
         glacier_km2 = round(glacier_km2)) %>% 
  unite(massbal_years, minyr, maxyr, sep="–") %>% 
  unite(coords, lon, lat, sep = ", ") %>% 
  unite(algal_melt_swe_cm, swe_cm, swe_cm_ci, sep = " ± ") %>% 
  relocate(reference, .after = sd_bal_cmwe)
view(melt_table_select_glaciers)

colnames(melt_table_select_glaciers) <- c("Glacier", "GLIMS ID", "Region", "Coords", "Glac area (km2)", 
                                          "Massbal yrs", "Mean summer massbal (cm SWE)", "Std. dev. summer massbal", "Reference",
                                          "Algal % cover", "Algal melt 2020 (cm SWE)", "Algal melt % of total")



melt_table_select_glaciers %>% 
  flextable() %>% 
  save_as_docx( path = here("figs/s2_temporal/select_glacier_snowmelt.docx"))


# SCRATCH----------------------------
## QC -------------------
minobs <- 6 # min obs per summer
glacier_algae3 <- glacier_algae2 %>% 
  drop_na(rgnd_algalmask) %>% 
  group_by(glac_id) %>% 
  mutate(max = max(count_algalmask), pct_pix_cover = count_algalmask/max) %>% 
  ungroup() %>%
  group_by(glac_id, year) %>% 
  mutate(nobs =n()) %>% 
  filter(nobs>minobs) %>% 
  ungroup()


set.seed(123)
glacier_algae3 %>% 
  group_by(glac_id, year) %>% 
  nest() %>% 
  ungroup() %>% 
  slice_sample(n= 9) %>% 
  unnest() %>% 
  ggplot(aes(x = date_dummy, y = rgnd_glacier, color = pct_pix_cover)) +
  facet_wrap(vars(glac_id, year)) +
  geom_line()+
  geom_point() +
  scale_color_distiller(palette = "PuBu", direction = 1)




# compare snowmelt estimates with my previous study (engstrom et al 2022)
algae5 %>% 
  left_join(my_glims2) %>%
  filter(glac_name %>% str_detect("Cat|Vow"), year==2020) %>%
  group_by(glac_name) %>% 
  filter(delta_albedo == min(delta_albedo))
# new vs old

 
  


# plot raw insol data ---------------------
insol2 %>% 
  ggplot(aes(x = date_dummy, y = insol)) +
  facet_grid(rows = vars(region), cols = vars(year), scales = "free_y") +
  geom_line()+
  geom_point() 

count(insol2, region, year)

set.seed(123)
glacier_insol2 %>% 
  group_by(glac_id, year) %>% 
  nest() %>% 
  ungroup() %>% 
  slice_sample(n= 9) %>% 
  unnest() %>% 
  ggplot(aes(x = date_dummy, y = insol)) +
  facet_wrap(vars(glac_id, year)) +
  geom_line()+
  geom_point() 


# remove low pixel obs ---------------
## regional dataset
colnames(algae2)


algae3 <- algae2 %>% 
  drop_na(rgnd_algalmask) %>% 
  filter(count_algalmask>20) %>% # ok if the glacierarea is masked
  group_by(region) %>% 
  mutate(max = max(count_algalmask), pct_pix_cover = count_algalmask/max) %>% 
  ungroup() %>% 
  group_by(region, year) %>% 
  ungroup()

algae3 %>% filter(count_algalmask<100) %>% select(region)

# Plot the raw data, colored by the number of cloudfree pixels in the algal bloom mask
algae3 %>% 
  ggplot(aes(x = date_dummy, y = rgnd_glacier, color = pct_pix_cover)) +
  facet_grid(rows = vars(region), cols = vars(year), scales = "free_y") +
  geom_line()+
  geom_point() +
  scale_color_distiller(palette = "PuBu", direction = 1)




rgnd_grid2 %>% 
  drop_na(rgnd_algalmask) %>% 
  filter(count_algalmask>20) %>% # ok if the glacierarea is masked
  group_by(region) %>% 
  mutate(max = max(count_algalmask), pct_pix_cover = count_algalmask/max) %>% 
  ungroup() %>% 
  group_by(region, year) %>% 
  ungroup()



  
