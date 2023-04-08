# use the thresh2 data once it runs (old stuff with 0.17)

library(tidyverse)
library(here)
library(janitor)
library(sf)
library(fs)
library(Kendall)
library(trend)
library(broom)
library(lubridate)
library(ggrepel)
library(lme4)

theme_set(theme_bw())


intensity_s2_raw <- read_csv(here("data/s2_temporal/annMaxRgnd.csv")) #within algae mask included
intensity_modis_raw <- read_csv(here("data/modis/meanMaxRgndRegionTerraDay100.csv")) # within Algae mask

doy_s2_raw <- read_csv(here("data/s2_temporal/meanMaxRgndDoyStats.csv"))# within per year classification mask
doy_modis_raw <- read_csv(here("data/modis/meanMaxRgndDoyStats100.csv"))# within Algae mask, excluding 0 pixels

s2ts_raw <- read_csv(here("data/s2_temporal/rgndTimeSeriesRegions.csv"))  # includes within algae mask
modts8_raw <- read_csv(here("data/modis/timeSeries8.csv")) # within the algae mask

# tidy--------------------------

my_clean <- function(df){
  df %>% 
    clean_names() %>% 
    rename(region = name) %>% 
    drop_na()
}

## intensity --------------------
intensity_s2 <- intensity_s2_raw %>% 
  select(name, year, contains("ALGAEMASK")) %>% 
  rename_with(~str_remove(.x, "ALGAL_RGND_ALGAEMASK_")) %>% 
  my_clean() %>% 
  rename(rgnd = mean)

intensity_modis <- intensity_modis_raw %>% 
  my_clean()%>% 
  rename(rgnd = mean)

intensity_long <- intensity_s2 %>% 
  add_column(platform = "s2") %>% 
  bind_rows(intensity_modis %>% add_column(platform = "modis"))

intensity_wide <- intensity_long %>% 
  select(-count, -std_dev) %>% 
  pivot_wider(names_from = platform, values_from = rgnd) %>% 
  drop_na()

## peak bloom date --------------------------

doy_s2 <- my_clean(doy_s2_raw) %>%
  rename(doy = mean) %>% 
  select(region, year, doy, std_dev, count) 
  
doy_modis <- my_clean(doy_modis_raw) %>% 
  rename(doy = mean) %>% 
  select(region, year, doy, std_dev, count) 

doy_long <- doy_s2 %>% 
  add_column(platform = "s2") %>% 
  bind_rows(doy_modis %>% add_column(platform = "modis")) %>% 
  mutate(date_dummy = as.Date(doy, "1900-01-01")) %>% 
  relocate(platform)

doy_wide <- doy_long %>% 
  select(-count, -std_dev, -doy) %>% 
  pivot_wider(names_from = platform, values_from = date_dummy) %>% 
  drop_na()

## time series ------------------
s2ts <- s2ts_raw %>% 
  select(name, date, contains("ALGAEMASK")) %>% 
  rename_with(~str_remove(.x, "ALGAL_RGND_ALGAEMASK_")) %>% 
  my_clean()

modts8 <- modts8_raw %>% 
  my_clean()

ts_long <- s2ts %>% 
  add_column(platform = "s2") %>% 
  bind_rows(modts8 %>% add_column(platform = "modis")) %>% 
  mutate(year = year(date), 
         yday =yday(date),
         date_dummy = as.Date(yday,"1900-01-01")) %>% 
  group_by(region, platform) %>% 
  mutate(percent_count = count / max(count)) %>% 
  ungroup() %>% 
  filter(year>=2019)
  
# ts_wide
  


# intensity -------------------------------------

intensity_wide %>% 
  ggplot(aes(y = modis, x = s2, color = region)) +
  geom_point() +
  coord_fixed(ratio=1) +
  scale_y_continuous(minor_breaks = seq(0,0.04, by = 0.01), breaks = seq(0,0.04, by = 0.01),limits = c(0,0.04))+
  scale_x_continuous(minor_breaks = seq(0,0.04, by = 0.01), breaks = seq(0,0.04, by = 0.01),limits = c(0,0.04)) +
  labs(y= "MODIS SRI", x = "Sentinel-2 SRI") +
  geom_smooth(method = "lm", aes(group = region), se=FALSE)
# the intercept is above zero, meaning MODIS had a higher "baseline" 
# the slope is less than 1:1, meaning MODIS is less sensitive to changes in RGND
ggsave(here("figs/modis/modis_v_s2_maxintensity_color_by_region.png"))

intensity_lmer1 <- lmer(modis~s2+(1|region), data = intensity_wide)
summary(intensity_lmer1) # slope = 0.53 ± 0.26 (95%CI), intercept = 0.019
intensity_lm <-  lm(modis~s2, data = intensity_wide)
summary(intensity_lm) # slope = 0.68 ± 0.37, intercept = 0.018
# LM overestimates the slopes slightly

# not colored by region
s2cor_panel_a <- ggplot(intensity_wide, aes(y = modis, x = s2)) +
  geom_point() +
  coord_fixed(ratio=1) +
  scale_y_continuous(minor_breaks = seq(0,0.04, by = 0.01), breaks = seq(0,0.04, by = 0.01),limits = c(0,0.04))+
  scale_x_continuous(minor_breaks = seq(0,0.04, by = 0.01), breaks = seq(0,0.04, by = 0.01),limits = c(0,0.04)) +
  labs(y= "MODIS max. annual SRI", x = "Sentinel-2 max. annual SRI", tag = "A")# +
  # geom_smooth(method = "lm")
s2cor_panel_a

cor(intensity_wide$s2, intensity_wide$modis) # r = 0.65

# check out which years are outliers for further investigation
ggplot(intensity_wide, aes(x = s2, y = modis, color = as.factor(year))) +
  geom_point() +
  labs(y = "MODIS max. annual SRI", x = "Sentinel-2 max. annual SRI") +
  facet_wrap(vars(region)) 
#* modis highly overestimates in:
#* interior south 2021
#* coast south 2019

cc <- intensity_long %>% 
  group_by(platform, year) %>% 
  summarise(rgnd = weighted.mean(rgnd, count)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = platform, values_from = rgnd) %>% 
  drop_na()
cor(cc$modis, cc$s2) # 0.59


# DOY --------------------------

s2cor_panel_b <- ggplot(doy_wide, aes(x = s2, y = modis)) +
  geom_point() +
  coord_fixed(ratio = 1) +
  scale_x_date(limits =c(as.Date("1900-07-01"), as.Date("1900-08-15"))) +
  scale_y_date(limits =c(as.Date("1900-07-01"), as.Date("1900-08-15"))) +
  labs(y = "MODIS date max. SRI", x = "Sentinel-2 date max. SRI", tag = "B") 
s2cor_panel_b
# the intercept is above zero, meaning MODIS had a later "baseline" 
# the slope is less than 1:1, meaning MODIS is less sensitive to changes in RGND
# Ie modis estimates are more closely clustered around Aug 1, 
# whereas S2 estimates are distributed into earlier into July
# likely because MODIS is picking up on noise ie non algal pixels later in the season

doy_wide_num <- doy_long %>% 
  select(-count, -std_dev, -date_dummy) %>% 
  pivot_wider(names_from = platform, values_from = doy) %>% 
  drop_na()

cor(doy_wide_num$s2, doy_wide_num$modis) # r = 0.81

# colored by region
ggplot(doy_wide, aes(y = modis, x = s2, color = region)) +
  geom_point() +
  coord_fixed(ratio = 1) +
  scale_x_date(limits =c(as.Date("1900-07-01"), as.Date("1900-08-15"))) +
  scale_y_date(limits =c(as.Date("1900-07-01"), as.Date("1900-08-15"))) +
  labs(y = "MODIS date max. SRI", x = "Sentinel-2 date max. SRI", tag = "B") +
  geom_smooth(method  ="lm", se=FALSE)
  # geom_smooth(method = "lm")
  # geom_errorbar(aes(ymin = date_s2 - s2_ti, ymax = date_s2 + s2_ti), width = 0.1)
# s2 predicts a bit earlier, perhaps because of the lower s2 period (if the peak happened Aug 1, and s2 got an image Aug3, s2 peak would be later)
# although the predictions in high bloom years are roughly comparable

# are the years with early onset low or high intensity?
doy_long %>% 
  select(platform, region, year, doy) %>% 
  filter(year>2018) %>% 
  left_join(select(intensity_wide, region, year, s2) %>% rename(rgnd =s2)) %>% 
  pivot_wider(names_from = platform, values_from = doy) %>% 
  ggplot(aes(y = modis, x = s2, color = rgnd)) +
  geom_point() +
  coord_fixed(ratio = 1) +
  facet_wrap(vars(region))
# low years have earlier date (makes sense, likely melt out earlier
# why then did our MODIS find a cor between intensity and EARLIER date?
# perhaps the MODIS relationship is correlation not causation: we would expect an earlier
# melt date due to climate change, and higher intensity due to more saturated snowpack. 

# combined s2 cor plot-------------------
ggpubr::ggarrange(s2cor_panel_a, s2cor_panel_b, ncol =2)
ggsave(here("figs/modis/s2cor.png"))


# time series --------------------------------

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


ts_splines <- ts_long %>% 
  group_by(platform, region, year) %>% 
  mutate(nobs =n(), df = 0.35*nobs) %>% 
  nest() %>% 
  mutate(spl = map(data, ~smooth.spline(.x$yday, .x$mean, .x$percent))) %>% 
  tidy_splines()

ts_splines %>% 
  ggplot(aes(x = date_dummy, y = glacier_mean, color = platform)) +
  geom_point(alpha = 0.5)+
  facet_grid(cols = vars(year), rows = vars(region), scales = "free") +
  geom_line(aes(y = .fitted), alpha = 0.5) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(y = "Snow redness index", x = "Date", color = "Platform")
































# SCRATCH ------------------------------------

## timeseries --------------------

#### interpolate  --------------------------
df_ratio <- 0.35 

# like augment, but includes predictions for each day
my_augment <- function(mod){
  data <- mod$data
  ydays <- min(mod$x):max(mod$x) # region specific range of days
  predlist <- predict(mod, ydays)
  tibble(yday = predlist$x, .fitted = predlist$y) %>% 
    left_join(tibble(mean = data$y, yday = data$x, w = data$w), by = "yday") %>%
    filter(.fitted>=0) # prediction should not be below 0
}
# predict_tbl(splmod) %>%
#   ggplot(aes(x = yday, y = rgnd_glacier, color = w)) +
#   geom_point(color = "red") +
#   geom_point(data = tt)+
#   scale_color_distiller(palette = "PuBu", direction = 1)


ts_splines <- ts %>% 
  group_by(platform, name, year) %>% 
  mutate(nobs =n(), df = df_ratio*nobs) %>% 
  nest() %>% 
  mutate(spl = map(data, ~smooth.spline(.x$yday, .x$mean, .x$percent_count, df =.x$df)),  
         splpreds  = map(spl, my_augment) ) %>% 
  unnest(cols = c(splpreds)) %>% 
  ungroup() %>% 
  mutate(date_dummy = as.Date(yday,"1900-01-01")) %>% 
  select(-data, -spl)

ts_splines %>% 
  ggplot(aes(x = date_dummy, y = mean, color = platform)) +
  geom_point()+
  facet_grid(rows = vars(region), cols = vars(year), scales = "free_y") +
  geom_line(aes(y = .fitted), alpha = 0.5, color = "red")+
  scale_color_distiller(palette = "PuBu", direction = 1)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(y = "Snow redness index", x = "Date", color = "Weight \n(% clear\nskies)")

ts %>% count(platform, year)


ts %>% 
  ggplot(aes(x = date_dummy, y = mean, color = platform))+
  geom_point(data = ts %>% filter(platform=="modis"))+
  geom_line(data = ts %>% filter(platform=="s2"))+
  facet_grid(rows = vars(name), cols = vars(year), scales = "free") +
  geom_smooth(data = ts %>% filter(platform=="modis"))


ts %>% 
  filter(platform=="modis") %>% 
  ggplot(aes(count))+
  geom_histogram()+
  facet_grid(cols = vars(name), rows = vars(year), scales = "free")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(x = "pixel count", y = "n images")

## spline modis --------------


ggplot(modts2, aes(x = date_dummy, y = rgnd)) +
  geom_line()+
  facet_grid(cols = vars(region), rows = vars(year))


# mydoy <- min(modts2$yday):max(modts2$yday)
# tt <- modts2 %>% filter(year==2018, region=="interiorSouth")
# splmod <- smooth.spline(x= tt$yday, y = tt$rgnd, w = tt$count)
# splmod$x
# predict(splmod, mydoy)

predict_tbl <- function(mod){
  data <- mod$data
  ydays <- min(mod$x):max(mod$x) # region specific range of days
  predlist <- predict(mod, ydays)
  tibble(x = predlist$x, .fitted = predlist$y) %>% 
    left_join(tibble(rgnd = data$y, x = data$x, w = data$w), by = "x") %>% 
    rename(yday = x) %>% 
    filter(.fitted>=0) # prediction should not be below 0
}
# predict_tbl(splmod) %>%
#   ggplot(aes(x = yday)) +
#   geom_point(aes(y = rgnd, color = w))+
#   geom_line(aes(y = .fitted))
#   scale_color_distiller(palette = "PuBu", direction = 1)



modis_spline <- modts2 %>% 
  group_by(region, year) %>% 
  mutate(df = 0.35*n()) %>% 
  nest() %>% 
  mutate(spl = map(data, ~smooth.spline(.x$yday, .x$rgnd, .x$count, df =.x$df)), # 
         splpreds  = map(spl, predict_tbl) ) %>% 
  unnest(cols = c(splpreds)) %>% 
  ungroup() %>% 
  select(region, year, yday, .fitted, rgnd, w) %>% 
  mutate(date_dummy = as.Date(yday,"1900-01-01")) %>% 
  rename(weight = w)

modis_spline %>% 
  ggplot(aes(x = date_dummy, y = rgnd, color = weight)) +
  geom_point()+
  facet_grid(cols = vars(region), rows = vars(year)) +
  geom_line(aes(y = .fitted), alpha = 0.5, color = "red")+
  scale_color_distiller(palette = "PuBu", direction = 1) +
  labs(x = "Date", y = "Snow Redness Index (SRI)")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(here("figs/modis/modis_spline.png"), height = 11, width = 8.5, units="in")


### AUC trend ---------------
auc <- modis_spline %>% 
  group_by(region, year) %>% 
  summarise(auc = sum(.fitted), max = max(.fitted)) %>% 
  ungroup()

sen <- function(..., weights = NULL) {
  mblm::mblm(...)
}

auc %>% 
  ggplot(aes(x = year, y = auc)) + 
  geom_point() +
  facet_wrap(vars(region)) +
  geom_smooth(method = sen)

auc %>% 
  group_by(region) %>% 
  nest() %>% 
  mutate(sens = map(data, ~sens.slope(.x$auc)), 
         sens_res= map(sens, tidy),
         sens_slope = map(sens, ~.x$estimates)) %>% 
  unnest(cols = c(sens_res, sens_slope)) %>% 
  select(region, sens_slope, p.value) %>% 
  arrange(p.value)

### max trend---------------------------
auc %>% 
  ggplot(aes(x = year, y = max)) + 
  geom_point() +
  facet_wrap(vars(region)) +
  geom_smooth(method = sen)

auc %>% 
  group_by(region) %>% 
  nest() %>% 
  mutate(sens = map(data, ~sens.slope(.x$max)), 
         sens_res= map(sens, tidy),
         sens_slope = map(sens, ~.x$estimates)) %>% 
  unnest(cols = c(sens_res, sens_slope)) %>% 
  select(region, sens_slope, p.value) %>% 
  arrange(p.value)

# spline S2 --------------
s2ts_splines <- s2ts2 %>%
  drop_na(mean) %>% 
  group_by(region, year) %>% 
  mutate(df = 0.35*n()) %>% 
  nest() %>% 
  mutate(spl = map(data, ~smooth.spline(.x$yday, .x$mean, .x$count, df =.x$df)), # 
         splpreds  = map(spl, predict_tbl) ) %>% 
  unnest(cols = c(splpreds)) %>% 
  ungroup() %>% 
  select(region, year, yday, .fitted, rgnd, w) %>% 
  mutate(date_dummy = as.Date(yday,"1900-01-01")) %>% 
  rename(weight = w)

s2ts_splines %>%  
  ggplot(aes(x = date_dummy, y = rgnd, color = weight)) +
  geom_point()+
  facet_grid(cols = vars(region), rows = vars(year)) +
  geom_line(aes(y = .fitted), alpha = 0.5, color = "red")+
  scale_color_distiller(palette = "PuBu", direction = 1) +
  labs(x = "Date", y = "Snow Redness Index (SRI)")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  

# compare timeseries modis/s2   ---------------------
spline_compare <- modis_spline %>% 
  rename(modis_sri=.fitted) %>% 
  select(region, year, yday, date_dummy, modis_sri) %>% 
  inner_join(rename(s2ts_splines, s2_sri = .fitted) %>% select(-weight, -rgnd) ) %>% 
  pivot_longer(cols = c(modis_sri, s2_sri)) %>% 
  mutate(name = str_remove(name, "_sri"))

spline_compare %>% 
  ggplot(aes(x = date_dummy, y = value, color = name)) +
  geom_line() +
  facet_grid(cols = vars(region), rows = vars(year))





# # scratch ---------------------
# glacier_area_raw <- dir_ls(here("data/s2_classifier_map/area/"), regexp = "T.kml$") %>%
#   map_df(st_read)
# # scale s2 so that is is within the algae mask (so its comparable with modis)
# algal_cover <- glacier_area_raw %>%
#   distinct(region, .keep_all=TRUE) %>% 
#   as_tibble() %>% 
#   mutate(glacier_km2 = as.numeric(glacier), algae_km2 = as.numeric(algae)) %>%
#   select(region, glacier_km2, algae_km2) %>% 
#   filter(region %in% c("alaskaRange", "coastNorth", "coastSouth", "interiorNorth", "interiorSouth")) %>% 
#   mutate(glac_algae_ratio = glacier_km2/algae_km2) # glacier km2 includes algae area
# 
# intensity_s2_scaled <- intensity_s2 %>% 
#   left_join(algal_cover) %>% 
#   mutate(mean = mean*glac_algae_ratio) %>% 
#   select(region, year, mean)
# 
# s2ts2 <- s2ts %>% 
#   left_join(algal_cover) %>% 
#   mutate(.fitted = .fitted*glac_algae_ratio) %>%
#   select(region, year, yday, .fitted) %>% 
#   mutate(date_dummy = as.Date(yday,"1900-01-01")) 
