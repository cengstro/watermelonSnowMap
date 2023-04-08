#* 
#* 
#* 
#* To do : landsat-8 max redness composites, visual correlation?
#*
#*
#*



# modis
library(tidyverse)
library(fs)
library(here)
library(janitor)
library(Kendall)
library(broom)
library(mblm)
library(snakecase)
library(GGally)
library(trend)
library(lme4)
library(lmerTest)
library(zyp)
library(lubridate)
library(officer)
library(flextable)

theme_set(theme_bw())

# annual intensity
terra8_raw <- read_csv(here("data/modis/meanMaxRgndRegionTerra8.csv")) 
terra_raw <- read_csv(here("data/modis/meanMaxRgndRegionTerraDay100.csv")) 
july_raw <- read_csv(here("data/modis/julyMeanMaxRgndRegionTerraDay100.csv"))
aqua_raw <- read_csv(here("data/modis/meanMaxRgndRegionAquaDay100.csv"))
aqua8_raw <- read_csv(here("data/modis/meanMaxRgndRegionAqua8100.csv"))

# annual intensity by elevation
terra_hi_raw <- read_csv(here("data/modis/meanMaxRgndRegionTerraDay100HighElevation.csv"))
terra_lo_raw <- read_csv(here("data/modis/meanMaxRgndRegionTerraDay100LowElevation.csv"))
terra8_hi_raw <- read_csv(here("data/modis/meanMaxRgndRegionTerra8100HighElevation.csv"))
terra8_lo_raw <- read_csv(here("data/modis/meanMaxRgndRegionTerra8100LowElevation.csv"))

# daily intensity time series
ts8_raw <- read_csv(here("data/modis/timeSeries8.csv")) 
ts_raw <- read_csv(here("data/modis/timeSeries.csv")) 

# DOY data
doy_raw <- read_csv(here("data/modis/meanMaxRgndDoyStats100.csv"))
july_doy_raw <- read_csv(here("data/modis/julyMeanMaxRgndDoyStats100.csv"))
aqua_doy_raw <- read_csv(here("data/modis/meanMaxRgndDoyStatsAqua100.csv"))


# temperatures stats etc
era_raw <- read_csv(here("data/modis/regionalEraStats1000.csv"))
daymet_raw <- read_csv(here("data/modis/regionalDaymetStats1000.csv"))

country <- tribble(
  ~name, ~country,
  "alaskaRange","Alaska",
  "coastNorth","Alaska",
  "coastSouth","Canada",
  "interiorNorth","Canada",
  "interiorSouth","Canada"
)


# tidy/wrangle ------------------

average_by_country <- function(df){
  df %>% 
    # filter(!(name %in% c("alaskaRange","coastNorth") & year==2009)) %>% # mt redoubt explosion
    left_join(country, by = c("name"="region")) %>% 
    group_by(country, year) %>% 
    summarise(mean = weighted.mean(mean, count)) %>% 
    ungroup()
}

# annual intensity
terra8 <-average_by_country(terra8_raw) 
terra <- average_by_country(terra_raw)
july <- average_by_country(july_raw)
aqua <- average_by_country(aqua_raw) 
aqua8 <- average_by_country(aqua8_raw)
terra_hi <- average_by_country(terra_hi_raw)
terra_lo <- average_by_country(terra_lo_raw)
terra8_hi <- average_by_country(terra8_hi_raw)
terra8_lo <- average_by_country(terra8_lo_raw)

# daily intensity time series

ts8 <- ts8_raw %>% 
  left_join(country, by = c("name"="region")) %>% 
  group_by(country, date) %>% 
  summarise(mean = weighted.mean(mean, count), count = sum(count)) %>% 
  ungroup() %>% 
  mutate(year = year(date), yday= yday(date)) %>% 
  drop_na() #%>% 
  # filter(!(country =="Alaska" & year==2009)) # mt redoubt explosion
ts <- ts_raw %>% 
  filter(name!="alaskaPeninsula") %>% 
  left_join(country, by = c("name"="region")) %>%  
  group_by(country, date) %>% 
  summarise(mean = weighted.mean(mean, count), count = sum(count)) %>% 
  ungroup() %>% 
  mutate(year = year(date), yday= yday(date)) %>% 
  drop_na() #%>% 
  # filter(!(country =="Alaska" & year==2009)) # mt redoubt explosion

## interpolate AUC and max spline ----------------------

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


ts_splines <- ts %>% 
  group_by(country, year) %>% 
  mutate(nobs =n(), df = 0.08*nobs) %>% 
  nest() %>% 
  mutate(spl = map(data, ~smooth.spline(.x$yday, .x$mean, w= .x$count, df = .x$df))) %>% 
  tidy_splines()

ts_splines %>% 
  ggplot(aes(x = date_dummy, y = glacier_mean, color = w)) +
  geom_point(alpha = 0.5)+
  facet_grid(cols = vars(year), rows = vars(country), scales = "free") +
  geom_line(aes(y = .fitted), alpha = 0.5, color = "red") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(y = "Snow redness index", x = "Date", color = "Weight")


# compute annual stats
ts_spl <- ts_splines %>% 
  # filter(!(country=="Alaska" & year==2009)) %>% # remove 2009 from AK data-- mt redoubt explorion
  group_by(country, year) %>% 
  summarise(max = max(.fitted), auc = sum(.fitted))

ts_spl_max <- ts_spl %>% 
  select(-auc)%>% 
  rename(mean = max)
ts_spl_auc <- ts_spl %>% 
  select(-max) %>% 
  rename(mean = auc)

#
ts8_splines <- ts8 %>% 
  group_by(country, year) %>% 
  mutate(nobs =n(), df = 0.3*nobs) %>% 
  nest() %>% 
  mutate(spl = map(data, ~smooth.spline(.x$yday, .x$mean, w= .x$count, df = .x$df))) %>% 
  tidy_splines()

ts8_splines %>% 
  ggplot(aes(x = date_dummy, y = glacier_mean, color = w)) +
  geom_point(alpha = 0.5)+
  facet_grid(cols = vars(year), rows = vars(country), scales = "free") +
  geom_line(aes(y = .fitted), alpha = 0.5, color = "red") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(y = "Snow redness index", x = "Date", color = "Weight")
# compute annual stats
ts8_spl <- ts8_splines %>% 
  # filter(!(country=="Alaska" & year==2009)) %>% # remove 2009 from AK data-- mt redoubt explorion
  group_by(country, year) %>% 
  summarise(max = max(.fitted), auc = sum(.fitted))

ts8_spl_max <- ts8_spl %>% 
  select(-auc)%>% 
  rename(mean = max)
ts8_spl_auc <- ts8_spl %>% 
  select(-max)%>% 
  rename(mean = auc)





## DOY ----------------
doy <- average_by_country(doy_raw)
july_doy <- average_by_country(july_doy_raw) 
aqua_doy <- average_by_country(aqua_doy_raw)

## temperatures stats etc--------------


colnames(daymet_raw) %>% 
  str_remove_all("_mean|_count|_stdDev") %>% 
  unique()

hh <- function(df){
  df %>% 
    rename(count = 3, mean = 4) %>% # dangerous way to do this
    select(-contains("stdDev")) %>% 
    left_join(country, by =c("name"="region")) %>% 
    group_by(country, year) %>% 
    summarise(mean = weighted.mean(mean, count)) %>% 
    ungroup()
}

jja_tmax_daymet <- daymet_raw %>% 
  select(name, year, contains("jja_tmax")) %>% 
  hh()
jja_tmin_daymet <- daymet_raw %>% 
  select(name, year, contains("jja_tmin")) %>% 
  hh()
jj_tmax_daymet <- daymet_raw %>% 
  select(name, year, contains("junjul_tmax")) %>% 
  hh()
jul_tmax_daymet <- daymet_raw %>% 
  select(name, year, contains("jul_tmax"), -contains("jun")) %>% 
  hh()
jul_tmin_daymet <- daymet_raw %>% 
  select(name, year, contains("jul_tmin"), -contains("jun")) %>% 
  hh()

colnames(era_raw) %>% 
  str_remove_all("_mean|_count|_stdDev") %>% 
  unique()

jja_temp_era <- era_raw %>% 
  select(name, year, contains("jja_temp")) %>% 
  hh() %>% 
  mutate(mean = mean - 273.15)
jj_temp_era <- era_raw %>% 
  select(name, year, contains("junjul_temp")) %>% 
  hh() %>% 
  mutate(mean = mean - 273.15)
jul_temp_era <- era_raw %>% 
  select(name, year, contains("jul_temp"), -contains("jun")) %>% 
  hh() %>% 
  mutate(mean = mean - 273.15)




# Trend in Intensity ---------------

# for ggplotting
sen <- function(..., weights = NULL) {
  mblm::mblm(...)
}

# like broom::tidy, but with added columns for sens.slope estimates
tidy_sens <- function(mod, nyrs=23){
  tidyres <- broom::tidy(mod) %>% 
    add_column(sens.slope = mod$estimates) %>% 
    # mutate(slope.x.nyrs = sens.slope*nyrs,
    #        conf.lo.nyrs = conf.low*nyrs,
    #        conf.hi.nyrs = conf.high*nyrs) %>% 
    select(p.value) #sens.slope, conf.low, conf.high, 
}

relative_increase <- function(df, nyrs=23){
  df <- df %>% 
    arrange(year) %>% 
    mutate(year = year-2000)
  zypmod <- zyp.sen(mean~year, dataframe = df)
  slp <- zypmod$coefficients[2]
  int <- zypmod$coefficients[1]
  frac_incr <- (slp*nyrs + int) / int
  (frac_incr - 1) * 100 # return percent increase
}



apply_mk <- function(df){
  df %>% 
    group_by(country) %>% 
    nest() %>% 
    mutate(mk = map(data, ~trend::sens.slope(.x$mean) %>% tidy_sens() ),
           pct_incr = map(data, relative_increase)) %>% 
    unnest(mk) %>% 
    unnest(pct_incr) %>% 
    select(-data) %>% 
    ungroup()
}


# Trends -----------------

annotate <- function(df, name){
  add_column(df, data=name, .before = 1)
}

annotate_pval <- function(df){
  df %>% 
    mutate(p.value = round(p.value,digits = 3),
         sig = case_when(p.value<0.001 ~ "***",
                         p.value<0.01 ~ "**",
                         p.value<0.05 ~ "*",
                         p.value<0.1~ ".",
                         .default = ""))
}


# intensity
mk_res <- annotate(apply_mk(terra8), "Terra 8-day composites (MOD09A1)") %>% 
  bind_rows(annotate(apply_mk(terra), "Terra (MOD09GA)")) %>% 
  bind_rows(annotate(apply_mk(july), "July Terra (MOD09GA)")) %>% 
  bind_rows(annotate(apply_mk(aqua), "Aqua (MYD09GA)")) %>% 
  bind_rows(annotate(apply_mk(aqua8), "Aqua 8-day composites (MYD09A1)")) %>% 
  bind_rows(annotate(apply_mk(terra_hi), "High-elevation Terra (MOD09GA)")) %>% 
  bind_rows(annotate(apply_mk(terra_lo), "Low-elevation Terra (MOD09GA)")) %>% 
  bind_rows(annotate(apply_mk(terra8_hi), "High-elevation Terra (MOD09A1)")) %>% 
  bind_rows(annotate(apply_mk(terra8_lo), "Low-elevation Terra (MOD09A1)")) %>% 
  bind_rows(annotate(apply_mk(ts_spl_max), "Max. spline Terra (MOD09GA)")) %>%  
  bind_rows(annotate(apply_mk(ts_spl_auc), "AUC spline Terra (MOD09GA)")) %>%  
  bind_rows(annotate(apply_mk(ts8_spl_max), "Max. spline Terra 8-day (MOD09A1)")) %>%  
  bind_rows(annotate(apply_mk(ts_spl_auc), "AUC spline Terra 8-day (MOD09A1)")) %>% 
  arrange(country) %>% 
  rename(region = country) %>% 
  mutate(pct_incr = round(pct_incr, digits=0)) %>% 
  annotate_pval()

print(mk_res, n=33)



colnames(mk_res) <- c("Data source", "Region", "P value", "% increase", "Significance")
mk_res %>% 
  flextable() %>% 
  save_as_docx( path = here("figs/modis/intensity_trend.docx"))

# Trend in Temp ---------------
# like broom::tidy, but with added columns for sens.slope estimates
tidy_sens2 <- function(mod, nyrs=23){
  tidyres <- broom::tidy(mod) %>% 
    add_column(sens.slope = mod$estimates) %>% 
    mutate(slope.x.nyrs = sens.slope*nyrs,
           conf.lo.nyrs = conf.low*nyrs,
           conf.hi.nyrs = conf.high*nyrs) %>%
    select(p.value, slope.x.nyrs, conf.lo.nyrs, conf.hi.nyrs) 
}

apply_mk2 <- function(df){
  df %>% 
    group_by(country) %>% 
    nest() %>% 
    mutate(mk = map(data, ~trend::sens.slope(.x$mean) %>% tidy_sens2() )) %>% 
    unnest(mk) %>% 
    select(-data) %>% 
    ungroup()
}

 
mk_res_temp <- annotate(apply_mk2(jja_tmax_daymet),"JJA tmax (DAYMET)") %>% 
  bind_rows(annotate(apply_mk2(jja_tmin_daymet),"JJA tmin (DAYMET)")) %>% 
  bind_rows(annotate(apply_mk2(jj_tmax_daymet),"Jun-Jul tmax (DAYMET)")) %>% 
  bind_rows(annotate(apply_mk2(jja_temp_era),"JJA temp (ERA5)")) %>% 
  bind_rows(annotate(apply_mk2(jj_temp_era ),"Jun-Jul temp (ERA5)")) %>% 
  arrange(country) %>% 
  relocate(data, country, slope.x.nyrs, conf.lo.nyrs, conf.hi.nyrs) %>% 
  rename(region = country, degree_increase = slope.x.nyrs, conf_lo = conf.lo.nyrs, conf_hi= conf.hi.nyrs) %>% 
  annotate_pval() %>% 
  mutate(across(c(degree_increase, conf_lo, conf_hi), ~round(.x, digits=1)))
  

mk_res_temp 
colnames(mk_res_temp) <- c("Data source", "Region", "Δ°C", "Conf. low", "Conf. high", "P value", "Significance")
mk_res_temp %>% 
  flextable() %>% 
  save_as_docx( path = here("figs/modis/temp_trend.docx"))




#SCRATCH-----------------------------


# intensity vs temp--------------------------

intensity_data <- annotate(terra8, "Terra 8-day composites (MOD09A1)") %>% 
  bind_rows(annotate(terra, "Terra (MOD09GA)")) %>% 
  bind_rows(annotate(july, "July Terra (MOD09GA)")) %>% 
  bind_rows(annotate(aqua, "Aqua (MYD09GA)")) %>% 
  bind_rows(annotate(aqua8, "Aqua 8-day composites (MYD09A1)")) %>% 
  bind_rows(annotate(terra_hi, "High-elevation Terra (MOD09GA)")) %>% 
  bind_rows(annotate(terra_lo, "Low-elevation Terra (MOD09GA)")) %>% 
  bind_rows(annotate(terra8_hi, "High-elevation Terra (MOD09A1)")) %>% 
  bind_rows(annotate(terra8_lo, "Low-elevation Terra (MOD09A1)")) %>% 
  bind_rows(annotate(ts_spl_max, "Max. spline Terra (MOD09GA)")) %>%  
  bind_rows(annotate(ts_spl_auc, "AUC spline Terra (MOD09GA)")) %>%  
  bind_rows(annotate(ts8_spl_max, "Max. spline Terra 8-day (MOD09A1)")) %>%  
  bind_rows(annotate(ts_spl_auc, "AUC spline Terra 8-day (MOD09A1)")) %>% 
  rename(algaedata = data)

tmp_data <- annotate(jja_tmax_daymet,"JJA tmax (DAYMET)") %>% 
  bind_rows(annotate(jja_tmin_daymet,"JJA tmin (DAYMET)")) %>% 
  bind_rows(annotate(jj_tmax_daymet,"Jun-Jul tmax (DAYMET)")) %>% 
  bind_rows(annotate(jja_temp_era,"JJA temp (ERA5)")) %>% 
  bind_rows(annotate(jj_temp_era,"Jun-Jul temp (ERA5)")) %>% 
  rename(tempdata = data, temp = mean)



##  daymet -----------------------

get_lm <- function(algae_df, temp_df, algae_name, temp_name){
  temp_df %>% 
    add_column(tempdata = temp_name) %>% 
    rename(temp = mean) %>% 
    left_join(algae_df %>% add_column(algaedata = algae_name)) %>% 
    nest(.by = c(country, algaedata, tempdata)) %>% 
    mutate(datasource = paste(algaedata, tempdata, sep=" v. "), .keep="unused") %>%
    mutate(lm_res = map(data, ~lm(mean~temp, data = .x) %>% tidy())) %>% 
    unnest(lm_res) %>% 
    filter(term=="temp") %>% 
    select(-term, -data) %>% 
    annotate_pval()
}
# get_lm(terra8, jja_temp_daymet, "JJA temp DAYMET","Terra 8")
  

## era-----------------------









# interaction 
dd <- terra_hi %>% 
  add_column(elev = "high") %>% 
  bind_rows(terra_lo %>% add_column(elev = "low")) %>% 
  filter(country=="Canada") %>% 
  select(-country) %>% 
  mutate(elev = as.factor(elev))

ggplot(dd, aes(x = year, y = mean, color = elev)) +
  geom_point()+
  geom_smooth(method = "lm")




##  per region ----------------------------

ggplot(doy, aes(x = year, y = doy)) +
  geom_point()+
  facet_wrap(vars(country))+
  geom_smooth(method = "lm")

doy %>% 
  group_by(region) %>% 
  nest() %>% 
  mutate(sens = map(data, ~pull(.x, doy) %>% sens.slope() %>% broom::tidy())) %>% 
  unnest(sens) %>% 
  arrange(p.value)

## North America ----------------------------
meandoy <- doy %>% 
  group_by(year) %>% 
  summarise(doy = weighted.mean(doy, count))
meandoy %>% 
  mutate(date_dummy= as.Date(doy, "1900-01-01")) %>% 
  ggplot(aes(x = year, y = date_dummy)) +
  geom_point()+
  geom_smooth(method = "lm") +
  scale_y_date(date_minor_breaks = "1 week") +
  labs(x = "Year", y = "Date of peak bloom")
MannKendall(meandoy$doy) # p =0.02
northamerica_doy_sens <- sens.slope(meandoy$doy) # -0.2
nyrs <- northamerica_doy_sens$parameter
northamerica_doy_sens$estimates*nyrs #-4.7 earlier peak bloom 2000 to 2022 
trend::sens.slope(meandoy$doy) %>% 
  tidy() %>% 
  mutate(across(c(conf.low, conf.high), ~.x*23))
# (95% CI: 1 to 9 days earlier)

## BC ----------------



doy %>% 
  left_join(country) %>% 
  group_by(country, year) %>% 
  summarise(doy = weighted.mean(doy, count)) %>% 
  nest() %>% 
  mutate(sens = map(data, ~sens.slope(.x$doy) %>% tidy_sens)) %>% 
  unnest(sens)




# Biomass trend --------------------


algae %>% 
  mutate(region =  to_title_case(as.character(region))) %>%
  ggplot(aes(x = year, y = rgnd)) +
  geom_point()+
  facet_wrap(vars(region)) +
  geom_smooth(method=sen) + # , data = t8 %>% filter(range %in% c("interior_north", "coast_south", "interior_south"))
  scale_x_continuous(breaks = seq(2000,2020,10), minor_breaks = seq(2000,2022,1)) +
  labs(y = "MODIS RGND")
algae %>% 
  group_by(region) %>% 
  nest() %>% 
  mutate(sens = map(data, ~pull(.x, rgnd) %>% sens.slope() %>% broom::tidy())) %>% 
  unnest(sens) %>% 
  arrange(p.value)



## BC ----------------------
bc_intensity <- algae2 %>% 
  filter(region %in% c("coastSouth", "interiorNorth","interiorSouth")) %>% 
  group_by(year) %>% 
  summarise(rgnd = weighted.mean(rgnd, count))

bc_intensity %>% 
  ggplot(aes(x = year, y = rgnd))+
  geom_point()+
  geom_smooth(method = sen)
bc_biomass_sens <- zyp.sen(rgnd ~ yr0, bc_intensity %>% mutate(yr0 = year-2000))
slp <- bc_biomass_sens$coefficients[2]
int <- bc_biomass_sens$coefficients[1]
(slp*nyrs + int) / int # a 77% increase over the 23 years
MannKendall(bc_intensity$rgnd) #p=0.001




# Fig 4: Tmp, SRI, DOY vs year -------------------------


doy %>% 
  left_join(algae %>% select(-count)) %>% 
  left_join(daymet) %>% 
  mutate(region = to_title_case(region)) %>% 
  pivot_longer(cols = c(doy, rgnd, junjul_tmax)) %>% 
  mutate(name = fct_relevel(name, c("junjul_tmax", "rgnd")) %>% 
           fct_recode(SRI = "rgnd")) %>% 
  ggplot(aes(x = year, y = value)) +
  geom_point() +
  facet_grid(rows = vars(name), cols = vars(region), scales="free") +
  geom_smooth(method = sen) +
  labs(y = "", x = "Year") +
  theme_bw() +
  scale_x_continuous(breaks = seq(2000,2020,10)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(here("figs/modis/fig4.pdf"), width = 6, height = 6, units = "in")
#* in Inkscape:
#* rename row variables, put on Left axis
#* replace DOY with actual dates
#* increase space between rows, 
#*  add tags (A, B ,C)






# biomass v phenology  --------------------------

doy %>% 
  select(-count) %>% 
  left_join(algae %>% select(-count)) %>% 
  mutate(date_dummy = as.Date(doy, "1900-01-01")) %>% 
  ggplot(aes(x = date_dummy, y = rgnd)) +
  geom_point()+
  geom_smooth(method = "lm", se=FALSE)+
  facet_wrap(vars(region), nrow = 1) +
  labs(x = "Date of maximum annual SRI", y = "Max. annual SRI")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
ggsave(here("figs/modis/rgnd_v_doy.png"))
# earlier onset, bigger bloom?

# alternative interp, it's just "noise" below say RGND = 0.015, eg dust and like, which tends to occur later in the season. 



# phenology predictors --------------------


## spring air temps ---------------------
# which is the best predictor?
doy_preds <- daymet %>% 
  select(region, year, mayjun_tmax, mayjun_tmin, jun_tmax, jun_tmin, jul_tmax, jul_tmin, junjul_tmax) %>% 
  left_join(doy %>% select(-count))

doy_preds %>% 
  pivot_longer(cols = c(-region, -year,-doy,-doy_sd)) %>% 
  ggplot(aes(x =value, y = doy )) +
  geom_point()+
  facet_wrap(vars(name), scales="free") +
  geom_smooth(method = "lm")
# all predictors suggest that a warmer spring/early summer leads to an earlier date of max RGND
doy_preds %>% 
  select(-year, -region, -doy_sd) %>% 
  ggpairs()
# tmax is signif correlated with onset date. Jul tmax is the strongest correlation




# plot of doy vs junjul tmax
doy_preds %>% 
  mutate(region = to_title_case(region)) %>% 
  mutate(date_dummy = as.Date(doy, "1900-01-01")) %>% 
  ggplot(aes(x =junjul_tmax, y = date_dummy )) +
  geom_point()+
  facet_wrap(vars(region), scales="free_x", nrow =1) +
  geom_smooth(method = "lm")  +
  scale_y_date(date_minor_breaks = "1 week") +
  labs(y = "Date of peak bloom", x = "June-July mean daytime air temperature (°C)")

# model coefficients, with set slope
lmer_doy_out <- lmer(doy~junjul_tmax+ (1|region), data = doy_preds)
lmer_doy_out %>% summary() 
# the date of peak bloom occurred 0.9 +/- 0.31 (95% CI) days earlier
# for each 1 degree increase in mean jun-jul temperature, 

# try allowing slopes to vary, see if there is any decrease in AIC
# lmer2_doy_out <- lmer(doy~junjul_tmax+ (junjul_tmax|region), data = doy_preds)
# AIC(lmer_doy_out,lmer2_doy_out) # no



# # all regions combined (psuedo replication, need to take weightd averages...)
# doy_temp_plot <- doy_preds %>% 
#   mutate(date_dummy= as.Date(doy, origin = "1900-01-01")) %>% 
#   ggplot(aes(x =junjul_tmax, y = date_dummy )) +
#   geom_point()+
#   # geom_errorbar(aes(ymin =doy - doy_sd, ymax = doy + doy_sd))+
#   geom_smooth(method = "lm") +
#   labs(x = "June-July mean tmax (°C)", y = "Date of peak bloom") +
#   scale_y_date(date_minor_breaks = "1 week")
# 
# doy_temp_plot
# 




## PDO -----

doy %>% 
  left_join(index_stats) %>% 
  pivot_longer(cols  = oni_summer:pdo_winter) %>% 
  ggplot(aes(x = value, y =doy)) +
  geom_point()+
  facet_grid(rows = vars(region), cols = vars(name), scales = "free") +
  geom_smooth(method="lm")
# higher summer PDO results in earlier bloom onset ?
# positive PDO -> warm and wet
# warmer = earlier melt, 
# in alaska, it might result in more snow in winter

# # full multivariate model maybe not useful, auto and temporal correlation
# doy2 %>% 
#   left_join(index_stats) %>% 
#   group_by(region) %>% 
#   nest() %>% 
#   mutate(lm_out = map(data, ~lm(doy~oni_summer+oni_winter+pdo_summer+pdo_winter, data = .x) %>% tidy())) %>% 
#   unnest(lm_out) %>% 
#   filter(term!="(Intercept)",
#          p.value<0.2) %>% 
#   arrange(p.value) 
# lmer(doy~oni_summer+oni_winter+pdo_summer+pdo_winter + (1|region), 
#      data = doy2 %>% left_join(index_stats)) %>% summary()
# 
doy %>% 
  left_join(index_stats) %>% 
  group_by(region) %>% 
  nest() %>% 
  mutate(lm_out  = map(data, ~lm(doy~pdo_summer, data = .x) %>% tidy())) %>% 
  unnest(lm_out) %>% 
  filter(term!="(Intercept)") %>% 
  arrange(p.value)
# BC earlier peak bloom is associated with positive PDO phase (warm and wet summers)-- likely earlier snowmelt


doy2 %>% 
  select(-doy_sd) %>% 
  left_join(index_stats %>% select(year, pdo_summer)) %>% 
  left_join(daymet2 %>% select(region, year, junjul_tmax)) %>% 
  ggpairs(columns  = 3:5, aes(color = region))
# higher PDO, warmer summer, 




# biomass predictors -------------------------------

algae_preds <- daymet %>% 
  left_join(algae)

algae_preds %>%
  group_by(country, year) %>% 
  summarise(across(where(is.numeric), ~weighted.mean(.x, count))) %>% 
  ungroup() %>% 
  select(country, year, rgnd, junjul_tmax, jja_tmin, jja_precip) %>% 
  GGally::ggpairs(aes(color =country))



algae_preds %>% 
  select(-max_snow_depth, -jul_tmax, -jul_tmin, -jun_tmax, -jun_tmin,-junjul_tmax,-mayjun_tmax,-mayjun_tmin, -winter_precip) %>% 
  pivot_longer(cols = c(-region, -year,-rgnd,-std_dev)) %>% 
  ggplot(aes(x =value, y =rgnd)) +
  geom_point()+
  facet_grid(cols = vars(name), rows = vars(region), scales="free") +
  geom_smooth(method = "lm")

algae_preds %>% 
  pivot_longer(cols = c(-region, -year,-rgnd,-std_dev)) %>% 
  ggplot(aes(x =value, y =rgnd )) +
  geom_point()+
  facet_grid(rows = vars(region), cols = vars(name), scales="free") +
  geom_smooth(method = "lm")
# tmax increases RGND, less so tmin
# coast south has higher RGND in warmer winters
# precip decreases RGND
# no effect of max "snowdepth"-- however this was calculated

algae_preds %>% 
  pivot_longer(cols = c(-region, -year,-rgnd,-std_dev)) %>% 
  ggplot(aes(x =value, y =rgnd )) +
  geom_point()+
  facet_grid(rows = vars(region), cols = vars(name), scales="free") +
  geom_smooth(method = "lm")



algae_preds %>% 
  select(-year, -std_dev) %>% 
  group_by(region) %>% 
  nest() %>% 
  mutate(res = map(data, ~lm(rgnd~., data =.x) %>% broom::tidy())) %>% 
  unnest(res) %>% 
  filter(term!="(Intercept)", p.value<0.2) %>% 
  mutate(abs = abs(estimate)) %>% 
  arrange(p.value)

lm(rgnd~junjul_tmax, data = algae_preds) %>% summary()
lm(rgnd~junjul_tmax, data = algae_preds) %>% BIC()
lm(rgnd~junjul_tmax+jja_precip, data = algae_preds) %>% BIC()


# wtd_mean <- algae_preds %>%
#   select(region, year, rgnd, junjul_tmax) %>%
#   left_join(select(algae, region, year, count)) %>%
#   group_by(year) %>%
#   summarise(rgnd = weighted.mean(rgnd, count),
#             junjul_tmax = weighted.mean(junjul_tmax, count)) %>%
#   ungroup()
# ggplot(wtd_mean, aes(x = junjul_tmax, y = rgnd)) +
#   geom_point() +geom_smooth(method = "lm")
# lm(rgnd ~ junjul_tmax, data = wtd_mean) %>% summary()
# # for an increase in temperature from 10 to 15 deg,
# # algal bloom intensity increased by ~25%

# BIC(lm(rgnd~., data = wtd_mean)) # very high---- not 100% sure the max likelihood is comparable between lm and lmer

lmer1 <- lme4::lmer(rgnd ~ junjul_tmax + jja_precip + (1|region), data = algae_preds) 
lmer2 <- lme4::lmer(rgnd ~ junjul_tmax + (1|region), data = algae_preds)
BIC(lmer1, lmer2) # the added param of percip is not worth it
summary(lmer2)
lmer2_tidy_res <- lmer2 %>% 
  broom.mixed::tidy() %>% 
  filter(effect =="fixed") %>% 
  select(term, estimate)
intercept <- lmer2_tidy_res[1,"estimate"] %>% pull()
slope <- lmer2_tidy_res[2,"estimate"] %>% pull()
simple_lmer2 <- function(tt) slope*tt + intercept

# for an increase from 10 to 15 degrees, algal bloom intensity increased by 20 %
(simple_lmer2(15)-simple_lmer2(10))/simple_lmer2(10)

# for each degree increase above 10 degrees, algal bloom intensity increased by 4 %
(simple_lmer2(11)-simple_lmer2(10))/simple_lmer2(10)


lmer3 <- lmerTest::lmer(rgnd ~ junjul_tmax + (1|region), data = algae_preds)
summary(lmer3) # p=0.006

# plot the relationship between precip and temperature
algae_preds %>% 
  ggplot(aes(x = jja_precip, y = jja_tmax)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(vars(region))
# wetter summers tend to be cooler


algae_preds %>% 
  group_by(region)


## PDO -----------


algae %>% 
  left_join(index_stats) %>% 
  pivot_longer(cols  = oni_summer:pdo_winter) %>% 
  ggplot(aes(x = value, y =rgnd)) +
  geom_point()+
  facet_grid(rows = vars(region), cols = vars(name), scales = "free") +
  geom_smooth(method="lm")

algae %>% 
  left_join(index_stats) %>% 
  group_by(region) %>% 
  nest() %>% 
  mutate(lm_out  = map(data, ~lm(rgnd~pdo_summer, data = .x) %>% tidy())) %>% 
  unnest(lm_out) %>% 
  filter(term!="(Intercept)") %>% 
  arrange(p.value)
# no effect of PDO on algal bloom intensity





# combined predictor plot (supp) ---------------------


pred_plot <- algae_preds %>% 
  select(region, year, rgnd, junjul_tmax) %>% 
  left_join(doy %>% select(region, year, doy)) %>% 
  mutate(region = to_title_case(region)) %>% 
  pivot_longer(cols = c(rgnd, doy )) %>% 
  mutate(name = fct_relevel(name, c("rgnd")) %>% 
           fct_recode(SRI = "rgnd",
                      day_of_year = "doy")) %>% 
  ggplot(aes(x = junjul_tmax, y=  value)) +
  labs(y = "", x = "June-July mean daily max. air temperature (°C)") 

pred_plot + 
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(rows = vars(name), cols = vars(region), scales="free_y") 
ggsave(here("figs/modis/algae_v_temp_regions.pdf"), width = 10, height = 5, units = "in")

# display all regions in a single panel 
pred_plot + 
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(rows = vars(name), scales="free_y")
ggsave(here("figs/modis/algae_v_temp.pdf"), width = 5, height = 5, units = "in")

# all regions in a single panel, coloured by region
pred_plot + 
  geom_point(aes(color = region)) +
  geom_smooth(aes(color = region), method = "lm") +
  facet_grid(rows = vars(name), scales="free_y")
ggsave(here("figs/modis/algae_v_temp_region_colors.pdf"), width = 5, height = 5, units = "in")
# chaotic




