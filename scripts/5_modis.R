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

theme_set(theme_bw())

# units: mean RGND of the max annual composite
# (dont attempt to apply the biomass index, hasnt been calibrated properly)
algae_raw <- read_csv(here("data/modis/meanMaxRgndRegionTerraDay100.csv")) 
doy_raw <- read_csv(here("data/modis/meanMaxRgndDoyStats100.csv"))
daymet_raw <- read_csv(here("data/modis/regionalDaymetStats1000.csv"))

#  helper f to read in ENSO data
read_index <- function(fn){
  read_table(fn, skip = 1, col_names = c("year",paste0("x",1:12)), col_types = cols(.default = "d")) %>% 
    add_column("fn"=fn)
}

# ENSO and PDO data from https://psl.noaa.gov/enso/data.html
index <- dir_ls(here("data/modis"), regexp = ".txt") %>% 
  map_df(read_index) 

country <- tribble(
  ~region, ~country,
  "alaskaRange","Alaska",
  "coastNorth","Alaska",
  "coastSouth","Canada",
  "interiorNorth","Canada",
  "interiorSouth","Canada"
)


# tidy/wrangle ------------------

algae_raw %>% 
  ggplot(aes(x = year, y = mean)) +
  geom_point() +
  facet_wrap(vars(name))

algae <- algae_raw %>% 
  drop_na() %>% 
  clean_names() %>% 
  rename(region = name, rgnd = mean) %>% 
  # filter(!(region %in% c("alaskaRange","coastNorth") & year==2009)) %>% # mt redoubt explosion
  left_join(country)

doy_raw %>% 
  ggplot(aes(x = year, y = mean)) +
  geom_point() +
  facet_wrap(vars(name))

doy <- doy_raw %>% 
  clean_names() %>% 
  select(name, year, mean, std_dev, count) %>% 
  rename(region = name, doy = mean, doy_sd = std_dev) %>% 
  # filter(!(region %in% c("alaskaRange","coastNorth") & year==2009)) %>% # mt redoubt explosion
  left_join(country)

daymet <- daymet_raw %>% 
  clean_names() %>% 
  select(name, year, contains("mean")) %>% 
  rename_with(~str_remove(.x, "_mean")) %>% 
  rename(region = name) %>% 
  left_join(country)

index_tidy <- index %>% 
  filter(year>=1999) %>%
  mutate(index = fn %>% basename() %>% str_extract("pdo|oni|nino"), .keep="unused") %>% 
  pivot_longer(x1:x12, names_to = "month") %>% 
  mutate(month = parse_number(month), 
         date  = parse_date_time(paste(year, month, sep = "-"), orders = "Ym"),
         # convert to hydrological years
         h_date = date %m+% months(3),
         h_month = month(h_date),
         h_year = year(h_date),
         value = if_else( (value > -10 & value < 10), value, NA_real_),
         season = case_when(h_month %in% 1:7~"winter", # OND JFMA accumulation season
                            h_month %in% 9:11~"summer")) %>% # JJA ablation season 
  # view() # sanity check
  select(index, h_year,season, value) %>% 
  rename(year = h_year) %>% 
  drop_na()

index_stats <- index_tidy %>% 
  # ONI and NINO3.4 basically the same, 
  filter(index != "nino") %>% 
  group_by(index,season, year) %>% 
  summarise(value = mean(value)) %>% 
  ungroup() %>% 
  unite(index, index, season) %>% 
  pivot_wider(names_from = "index")

# helper functinos ---------------

# for ggplotting
sen <- function(..., weights = NULL) {
  mblm::mblm(...)
}

# like broom::tidy, but with added columns for sens.slope estimates
tidy_sens <- function(mod, nyrs=23){
  tidy(mod) %>% 
    add_column(sens.slope = mod$estimates) %>% 
    mutate(slope.x.nyrs = sens.slope*nyrs)
}


# Phenology trend -----------------



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
MannKendall(meandoy$doy) # p =0.01
northamerica_doy_sens <- sens.slope(meandoy$doy) # -0.2
nyrs <- northamerica_doy_sens$parameter
northamerica_doy_sens$estimates*nyrs #-4.6 earlier peak bloom 2000 to 2022 
trend::sens.slope(meandoy$doy) %>% 
  tidy() %>% 
  mutate(across(c(conf.low, conf.high), ~.x*23))
# (95% CI: 1 to 8 days earlier)

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
bc_intensity <- algae %>% 
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



# JJA
doy %>% 
  left_join(algae %>% select(-count)) %>% 
  left_join(daymet) %>% 
  mutate(region = to_title_case(region)) %>% 
  pivot_longer(cols = c(doy, rgnd, jja_tmax)) %>% 
  mutate(name = fct_relevel(name, c("jja_tmax", "rgnd")) %>% 
           fct_recode(SRI = "rgnd")) %>% 
  ggplot(aes(x = year, y = value)) +
  geom_point() +
  facet_grid(rows = vars(name), cols = vars(region), scales="free") +
  geom_smooth(method = sen) +
  labs(y = "", x = "Year") +
  theme_bw() +
  scale_x_continuous(breaks = seq(2000,2020,10)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(here("figs/modis/fig4_jja.pdf"), width = 6, height = 6, units = "in")
#* in Inkscape:
#* rename row variables, put on Left axis
#* replace DOY with actual dates
#* increase space between rows, 
#*  add tags (A, B ,C)



algae %>% 
  left_join(daymet) %>% 
  left_join(doy %>% select(-count)) %>% 
  mutate(region = to_title_case(region)) %>% 
  pivot_longer(cols = c(doy, rgnd, junjul_tmax)) %>% 
  mutate(name = fct_relevel(name, c("junjul_tmax", "rgnd")) %>% 
           fct_recode(SRI = "rgnd")) %>% 
  group_by(country, year, name) %>% 
  summarise(value = weighted.mean(value, count)) %>% 
  ggplot(aes(x = year, y = value)) +
  geom_point() +
  facet_grid(rows = vars(name), cols = vars(country), scales="free") +
  geom_smooth(method = sen) +
  labs(y = "", x = "Year") +
  theme_bw() +
  scale_x_continuous(breaks = seq(2000,2020,10)) #+
  # theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(here("figs/modis/fig4_bycountry.pdf"), width = 6, height = 6, units = "in")




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
  select(region, year, mayjun_tmax, mayjun_tmin, jun_tmax, jun_tmin, jul_tmax, jul_tmin, junjul_tmax, jja_tmax) %>% 
  left_join(doy %>% select(-count))

doy_preds %>% 
  pivot_longer(cols = c(-region, -country, -year,-doy,-doy_sd)) %>% 
  ggplot(aes(x =value, y = doy )) +
  geom_point()+
  facet_wrap(vars(name), scales="free") +
  geom_smooth(method = "lm")
# all predictors suggest that a warmer spring/early summer leads to an earlier date of max RGND
doy_preds %>% 
  select(-year, -region, -doy_sd) %>% 
  ggpairs()
# tmax is signif correlated with onset date. Jul tmax is the strongest correlation




# plot of doy vs jja tmax
doy_preds %>% 
  mutate(region = to_title_case(region)) %>% 
  mutate(date_dummy = as.Date(doy, "1900-01-01")) %>% 
  ggplot(aes(x =jja_tmax, y = date_dummy )) +
  geom_point()+
  facet_wrap(vars(region), scales="free_x", nrow =1) +
  geom_smooth(method = "lm")  +
  scale_y_date(date_minor_breaks = "1 week") +
  labs(y = "Date of peak bloom", x = "June-July mean daytime air temperature (°C)")

# model coefficients, with set slope
lmer_doy_out <- lmer(doy~jja_tmax+ (1|region), data = doy_preds)
lmer_doy_out %>% summary() 
# the date of peak bloom occurred 0.8 +/- 0.35 (95% CI) days earlier
# for each 1 degree increase in mean jun-jul temperature, 
# i.e. 8 days earlier per 10 degree increase in mean summer temp

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
  select(-year, -country, -std_dev) %>% 
  group_by(region) %>% 
  nest() %>% 
  mutate(res = map(data, ~lm(rgnd~., data =.x) %>% broom::tidy())) %>% 
  unnest(res) %>% 
  filter(term!="(Intercept)", p.value<0.2) %>% 
  mutate(abs = abs(estimate)) %>% 
  arrange(p.value)

lm(rgnd~jja_tmax, data = algae_preds) %>% summary()
lm(rgnd~jja_tmax, data = algae_preds) %>% BIC()
lm(rgnd~jja_tmax+jja_precip, data = algae_preds) %>% BIC()


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

lmer1 <- lme4::lmer(rgnd ~ jja_tmax + jja_precip + (1|region), data = algae_preds) 
lmer2 <- lme4::lmer(rgnd ~ jja_tmax + (1|region), data = algae_preds)
BIC(lmer1, lmer2) # the added param of percip is not worth it
summary(lmer2)
lmer2_tidy_res <- lmer2 %>% 
  broom.mixed::tidy() %>% 
  filter(effect =="fixed") %>% 
  select(term, estimate)
intercept <- lmer2_tidy_res[1,"estimate"] %>% pull()
slope <- lmer2_tidy_res[2,"estimate"] %>% pull()
simple_lmer2 <- function(tt) slope*tt + intercept

# for an increase from 10 to 15 degrees, algal bloom intensity increased by 25%
(simple_lmer2(15)-simple_lmer2(10))/simple_lmer2(10)

# for each degree increase above 10 degrees, algal bloom intensity increased by 5%
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
  select(region, year, rgnd, jja_tmax) %>% 
  left_join(doy %>% select(region, year, doy)) %>% 
  mutate(region = to_title_case(region)) %>% 
  pivot_longer(cols = c(rgnd, doy )) %>% 
  mutate(name = fct_relevel(name, c("rgnd")) %>% 
           fct_recode(SRI = "rgnd",
                      day_of_year = "doy")) %>% 
  ggplot(aes(x = jja_tmax, y=  value)) +
  labs(y = "", x = "JJA mean daily max. air temperature (°C)") 

pred_plot + 
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(rows = vars(name), cols = vars(region), scales="free_y") +
  lims(x = c(5,20))
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




