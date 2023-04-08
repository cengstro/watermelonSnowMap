#* glacier mass balance vs RGND
#* using published mass balance data for the 
#* Eklutna, Gulkana glaciers in Alaska
#* G211009E61198N, G214576E63274N   respectively
#* 
#* did not include Kahiltna because it has less than 10% algal cover (G208771E62849N)
#* 
#* UPSTREAM
#* 09.0_trend.js (rgnd data per glacier)
#* png snapshots from publications containing mass balance data
#*         geck_2021_eklutna, o'neel_2019_gulkana and young_2018_kahiltna
#* 
#* 
#* 
#* 
library(tidyverse)
library(here)
library(fs)
library(lubridate)
library(corrr)




# Digitize the plots --------------
# 
# library(digitize)
# 
# EKLUTNA GLACIER MASS BALANCE DATA
# dd <- ReadAndCal(here("data/mass_bal/geck_2021_eklutna.png"))
# data.points = DigitData(col = 'red')
# df = Calibrate(data.points, dd, 2000, 2015, -4, 0)
# df
# plot(df)
# eklutna <- df %>% 
#   mutate(year = x %>% round(), .keep="unused")
# 
# eklutna %>% 
#   rename(summer_mass_bal_mwe = y) %>% 
#   write_csv(here("data/mass_bal/eklutna.csv"))
# 
# 
# 
# 
# # GULKANA GLACIER MASS BALANCE DATA
# dd2 <- ReadAndCal(here("data/mass_bal/o'neel_2019_gulkana.jpeg"))
# data.points2 = DigitData(col = 'red')
# df2 = Calibrate(data.points, dd2, 2000, 2010, -5, 0)
# df2
# plot(df2)
# gulkana <- df2 %>% 
#   mutate(year = x %>% round())
# 
# gulkana %>% 
#   rename(summer_mass_bal_mwe = y) %>% 
#   write_csv(here("data/mass_bal/gulkana.csv"))
# 
# 
# 
# # KAHILTNA GLACIER MASS BALANCE DATA
# dd3 <- ReadAndCal(here("data/mass_bal/young_2018_kahiltna.jpeg"))
# data.points3 = DigitData(col = 'red')
# df3 = Calibrate(data.points3, dd3, 2004, 2013, -2.5, 0)
# df3
# plot(df3)
# kahiltna <- df3 %>% 
#   mutate(year = x %>% round(), .keep="unused")
# 
# kahiltna %>% 
#   rename(summer_mass_bal_mwe = y) %>% 
#   write_csv(here("data/mass_bal/kahiltna.csv"))
# 



# # Black rapids
# dd2 <- ReadAndCal(here("data/mass_bal/black_rapids_kienholz_2017.jpg"))
# data.points2 = DigitData(col = 'red')
# df2 = Calibrate(data.points2, dd2, 1980, 2017, -6, 0)
# df2 %>% 
#   mutate(year =round(x)) %>% 
#   filter(year>2010) %>% 
#   mutate()


# read in CSVs -----------------------


peytos <- dir_ls(here("data/mass_bal/peyto_massbal_pradhananga_etal_2021/")) %>%
  map_df(read_tsv, skip = 3, id = "filename")
  


# the data from the publications
gulkana_mb <- read_csv(here("data/mass_bal/gulkana.csv")) 
eklutna_mb <- read_csv(here("data/mass_bal/eklutna.csv")) 
# modis RGND data
rgnd <- read_csv("data/time_series/modis/meanAnnSDIPerGlacier.csv") 
daymet <- read_csv(here("data/mass_bal/gulkana_eklutna_daymet_daily.csv"))



#  the combined data --------------
smb <- read_csv(here("data/mass_bal/summer_mass_balance_data.csv"))
distinct(smb, glac_name)
distinct(smb, reference)
# tidy --------------------------

peyto2 <- peytos %>% 
  rename(datetime = 2, ice_ht_m = sr50.1) %>% 
  mutate(site = basename(filename) %>% str_sub(6,10),
         datetime = parse_date_time(datetime, "Y m d h s"),
         date = as_date(datetime),
         hour = hour(datetime), .keep = "unused")

#  summer mass balance
peyto_smb <- peyto2 %>% 
  group_by(site, date) %>% 
  summarise(ice_ht_m = mean(ice_ht_m)) %>% 
  ungroup() %>% 
  filter(month(date) %in% 5:9) %>% 
  group_by(site, year = year(date)) %>% 
  summarise(max = max(ice_ht_m), min = min(ice_ht_m)) %>% 
  ungroup() %>% 
  mutate(smb = min - max) %>% 
  drop_na()

peyto_smb %>% 
  ggplot(aes(x  = year, y = smb, color = site)) +
  geom_point()

# simply estimating mass balance as average of upper and lower, and extrapolating to entire glacier surface
peyto_smb %>% 
  filter(year %in% c(2011, 2012, 2013)) %>% 
  summarise(mean = mean(smb))








mb <- gulkana_mb %>% 
  add_column(glac_id="G214576E63274N", glac_name = "Gulkana") %>% 
  bind_rows(add_column(eklutna_mb, glac_id="G211009E61198N", glac_name ="Eklutna"))

mb %>% 
  filter(year %in% 2000:2018) %>% 
  summarise(mean = mean(summer_mass_bal_mwe), ci = sd(summer_mass_bal_mwe)*1.96)
  


rgnd2 <- rgnd  %>% 
  mutate(rgnd_mean_ja = sdi_jja_mean) %>% # correct a typo
  select(glac_id, year, rgnd_mean_ja) %>% 
  filter(glac_id %in% c("G211009E61198N", "G214576E63274N", "G208771E62849N"))

daymet2 <- daymet %>% 
  select(glac_id, date, prcp, tmax, tmin, insol) %>% 
  mutate(year = year(date)) %>% 
  filter(year>2000) %>% 
  group_by(glac_id, year) %>% 
  summarise(pdds_tmax = sum(tmax[tmax>0]), 
            pdds_tmin = sum(tmin[tmin>0]),
            prcp = sum(prcp),
            insol = sum(insol)) %>% 
  ungroup()


# join the two datasets
cc <- mb %>% 
  inner_join(rgnd2) %>% 
  inner_join(daymet2)


# plot --------------

cc %>% 
  select(summer_mass_bal_mwe, rgnd_mean_ja, pdds_tmax, pdds_tmin, prcp, insol, glac_name) %>% 
  GGally::ggpairs(aes(color = glac_name))

cc %>% 
  ggplot(aes(y = summer_mass_bal_mwe, x = rgnd_mean_ja, color = glac_id)) +
  geom_point()

# model ------------------
library(lme4)
mod1 <- lm(summer_mass_bal_mwe ~ rgnd_mean_ja + pdds_tmin + prcp + insol, data = cc)
mod2 <- lm(summer_mass_bal_mwe ~ scale(rgnd_mean_ja) + scale(pdds_tmin) + scale(prcp) + scale(insol), data = cc)
mod3 <- lm(summer_mass_bal_mwe ~ log(rgnd_mean_ja) + log(pdds_tmin) + log(prcp) + log(insol), data = cc)
mod4 <- lmer(summer_mass_bal_mwe ~ rgnd_mean_ja + pdds_tmin + prcp + insol + (1|glac_id), data = cc)
summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)

sen <- function(..., weights = NULL) {
  mblm::mblm(...)
}

cc %>% 
  filter(!year %in% c(2004,2005, 2009)) %>% 
  ggplot(aes(x = mean_rgnd, y = summer_mass_bal_mwe)) +
  geom_point() +
  facet_wrap(vars(glac_name)) +
  geom_smooth(method = sen) +
  labs(x = "mean RGND (max ann comp)") +
  ggrepel::geom_text_repel(aes(label = year))

# Gulkana is the one in the Hayes range-- smallest algae area
# Eklutna has the highest percent algae cover
# kahiltna is massive, although on 




cc %>% 
  filter(glac_name=="Eklkutna Glacier") %>% 
  select(summer_mass_bal_mwe, mean_rgnd) %>% 
  correlate()


# try with algae3 masked data-- is the effect stronger or weaker?

# try with the JA mean RGND (expect mean to be better predictor than mean of max composite)
# daily_algaemasked_rgnd_select_glaciers, and daily_rgnd_select_full_glacier




# same as above, but with mean Jul-Aug RGND ---------------------------------------
# ie downloading the daily timeseries

# V3 algae data
# Gulkana-- 5/17.5=29% coverage
# Eklutna 8.5/29.5= 29% coverage
# Kahiltna 42/485 = 9%



# the duration of the RGND should impact the bloom as well, 
# makes sense to compare mean JA with mean JA
# data computed using 09_modis_rgnd
algae_rgnd_raw <- read_csv(here("data/time_series/modis/daily_algaemasked_rgnd_select_glaciers.csv"))
glac_rgnd_raw <- read_csv(here("data/time_series/modis/daily_rgnd_select_full_glacier.csv"))

clean_rgnd <- function(df){
  df %>%   janitor::clean_names() %>% 
    select(glac_name, date, mean, count, std_dev) %>% 
    mutate(glac_name = if_else(glac_name=="None","Kahiltna Glacier",glac_name)) %>% 
    rename(mean_rgnd = mean, std_dev_rgnd = std_dev, count_rgnd = count) %>% 
    # mutate(mean_rgnd = mean_rgnd/1000) %>% # convert back to original RGND units?? I think I had gotten rid of the divide by 1000 when i re-ran this... but it was definitely slow, so might want to bring it back
    filter(count_rgnd> 10, month(date) %in% c(7,8))
}
algae_rgnd <- clean_rgnd(algae_rgnd_raw)
glac_rgnd <- clean_rgnd(glac_rgnd_raw)

mean_rgnd <- function(df){
  algae_rgnd %>% 
    group_by(glac_name, year = year(date)) %>% 
    summarise(max_rgnd = max(mean_rgnd), mean_rgnd = mean(mean_rgnd)) %>% 
    ungroup()
}

algae_mean <- mean_rgnd(algae_rgnd) %>% inner_join(mb)
glac_mean <- mean_rgnd(glac_rgnd) %>% inner_join(mb) %>% print()

library(ggrepel)
algae_mean %>% 
  filter(!year %in% c(2004,2005, 2009)) %>% 
  mutate(glac_name = case_when(glac_name=="Gulkana Glacier" ~ "Gulkana Glacier (29%)",
                               glac_name=="Eklutna Glacier" ~ "Eklutna Glacier (29%)",
                               glac_name=="Kahiltna Glacier" ~ "Kahiltna Glacier (9%)")) %>% 
  ggplot(aes(x = mean_rgnd, y = summer_mass_bal_mwe)) +
  geom_point() +
  facet_wrap(vars(glac_name), ncol = 1) +
  # geom_smooth(method = "lm") +
  labs(x = "RGND (Jul-Aug mean)", y = "Summer mass balance m.w.e.") #+
  # geom_text_repel(aes(label = year))

# 2004, 2005, 2009... hmm, wasn't 2009 the Mt Redoubt explosion??

# how do I know what's algae, and what's mineral dust?
ggsave(here("figs/glac_mb_clean.png"))

# data from geck_2021_eklutna (2000-2019), o'neel_2019_gulkana(2000-2018) and young_2018_kahiltna(2004-2013)
algae_mean %>% 
  group_by(glac_name) %>% 
  summarise(min(year), max(year))
# lm coefs
algae_mean %>% 
  group_by(glac_name) %>% 
  group_map(~ broom::tidy(lm(summer_mass_bal_mwe ~ mean_rgnd, data = .x)))

# weaker assoc w max, but I havent cleaned the data, likely lots of noisy observations, 
# and from a RF perspective it makes more sense that the MEAN would be most important (duration)
# as a bonus, we did not filter out clouds in this iteration of the modis rgnd0, so it takes into account cloudy spells
algae_mean %>% 
  ggplot(aes(x = max_rgnd, y = summer_mass_bal_mwe)) +
  geom_point() +
  facet_wrap(vars(glac_name)) +
  geom_smooth(method = "lm") +
  labs(x = "max Jul-Aug RGND")


















# peyto stream flow -------------------------
pp <- read_csv(here("data/mass_bal/Peyto_mass_balance_wgms_1966_2017.csv"))
glimpse(pp)
ppp <- pp %>% 
  janitor::clean_names() %>% 
  select(name, year, summer_balance) %>% 
  rename(glac_name = name, summer_mass_bal_mwe = summer_balance) %>% 
  filter(glac_name=="PEYTO", year>1999) %>% 
  mutate(glac_name = glac_name %>% str_to_title() %>% paste("Glacier"))
view(ppp)




# data from https://essd.copernicus.org/articles/13/2875/2021/

library(lubridate)
ss <- read_csv(here("data/mass_bal/stream/PeytoOutlet_streamflow_dly_2013_2018.csv"))
sss <- ss %>% 
  janitor::clean_names() %>% 
  rename(date = datetime)

ss1 <- ss %>% 
  janitor::clean_names() %>% 
  filter(month(datetime) %in% c(7, 8)) %>% 
  group_by(year = year(datetime)) %>% 
  summarise(mean_flow = mean(discharge_m3_s), sum_flow = sum(discharge_m3_s))
peyto_stream <- rr %>% 
  filter(glac_name=="Peyto Glacier") %>% 
  inner_join(ss1) 

peyto_stream %>% 
  ggplot(aes(x=sum, y = mean_flow)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "sum RGND (max ann comp method)", y = "mean JJA streamflow (m3/s)",
       title="Summer streamflow below the Peyto Gl, AB versus total algal biomass")
lm(mean_flow~mean, data = peyto_stream) %>% broom::tidy()


# NSD Jul, or JA data. signif JJA data. 
# try looking at daily RGND on glacier, mask cloudy days and days with rainfall, 
# is there a connection?


library(janitor)

# check peyto stream flow again, this time looking at daily rgnd data
peyto_rgnd <- read_csv(here("data/time_series/modis/peyto_daily_terra_rgnd.csv"))
peyto <- peyto_rgnd %>% 
  janitor::clean_names() %>%
  filter(rgnd0_count ==3) %>% # the area hsa 3 MODIS grid cells
  drop_na(rgnd0_mean) %>% 
  inner_join(sss) %>% print()

peyto %>% 
  ggplot(aes(x = date)) + 
  geom_line(aes(y = rgnd0_mean*200), color = "red") +
  geom_line(aes(y = discharge_m3_s), color = "blue") +
  facet_wrap(vars(year(date)), scales="free_x")

peyto %>% 
  filter(month(date)==7) %>% 
  ggplot(aes(x = date)) + 
  geom_line(aes(y = rgnd0_mean*200), color = "red") +
  geom_line(aes(y = discharge_m3_s), color = "blue") +
  facet_wrap(vars(year(date)), scales="free_x")

peyto %>% 
  filter(month(date)==7) %>% 
  group_by(year(date)) %>% 
  summarise(sum(discharge_m3_s))

# which events are due to rainfall?
rain <- read_tsv(here("data/mass_bal/stream/ERAI_tRHuQsiQlip_BiasCorrect2PeytoBow_Jan1979_Aug2019.obs.txt"), skip=8)
glimpse(rain)
rain2 <- rain %>% 
  clean_names() %>% 
  rename(date = 1) %>% 
  mutate(date = date %>% str_sub(1,10) %>% parse_date(format = "%Y %m %d")) %>% 
  filter(year(date) %in% 2013:2018, month(date) %in% c(6,7,8,9)) %>% 
  group_by(date) %>% 
  summarise(p_1 = sum(p_1)) %>% 
  ungroup()
rain2
peyto2 <- peyto %>% 
  left_join(rain2) %>% 
  select(date, rgnd0_mean, discharge_m3_s, p_1)
peyto2 %>% 
  filter(month(date) %in% 7:8) %>% 
  ggplot(aes(x = date)) + 
  geom_line(aes(y = rgnd0_mean*200), color = "red") +
  geom_line(aes(y = discharge_m3_s), color = "blue") +
  geom_col(aes(y = p_1/5), color = "green") +
  facet_wrap(vars(year(date)), scales="free_x")

# even if we find a correlation between a streamflow peak and a big algal bloom
# who is to say that it's not due to the prolonged snowcover associated with the algal bloom?
# or the high temps that co-occur with peak algae season?
# would need accurate covs to disentangle this effect... lacking this data, there's not much we can conclude







# SCRATCH -----------------



# read_csv(here("data/time_series/modis/glacier_mass_bal_max_comp_algae3.csv")) #modis_rgnd.js for Eklutna, Gulkana, and Kahiltna RGND data
