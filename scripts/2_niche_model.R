library(tidyverse)
library(here)
library(sf)
library(broom)
library(GGally)
library(rnaturalearth)
library(glmmfields)
library(ggpubr)
library(tidymodels)
library(vip)
library(MuMIn)

theme_set(theme_bw())


# read and tidy data ------------------
world <- ne_countries(scale = "medium", returnclass = "sf")
# note that "jja" denotes jan-feb in the southern hemisphere

cov_raw <- st_read(here("data/s2_classifier_map/covariates/covariates50.kml")) 

cov <- cov_raw %>% 
  janitor::clean_names() %>% 
  select(-description:-icon) %>% 
  mutate(across(c(-name,-subregion_of, -geometry), as.numeric)) %>% 
  mutate(algae = algae_km2/glacier_km2, .keep = "unused") %>% 
  rename(region = name) %>% 
  st_centroid() %>% 
  mutate(lon= st_coordinates(geometry)[,1],
         lat= st_coordinates(geometry)[,2]) %>% 
  as_tibble() %>% 
  select(-geometry, -system_index) %>% 
  drop_na() 

colnames(cov)
nrow(cov)
  # mutate(across( contains("tmean"), ~.x-273.15  )) %>%  # convert to celsius
  # mutate(across( contains("temp"), ~.x-273.15  )) %>% 
  # mutate(across( contains("tmin"), ~.x-273.15  )) 



coast_covs <- here("data/s2_classifier_map/covariates/coastCovsPlusDist25.kml") %>% 
  st_read() %>% 
  janitor::clean_names() %>% 
  as_tibble() %>% 
  select(-description:-icon, -name) %>% 
  mutate(across(c(-geometry), as.numeric)) %>% 
  select(-system_index) %>% 
  mutate(algae = algae_km2/glacier_km2, .keep = "unused") %>% 
  mutate(across( contains("tmean"), ~.x-273.15  )) %>%  # convert to celsius
  mutate(across( contains("temp"), ~.x-273.15  )) %>% 
  mutate(across( contains("tmin"), ~.x-273.15  ))



# plot data ---------------------


cov %>% 
  as_tibble() %>% 
  select(where(is.numeric), -lon, -lat) %>% 
  pivot_longer(everything()) %>% 
  ggplot(aes(value)) +
  geom_histogram()+
  facet_wrap(vars(name), scales="free")



pd <- cov %>% 
  as_tibble() %>% 
  mutate(algae = if_else(algae>0.01, "algae", "no algae")) %>% 
  mutate(across( contains("temp"), ~.x-273.15  )) %>%  # convert to celsius
  mutate(across(contains("pres"), ~.x/1000)) %>% #convert to kPa
  # select(-contains("tmean")) %>% 
  rename_with( ~str_replace(.x, "jja", "summ") ) %>%
  select(algae, where(is.numeric), -lon, -lat) %>% 
  pivot_longer(-algae ) %>% 
  mutate(name=fct_relevel(name, "summ_prcp", after=1) %>% fct_relevel(rev)) %>% #pull(name) %>% levels()
  ggplot(aes(value, color = algae)) +
  geom_density()+
  facet_wrap(vars(name), scales="free", nrow=3,ncol =5) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(size=9) ) + 
  labs(tag="D")

pd

my_hist <- function(df){
  df %>% 
    as_tibble() %>% 
    select(where(is.numeric)) %>% 
    pivot_longer(everything()) %>% 
    ggplot(aes(value)) +
    geom_histogram() +
    facet_wrap(vars(name), scales = "free")
}
my_hist(cov)

mylog <- function(v) log(v+0.000001)

my_pp <- function(df, color=FALSE){
  df <- df %>% 
    mutate(subregion_of = if_else(subregion_of %in% c("newZealand","caucasus", "kamchatka"),"other", subregion_of)) %>% 
    mutate(log_jja_prcp = mylog(jja_prcp), .keep="unused") %>% 
    mutate(log_wint_prcp = mylog(wint_prcp), .keep="unused") %>% 
    rename(log_summ_prcp =log_jja_prcp,
           summ_pres =jja_pres,
           summ_temp =jja_temp) %>% # summer indicates jan feb in southern hemi
    select(subregion_of, where(is.numeric)) %>% 
    pivot_longer(cols = c(-subregion_of, -algae)) 
  if(color){
    df %>% 
      ggplot(aes(x = value, y = algae, color = subregion_of)) +
      geom_point(alpha = 0.3) +
      facet_wrap(vars(name), scales = "free", ncol = 3) +
      geom_smooth(method ="lm", data = df %>% filter(subregion_of !="other")) +
      scale_y_log10()+
      scale_color_brewer(palette = "Accent") +
      theme(legend.position = "bottom") +
      labs(y = "log(algal cover)", x = "")
  }else{
    df %>%   
      ggplot(aes(x = value, y = algae)) +
      geom_point(alpha = 0.3) +
      facet_wrap(vars(name), scales = "free", ncol = 3) +
      geom_smooth(method ="lm") +
      scale_y_log10() +
      labs(y = "log(algal cover)", x = "")
  }
}
my_pp(cov)
# log trans prcp vars
my_pp(cov, color = TRUE)
ggsave(here("figs/distribution/spatial_niche_algae_v_all.png"), width = 8, height = 10)
# jja pres: higher algae interior, but less stable weather, could explain the negative cor between algae and pressure in north america
# temp, obvious
# 
# #  the strong performance of winter temp is likely because in the arctic its really cold, and theres also no algae there...




## supp fig: algae vs summer temp----------------- 


## where do the various continental regions fall in this plot?
cov %>% 
  ggplot(aes(x = jja_temp, y = algae)) +
  geom_point(alpha = 0.3) +
  facet_wrap(vars(region)) 
ggsave(here("figs/distribution/covariates/algae_v_temp_region_facets.png"))


cov %>% 
  ggplot(aes(x = jja_temp, y = algae)) +
  geom_point(alpha = 0.3) +
  facet_wrap(vars(region)) +
  geom_smooth(span=1)

cov_filt <- cov %>% 
  filter(!subregion_of %in% c("caucasus", "newZealand")) %>% 
  mutate(subregion_of = fct_relevel(subregion_of, c("northAmerica", "greenland", "arctic", "europe","highMtnAsia","kamchatka","andes", "antarctica")),
         algae = algae*100)

cov_filt %>% 
  ggplot(aes(x = jja_temp, y = algae)) +
  geom_point(alpha = 0.3) +
  facet_wrap(vars(subregion_of)) +
  # geom_smooth(span=1) +
  labs(y = "Percent cover", x = "JJA temperature")


# precipitation an important factor in patagonia? 
# units: mean monthly precip in metres

# rationale for region random effect: data is suspect
cov %>% 
  # # remove regions with few positive observations
  # mutate(present_count = if_else(algae>0, 1,0)) %>% 
  # group_by(region) %>% 
  # mutate(n_pres= sum(present_count), max_algae = max(algae)) %>% 
  # filter(n_pres > 10, max_algae>0.01) %>% 
  # ungroup() %>% 
  mutate(region = fct_relevel(region, c("ellesmere", "northwestGreenland", "northeastGreenland", "westCentralGreenland", "eastCentralGreenland","baffin", "southwestGreenland", "southeastGreenland","svalbard","russianArctic",
    "alaskaRange", "coastNorth", "coastSouth"))) %>% 
  filter(algae>0.001) %>% 
  ggplot(aes(x = jja_prcp, y = algae)) +
  geom_point(alpha = 0.3) +
  facet_wrap(vars(region), scales = "free") +
  # scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method="lm")
ggsave(here("figs/distribution/4_covariates/precip.png"))
# where present
# in the arctic, except baffin and svalbard and northeastGreenland percent cover increased with summer precipitation 
# outside of the arctic, except for southwest greenland and interiorSouth percent cover declines with increasing summer precip

ggplot() +
  geom_sf(data = cov1 %>% st_as_sf, aes(fill = jja_prcp, color=jja_prcp)) +
  scale_color_viridis_c() +
  scale_fill_viridis_c()

cov %>% 
  mutate(algae_present = as.factor(if_else(algae>0.01, 1,0)) %>% fct_relevel("1")) %>% 
  ggplot(aes(x =region, y = jja_prcp, color = algae_present))+
  geom_boxplot() +
  coord_flip()
  
  

# pressure?
cov %>% 
  ggplot(aes(x = jja_pres, y = algae)) +
  geom_point(alpha = 0.3) +
  facet_wrap(vars(region)) +
  scale_x_log10() +
  scale_y_log10()
# coast range, ak peninsula have the opposite trend as most other regions

# temperature and precipitation
cov_filt %>% 
  mutate(jja_prcp= jja_prcp*100) %>% 
  filter(region!="kamchatka") %>% 
  arrange(algae) %>% 
  ggplot(aes(x = jja_prcp, y = jja_temp, color = algae)) +
  geom_point(alpha = 0.5) +
  scale_color_distiller(palette = "YlOrRd", direction=1) +
  theme_gray() +
  labs(color = "% cover", y = "Mean summer daily max. temperature (°C)", x = "Mean summer monthly precip (cm)") +
  facet_wrap(vars(subregion_of))
ggsave(here("figs/distribution/4_covariates/algae_tmp_pcp.png"))




# pressure
cov_filt %>% 
  filter(region!="kamchatka") %>% 
  arrange(algae) %>% 
  ggplot(aes(x = jja_pres, y = jja_temp, color = algae)) +
  geom_point(alpha = 0.5) +
  scale_color_distiller(palette = "YlOrRd", direction=1) +
  theme_gray() +
  labs(color = "% cover") +
  facet_wrap(vars(subregion_of))




# models -----------------------

# OLS regression
# quadratic LM with orthogonal polynomials
lm1 <- lm(algae~poly(jja_temp, 2) + poly(jja_prcp,2)+ poly(jja_pres,2), data = cov %>% drop_na(jja_temp, jja_prcp, jja_pres))
summary(lm1)
# if we set raw=TRUE the first order precip becomes "significant"



my_preprocess <- function(df, logtrans=TRUE){
  sel <- select(df, where(is.numeric))
  if(logtrans){
    mutate(sel, 
           algae = mylog(algae),
           jja_prcp = mylog(jja_prcp),
           wint_prcp = mylog(wint_prcp))
  }else{sel}
}
my_resid_plot <- function(mod) ggplot(augment(mod), aes(x = .fitted, y = .resid)) + geom_point()

# in order of least to most variables
mod0 <- lm(algae ~ jja_temp, data = my_preprocess(cov, logtrans=FALSE))
my_resid_plot(mod0)
mod1 <- lm(algae ~ jja_temp, data = my_preprocess(cov))
my_resid_plot(mod1)
mod2 <- lm(algae ~ ., data = my_preprocess(cov) )
my_resid_plot(mod2)
BIC(mod1, mod2)

# interp of the 10 var regression:

my_tidy <- function(mm) {
  tidy(mm) %>% 
    filter(term!="(Intercept)") %>% 
    arrange(p.value)
}
my_tidy(mod0) # 0.2 % increase in algal cover per degree increase
my_tidy(mod1) # univariate regression
# https://data.library.virginia.edu/interpreting-log-transformations-in-a-linear-model/
exp(0.17)- 1 # 18% increase in algal cover for every 1 degree increase in jja temp
my_tidy(mod2) # multivariate regression

# the jja max coef increases as I add in variables, although not by much with the final model
# why is jja_mean a negative coef while max and min are positive?? 




# bio01	Annual mean temperature	-290	320	°C	0.1
# bio02	Mean diurnal range (mean of monthly (max temp - min temp))	9	214	°C	0.1
# bio03	Isothermality (bio02/bio07)	7	96	%	0
# bio04	Temperature seasonality (Standard deviation * 100)	62	22721	°C	0.01
# bio05	Max temperature of warmest month	-96	490	°C	0.1
# bio06	Min temperature of coldest month	-573	258	°C	0.1
# bio07	Temperature annual range (bio05-bio06)	53	725	°C	0.1
# bio08	Mean temperature of wettest quarter	-285	378	°C	0.1
# bio09	Mean temperature of driest quarter	-521	366	°C	0.1
# bio10	Mean temperature of warmest quarter	-143	383	°C	0.1
# bio11	Mean temperature of coldest quarter	-521	289	°C	0.1
# bio12	Annual precipitation	0	11401	mm	0
# bio13	Precipitation of wettest month	0	2949	mm	0
# bio14	Precipitation of driest month	0	752	mm	0
# bio15	Precipitation seasonality	0	265	Coefficient of Variation	0
# bio16	Precipitation of wettest quarter	0	8019	mm	0
# bio17	Precipitation of driest quarter	0	2495	mm	0
# bio18	Precipitation of warmest quarter	0	6090	mm	0
# bio19	Precipitation of coldest quarter
# 
# rename(bio01_ann_mean_temp = bio01,
#        bio02_mean_diurnal_range = bio02,
#        bio03_isothermality = bio03,
#        bio04_temp_seasonality = bio04,
#        bio05_max_temp_warmest_mo = bio05,
#        bio06_min_temp_coldest_mo = bio06,
#        bio07_temp_ann_range = bio07,
#        bio08_mean_temp_wettest_q = bio08,
#        bio09_mean_temp_driest_q = bio09,
#        bio10_mean_temp_warmest_q = bio10,
#        bio11_mean_temp_coldest_q = bio11,
#        bio12_ann_precip = bio12,
#        bio13_precip_wettest_mo = bio13,
#        bio14_precip_driest_mo = bio14,
#        bio15_precip_seasonality = bio15,
#        bio16_precip_wettest_q = bio16,
#        bio17_precip_driest_q = bio17,
#        bio18_precip_warmest_q = bio18,
#        bio19_precip_coldest_q = bio19,
#        region = name)


# interior vs maritime ------------------
# 
# in case i want to map things later in native proj
# aea_string <- function(center_on, parallels){
#   lon_0 <- center_on[1]
#   lat_0 <- center_on[2]
#   lat_1 <- parallels[1]
#   lat_2 <- parallels[2]
#   paste0("+proj=aea +lat_0=",lat_0," +lon_0=",lon_0," +lat_1=",lat_1," +lat_2=",lat_2,
#          " +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +type=crs")
# }
# 
# america_proj <- aea_string(c(-135, 57),c(49, 63))

cc <- coast_covs %>% 
  st_as_sf() %>% 
  # st_transform(america_proj) %>% 
  st_centroid() %>% 
  mutate(lat = st_coordinates(geometry)[,"Y"]) 

cc %>% 
  ggplot(aes(x = distance, y = algae*100,)) +
  geom_point(alpha = 0.5)+
  labs(y ="Algal cover %", x = "Distance to Pacific Ocean (km)") 
# algal cover increases with distance from coast
ggsave(here("figs/distribution/alga_v_dist.png"))

cc %>% 
  ggplot(aes(x = distance, y = jja_temp, color = lat)) +
  geom_point() +
  scale_color_viridis_c()

cc %>% 
  ggplot(aes(x = distance, color = jja_temp, y = lat)) +
  geom_point() +
  scale_color_viridis_c()


cc %>% 
  ggplot(aes(y = algae*100, x = jja_temp)) +
  geom_point() +
  labs(y ="Algal cover %", x = "June-July temp")
# plenty of warm places have no algae
# the long tail is just a few locations in alaska





# scratch --------------------------


# logistic regression
ppa <- cov %>%
  mutate(algae = if_else(algae>0.001, 1,0)) %>%
  ggplot(aes(x = jja_temp, y = algae)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method ="glm", method.args = list(family = "binomial")) +
  labs(y = "P(algae present)", x = "", tag = "A") +
  theme_bw()
ppa

# linear regression
ppb <- cov %>%
  # filter(!region %in% c("antarctica","newZealand","caucasus", "kamchatka")) %>%
  ggplot(aes(x = jja_temp, y = algae)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method ="lm") +
  scale_y_log10() +
  labs(y = "log(Algal cover)", x = "", tag = "B") +
  theme_bw()
ppb

# linear regression, colored by region
cov_refactored <- cov %>%
  mutate(subregion_of = if_else(subregion_of %in% c("newZealand","caucasus", "kamchatka"),"other", subregion_of))
ppc <- cov_refactored %>%
  ggplot(aes(x = jja_temp, y = algae, color = subregion_of)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method ="lm", data = cov_refactored %>% filter(subregion_of !="other")) +
  scale_y_log10() +
  scale_color_brewer(palette = "Accent") + # for filtered dataset
  # scale_color_brewer(palette = "Set3") + # for full dataset
  labs(y = "log(Algal cover)", x = "Mean summer daytime air temp (°C)", tag = "C", color = "Region") +
  theme_bw() +
  theme(legend.position = "bottom")
ppc
ggpubr::ggarrange(ppa, ppb, ppc,nrow =3, heights = c(0.7, 0.7, 1))
ggsave(here("figs/distribution/spatial_niche_algae_v_temp.png"), width = 5, height = 10)




# north america only 
nacov1 <- cov_select_era %>% 
  filter(subregion_of=="northAmerica") %>% 
  select(where(is.numeric))


my_pp(nacov2)







nacov %>% 
  select(where(is.numeric)) %>% 
  pivot_longer(everything()) %>% 
  ggplot(aes(value)) +
  geom_histogram() +
  facet_wrap(vars(name), scales = "free")


na_mod1 <- lm(algae ~ ., data = nacov)
ggplot(augment(na_mod1), aes(x = .fitted, y = .resid)) + geom_point()
na_mod2 <- lm(algae ~ ., data = mutate(nacov, across(where(is.numeric), ~log(.x+0.0001))))
ggplot(augment(na_mod2), aes(x = .fitted, y = .resid)) + geom_point()

na_coefs <- tidy(na_mod1) %>% 
  filter(term!="(Intercept)") %>% 
  arrange(p.value)

na_sig_predictors <- na_coefs %>% 
  filter(p.value<0.05) %>% 
  pull(term)

nacov %>% 
  pivot_longer(cols = c(-algae)) %>% 
  filter(name %in% na_sig_predictors) %>% 
  ggplot(aes(x = algae, y =value)) +
  geom_point(alpha = 0.1) +
  facet_wrap(vars(name), scales = "free") +
  geom_smooth(method = "lm") +
  scale_x_log10()





# gamma GLM --------------------------
range(cov$algae) # cannot have zeros in gamma model

d <- cov %>% 
  as_tibble() %>% 
  select(-region, -subregion_of) %>% 
  # assume that there is no difference between a zero, and .1 %
  # reasonable given detection limits of our model, and false positives in some regions
  mutate(algae = algae + 0.001,
         slope = slope+0.001) %>% 
  drop_na() #%>% 
  # select(-contains("mean"))

# confirm zeros removed
d %>% 
  summarise(across(everything(), min))
d %>% 
  filter(if_any(everything(), ~.x==0.001))
d %>% 
  filter(if_any(everything(), ~.x==0)) #  no zeros

d %>% 
  select(where(is.numeric), -lon, -lat) %>% 
  pivot_longer(everything()) %>% 
  ggplot(aes(value)) +
    geom_histogram()+
    facet_wrap(vars(name), scales="free")


# # specify the regression type and parameters
# regression_spec <- linear_reg(penalty = 0.01) %>%
#   set_engine("glmnet", family = Gamma(link = "log")) # penalized linear regression
# 
# # specify the preprocessing recipe
# recipe_spec <- recipe(algae~., data = d) %>%
#   update_role(lat,lon, new_role = "ID") %>%
#   step_zv(all_predictors()) %>%
#   step_normalize() %>%
#   step_poly()
# prep(recipe_spec)# no polynomials added
# 
# wf <- workflow() %>%
#   add_recipe(recipe_spec) %>%
#   add_model(regression_spec)
# 
# glm_fit <- fit(wf, data = d) 
# 
# d_aug <- d %>%
#   bind_cols(predict(glm_fit$fit, d)) %>%
#   mutate(.resid = algae - .pred,
#          .rand = runif(nrow(d)))
# 
# broom:::tidy.glmnet(glm_fit$fit) %>% 
#   filter(lambda>0.01) %>% 
#   filter(lambda == min(lambda))
# # # A tibble: 10 × 5
# # term         step   estimate lambda dev.ratio
# # <chr>       <dbl>      <dbl>  <dbl>     <dbl>
# # 1 (Intercept)    45  -33.9      0.0103     0.482
# # 2 jja_temp        45   0.000690 0.0103     0.482
# # 3 wint_tmean     45   0.574    0.0103     0.482
# # 4 jja_pres        45   0.000113 0.0103     0.482
# # 5 jja_prcp        45  -3.32     0.0103     0.482
# # 6 slope          45   0.0445   0.0103     0.482
# # 7 wint_prcp      45  -2.77     0.0103     0.482
# # 8 jja_tmean       45  -0.0150   0.0103     0.482
# # 9 jja_tmin        45  -0.0751   0.0103     0.482
# # 10 wint_tmin      45  -0.424    0.0103     0.482
# 
# # nonsense...



my_formula_raw <- algae ~ jja_temp + I(jja_temp^2) + # summer temp
  wint_pres +
  jja_pres +
  slope 

m_gamma_raw <- glm(my_formula_raw, data = d, family = Gamma(link = "log"))

dredge(m_gamma_raw, rank = "AIC")



my_formula_poly <- algae ~ poly(jja_temp, 2) + # summer temp
  wint_temp +
  wint_prcp + jja_prcp +
  #wint_pres + jja_pres +
  slope 
m_gamma <- glm(my_formula_poly, data = d, family = Gamma(link = "log"))
summary(m_gamma)

# removing the worst predictor variables does not improve AIC
m_gamma2 <- glm(algae ~ . + I(jja_temp^2), 
                data = d %>% select(-lon, -lat, -wint_prcp, -wint_temp,-slope), 
                family = Gamma(link = "log"))
# m_gamma2$aic
# m_gamma3 <- glm(algae ~ ., data = d %>% select(-lon, -lat, -jja_prcp, -jja_tmin, -jja_temp), family = Gamma(link = "log"))
# m_gamma3$aic

var_importance <- vi(m_gamma2)

var_importance %>% 
  mutate(Variable = fct_reorder(Variable, Importance)) %>% 
  ggplot(aes(x = Importance, y=Variable, fill = Sign)) +
  geom_col()

d_aug <- d %>% 
  mutate(.pred = predict(m_gamma2, type="response"),
         .resid = residuals(m_gamma2, type="response"))



## plot the model predictions ---------------------------
bins <- c("0", "1–10","10–20","20–30","30–40", ">40")

# loosly based on 9 class brewer YlOrRd
# , "#ffffcc",
my_pal <- c("#FFFFFF", "#fed976", "#fd8d3c", "#fc4e2a", "#bd0026", "#800026") #"#810f7c" #  "#bd0085"

bin_percent <- function(var){
  case_when(var < 0.01 ~ bins[1],
            (var >= 0.01 & var < 0.1) ~ bins[2],
            (var >= 0.1 & var < 0.2) ~ bins[3],
            (var >= 0.2 & var < 0.3) ~ bins[4],
            (var >= 0.3 & var < 0.4) ~ bins[5],
            var >= 0.4 ~ bins[6] ) 
}

wag_proj_string <-"+proj=wag4 +lon_0=0 +datum=WGS84 +units=m +no_defs" # Wagner IV proj

d_sf <- d_aug %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% 
  st_transform(wag_proj_string) %>% 
  mutate(algae_discrete = bin_percent(algae) %>% 
           fct_relevel(bins),
         .pred_discrete = bin_percent(.pred) %>% 
           fct_relevel(bins),
         .resid_constrained = case_when(.resid< -0.5 ~ -0.5,
                                     .default = .resid),
         .resid_order = abs(.resid)
         )

world <-st_transform(world, wag_proj_string)

plot_glm <- function(color, order){
  ggplot()+
    geom_sf(data = world, fill="gray20", color=NA) +
    geom_sf(data = d_sf %>% arrange({{order}}), aes(color = {{color}}), size = 0.5) +
    theme(panel.background = element_rect(fill = "grey50"),
          panel.grid.major = element_line(colour="gray30", size=0.1)) 
}

# actual
pa <- plot_glm(algae_discrete, algae) +
  labs(color = "actual", tag="A") +
  scale_color_manual(values = my_pal) 
pa
# predicted
pb <- plot_glm(.pred_discrete, .pred) +
  labs(color = "predicted", tag="B") +
  scale_color_manual(values = my_pal) 
pb
# resids
pc <- plot_glm(.resid_constrained, .resid_order) +
  scale_color_gradient2() +
  labs(tag="C", color = "residual")
pc
# resid = algae - .pred
# blue is under, red is over prediction

ggarrange(pa, pb, pc,ncol =1)
ggsave(here("figs/distribution/4_covariates/pred.png"), height = 12, width =6, units="in")
# m_binom <- glm(is_algae ~., data = d %>% select(-lat, -lon) %>% mutate(is_algae>))







# glmmfields ----------------------------



# 
# Chain 1: Initialization between (-2, 2) failed after 100 attempts. 
# Chain 1:  Try specifying initial values, reducing ranges of constrained values, or reparameterizing the model.
# [1] "Error in sampler$call_sampler(args_list[[i]]) : Initialization failed."
# error occurred during calling the sampler; sampling not done
m_spatial <- glmmfields(algae ~ ., data = d, 
                        family = Gamma(link = "log"),
                        lat = "lat", lon = "lon", nknots = 12, iter = 100, chains = 1,
                        prior_intercept = student_t(3, 0, 10), 
                        prior_beta = student_t(3, 0, 3),
                        prior_sigma = half_t(3, 0, 3),
                        prior_gp_theta = half_t(3, 0, 10),
                        prior_gp_sigma = half_t(3, 0, 3),
                        seed = 123 # passed to rstan::sampling()
)
