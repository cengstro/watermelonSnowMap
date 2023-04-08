# updated 2023-02-06
#*
#*
#* TO DO
#* use F1 (or pr auc?) instead of accuracy for VI metric, are all variables positive?
#*
#*
#*
#*
#*
#*
library(tidyverse)
library(here)
library(janitor)
library(fs)
library(umap)
# library(ggrepel)
library(sf)
library(tidymodels)
library(probably)
library(vip)
library(recipeselectors)
library(lubridate)
library(rnaturalearth)
library(cowplot)

# A) download most recent data from G drive
# # SLOW ~5 min download the files
# drive_files <- drive_ls(path = "training", pattern = ".csv", recursive = TRUE) %>% print()
# local_dir <- here(paste0("data/s2_classifier_map/training/combined/"))#, mapid, "/"
# local_paths <- paste0(local_dir, drive_files %>% pull(name))
# drive_files %>%
#   pull(id) %>%
#   walk2(local_paths, ~drive_download(as_id(.x), path = .y, overwrite = TRUE))#

# B) load current local data
local_paths <- dir_ls(here("data/s2_classifier_map/training/combined/"), regexp="csv")

# I defined regions differently in the training data collection stage,
# this dictionary equates the two "region" collections
region_dict <- read_csv(here("data/s2_classifier_map/training_region_name_dict.csv"))
regions <- st_read(here("data/s2_classifier_map/glacierRegionsV15.kml"))

train_raw <- local_paths %>% 
  map_df(read_csv, id = "path") 

# tidy -----------------
train_full <- train_raw %>% 
  janitor::clean_names() %>% 
  mutate(basename = basename(path) %>% str_remove(".csv"), 
         region = str_sub(basename, 1,2),
         hrs = str_split_fixed(basename, "_", 4)[,3],
         class_name = as.factor(class) %>% 
           fct_recode(algae = "1",
                      localDust = "2",
                      windDust = "3",
                      glare = "4",
                      rainbow = "5",
                      shadow = "6",
                      haze = "7",
                      ice = "8",
                      crevasse = "9",
                      other = "10",
                      whiteSnow = "11") %>% 
           as.character(), 
         coords = str_extract(geo, "(?<=\\[)[:graph:]+(?=\\])"),#[:digit:]+ (?=\\,)
         lon = str_split_fixed(coords, ",", 2)[,1] %>% parse_number(),
         lat = str_split_fixed(coords, ",", 2)[,2] %>% parse_number(), #%>% as.numeric(length=10)
         class = if_else(class==1,1,0) %>% as.factor(), # recode to 0 or 1
         .keep="unused") %>% 
  select(-system_index, -coords)
glimpse(train_full)

#  join the training region id as property of the regions featre collection
region_dict2 <- region_dict %>% 
  select(-`system:index`)
regions2 <- regions %>% 
  select(Name,subregionOf,geometry) %>% 
  rename(name = Name) %>% 
  left_join(region_dict2)


# data frequency -----------------------

# n pts
train_full %>% distinct(basename) %>% nrow() # 225 images
nrow(train_full) # total = 9312
count(train_full, class) # per class
# # A tibble: 2 × 2
# class     n
# <fct> <int>
#   1 0      8199
# 2 1      1115
count(train_full, region)  #per region
count(train_full, region, class) %>% pivot_wider(names_from = "class", values_from = "n") # per class, per region

# number of "images" (actually mosaics)
distinct(train_full, basename) %>% nrow() # total
distinct(train_full, basename, region) %>% 
  count(region) # per region

## temoporal distribution--------------------
do_date_dummy <- function(date) as_date(str_replace(date, "20[:digit:]{2}", "1900"))

dates <- train_full %>% 
  select(basename, class_name, region, hrs, lat, lon) %>% 
  mutate(date = as_datetime(as.numeric(hrs)*60*60), 
         year =year(date),
         date_dummy = do_date_dummy(date),
         hemi = if_else(lat>0, "north", "south"))

dates %>% 
  filter(month(date)!=12) %>% # some late december in S hemisphere due, actually in early jan but mosaics labelled by start day of week
  ggplot(aes(date_dummy)) +
  geom_histogram() +
  facet_wrap(vars(hemi), scales="free_x", ncol=1)

classname_colors <-
  c(algae = "#ff0000", # red
    localDust = "#b87524", # brown
    windDust = "#5a1b95", # purple
    glare = "#ce7e00", #orange
    rainbow = "#fffd0e", #yellow
    haze = "#1ace00",#green
    shadow = "#1100ce",#blue
    crevasse = "#ff00ff", # pink
    whiteSnow = "#adadad", # grey
    ice = "black", # black
    other = "pink")

dates %>% 
  filter(month(date)!=12) %>% # some late december in S hemisphere due, actually in early jan but mosaics labelled by start day of week
  ggplot(aes(x = date_dummy, fill =class_name )) +
  geom_histogram() +
  facet_grid(cols = vars(hemi), 
             rows = vars(year), scales="free_x")  +
  scale_fill_manual(values = classname_colors) 
  




## spatiotemporal distribution of training points-----


train_sf <- train_full %>%
  mutate(date = as_datetime(as.numeric(hrs)*60*60), 
         year =year(date),
         date_dummy = do_date_dummy(date),
         # potentially misleading december dates in southern hemisphere-- manually shift to jan 1
         date_dummy = if_else(month(date_dummy)==12,
                              as_date("1900-01-01"),
                              date_dummy),
         longi = lon, lati = lat) %>% # filter(is.na(date)) %>% view() # ok, missing some meta from a few in north america
  group_by(basename, date, class_name) %>%
  mutate(n = n()) %>% 
  nest() %>%
  mutate(data = map(data, slice_sample, n=1)) %>%
  unnest() %>%
  ungroup() %>%
  st_as_sf(coords = c("lon", "lat"),
           crs ="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") # default mercator

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

train_sf %>% 
  ggplot() +
  geom_sf(data = world, fill="gray20") +
  geom_sf(data = train_sf, aes(color = region))

aea_string <- function(center_on, parallels){
  lon_0 <- center_on[1] # map centered on longitude
  lat_0 <- center_on[2] # ... latitude
  lat_1 <- parallels[1] # the parallels at which area is un-distorted
  lat_2 <- parallels[2]
  paste0("+proj=aea +lat_0=",lat_0," +lon_0=",lon_0," +lat_1=",lat_1," +lat_2=",lat_2,
         " +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +type=crs")
}
america_proj <- aea_string(c(-135, 57), c(50,70)) #49, 63
gris_proj <- aea_string(c(-33, 66), c(50,70)) #60, 75
europe_proj <- aea_string(c(10,52), c(46, 61))
kamchatka_proj <- aea_string(c(159,56), c(50,60))
arctic_proj <- sf::st_crs("ESRI:102016")
hma_proj <- aea_string(c(75,36), c(30,40))
andes_proj <- aea_string(c(-68,-25), c(-20, -50))
antarctica_proj <- aea_string(c(-63,-68), c(-65,-75))
nz_proj <- aea_string(c(172,-42), c(-40,-46))


states <- rnaturalearth::ne_states(c("united states of america", "canada"), returnclass = "sf") # canada and USA states

my_aea_plot <- function(title, code, proj, xlim = c(-1e6, 1e6), ylim = c(-1e6, 1e6), border_sf = world, minlon=-180, maxlon=180, minlat = -90, maxlat = 90, orient = "h"){
  train_sf2 <- train_sf %>% 
    filter(region %in% code) %>% 
    filter(longi>minlon, longi<maxlon) %>% 
    filter(lati>minlat, lati<maxlat)
  
  p1 <- ggplot() +
    geom_sf(data = border_sf %>% st_transform(proj), fill="gray20") +
    geom_sf(data = train_sf2 %>% st_transform(proj), aes(color = class_name)) +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) + # for america_proj
    scale_color_manual(values = classname_colors) +
    ggtitle(title)+
    theme(legend.position = "none") +
    ggspatial::annotation_scale(location = 'bl', width_hint = 0.5)
  
  p2 <- train_sf2 %>%
    drop_na(year) %>% 
    ggplot(aes(x = date_dummy, y = n, fill =class_name )) +
    geom_bar(position="stack", stat="identity", width = 6) +
    scale_fill_manual(values = classname_colors) 
    
  
  if(orient == "h"){
    p2 <-  p2 + facet_grid(cols = vars(year))
    plot_grid(p1,p2 , nrow = 2)
  }else{
    p2 <-  p2 + facet_grid(rows = vars(year))
    plot_grid(p1,p2 , ncol = 2)
  }
  
}

my_aea_plot("north america",c("na", "vo", "no"),america_proj, 
            xlim = c(-1.5e6, 1.5e6), ylim = c(-1e6, 0.9e6), 
            border_sf = states, 
            minlon = -180, maxlon = 0,maxlat = 66) # remove kamchatka from "vo"
my_aea_plot("subarctic greenland", c("gr","eu"), gris_proj, 
            xlim = c(-1.1e6, 1e6), ylim = c(-1e6, 0.3e6),
            maxlon = 0, maxlat = 66) # remove rest of europe outside iceland
my_aea_plot("europe", c("eu","no"), europe_proj, 
            xlim = c(-0.5e6, 0.5e6), ylim = c(-1e6, 1.5e6),
            minlon = 0, orient = "v", maxlat=66) # minus iceland
my_aea_plot("kamchatka", "vo", kamchatka_proj, 
            xlim = c(-0.5e6, 0.3e6), ylim = c(-0.5e6, 0.8e6),
            minlon = 0, orient = "v") 
my_aea_plot("arctic", c("ar", "gr", "eu", "na", "no"), arctic_proj, 
            xlim = c(-2.6e6, 2e6), ylim = c(-2.7e6, 2.3e6),
            minlat=66) 
my_aea_plot("high mtn asia", "hi", hma_proj, 
            xlim = c(-3e6, 2e6), ylim = c(-1.5e6, 2.5e6)) # HMA +caucasus
my_aea_plot("andes","pa",  andes_proj, 
            xlim = c(-1e6, 1e6), ylim = c(-4e6, 3e6),
            orient ="v") 
my_aea_plot("antarctic peninsula", "an", antarctica_proj, 
            xlim = c(-0.5e6, 1.85e6), ylim = c(-0.4e6, 1.6e6),
            orient = "v") 
my_aea_plot("new zealand", "nz", nz_proj, 
            xlim = c(-0.5e6, 0.7e6), ylim = c(-0.7e6, 1e6),
            orient = "v") 

  

## check for duplicate band values ----------------
# (this could happen when points randomly overlap on same pixel)
ndistinct <- train_full %>% 
  select(where(is.numeric)) %>%
  distinct() %>% 
  nrow() %>% print()

nrow(train_full) - ndistinct # 0




# compute normalized differences ----------------

nd <- function(a,b) (a-b)/(a+b)
# nd(1,-1) # returns Inf if 0 in denominator

# normalize each band by B2
# additional 19 features
train_nds <- train_full %>% 
  mutate(nd2_1 = nd(b2, b1),
         nd2_7 = nd(b2, b7),
         nd2_8 = nd(b2, b8),
         nd2_8a = nd(b2, b8a),
         nd2_9 = nd(b2, b9),
         nd2_11 = nd(b2, b11),
         nd2_12 = nd(b2, b12),
         nd4_8 = nd(b4, b8), # reverse NDVI
         nd5_4 = nd(b5, b4),
         nd5_6 = nd(b5, b6),
         nd4_6 = nd(b4, b6),
         nd5_3 = nd(b5, b3),
         nd6_3 = nd(b6, b3),
         nd4_6 = nd(b4, b6),
         nd4_2 = nd(b4, b2),
         nd5_2 = nd(b5, b2),
         nd6_2 = nd(b6, b2),
         rgnd = nd(b4, b3),
         gbnd = nd(b3, b2),
         ndsi = nd(b3, b11))
colnames(train_nds)


# !generate folds, splits--------------
set.seed(123)
split <- train_nds %>% 
  select(class, where(is.numeric), -lat, -lon) %>%  # set aside ID vars
  mutate(class = fct_relevel(class, "1")) %>% 
  initial_split(strata = class) # since algae are minority class, want to ensure a high number in the test set
train <- training(split)
test <- testing(split)
folds <- vfold_cv(train, strata = class)
rf_rec_raw <- recipe(class~., data = train)

rf_rec_raw %>% prep() 

# feature EDA ------------------------
## check for NZV or lincomb -------------------------
# any NZV or lincomb? 
rf_rec_raw %>% 
  step_nzv(all_numeric_predictors()) %>% 
  step_lincomb(all_numeric_predictors()) %>% 
  prep()
# no

# are there highly correlated features?
rf_rec_raw %>% 
  step_corr(all_numeric_predictors(), threshold = 0.9) %>% 
  prep() %>% juice() %>% ncol()
# removes 10 features. can tune on this 

# which features are most highly correlated?
cor_tbl <- rf_rec_raw %>% 
  prep() %>% 
  juice() %>% 
  select(where(is.numeric), -lon, -lat) %>% 
  cor() %>% 
  as_tibble(rownames = "features") %>% 
  pivot_longer(-features) %>% 
  filter(features > name) %>% # removing redundant pairs
  drop_na() %>% 
  arrange(desc(abs(value))) %>% 
  print(n = 20)

cor_tbl %>% 
  ggplot(aes(x = value)) + 
  geom_histogram(color = "white") +
  labs(x = "Pearson's correlation coef")




## variable importance--------------------
# what metric does permutation evaluate, accuracy?
# TO DO can we change this to evaluate our metric, PR AUC?
perm_rf_mod <- rand_forest(trees = 500) %>% 
  set_mode("classification") %>% 
  set_engine("ranger", importance = "permutation")


ff <- perm_rf_mod %>% 
  fit(class~., data = train)

ff %>% 
  vi_permute(target = "class", 
             metric = ~pr_auc(truth = actual, estimate = predicted),
             smaller_is_better = FALSE,
             pred_wrapper = function(object, newdata) predict(object, newdata))

fold_vi <- folds %>% 
  mutate(tt = map(splits, training),
         mod = map(tt, ~fit(perm_rf_mod, class~., data = .x)),
         vi = map(mod, vi)) %>% 
  unnest(vi) %>% 
  select(id, Variable, Importance)

fold_vi %>% 
  mutate(Variable = Variable %>% fct_reorder(Importance)) %>% 
  ggplot(aes(y = Importance, x = Variable)) +
  geom_boxplot() +
  labs(y = "Decrease in accuracy by permuting variable") +
  coord_flip()
ggsave(here("figs/map_qc/global_vi_permutation.png"))
# all features have positive importance, so not removing any features


all_features <- fold_vi %>% 
  mutate(Variable = Variable %>% fct_reorder(Importance)) %>% 
  pull(Variable) %>% 
  levels()

# # Use recursive feature elimination NOT WORKING
# rfe_model <- rand_forest(mode = "classification") %>% set_engine("ranger", importance = "permutation")
# rfe_rec <- recipe(class~., data = train) %>% 
#   step_select_vip(all_predictors(), outcome = "class", model = rfe_model, threshold = 0.9)
# # SLOW
# rfe_rec %>% 
#   prep() %>% 
#   juice()
# 




# tuning  --------------------------------------------

my_metrics <- metric_set(pr_auc)

rf_tune_spec <- rand_forest(trees=500, mtry = tune(), min_n = tune()) %>% # use 1000 trees for my final model
  set_mode("classification") %>% 
  set_engine("ranger", probability = TRUE)

rf_tune_wf <- workflow() %>%
  add_recipe(rf_rec_raw) %>%
  add_model(rf_tune_spec)

# # SLOW ~15 min
# set.seed(345)
# coarse_tune_res <- tune_grid(
#   rf_tune_wf,
#   resamples = folds,
#   metrics = my_metrics, 
#   grid = 10
# )
# 
# # coarse search results
# coarse_tune_res %>%
#   collect_metrics() %>%
#   select(mean, min_n, mtry, std_err) %>%
#   pivot_longer(min_n:mtry,
#                values_to = "value",
#                names_to = "parameter"
#   ) %>%
#   ggplot(aes(x = value, y = mean)) +
#   geom_point(show.legend = FALSE) +
#   geom_linerange(aes(ymax = mean + std_err,
#                      ymin = mean - std_err)) +
#   facet_wrap(~parameter, scales = "free_x") +
#   labs(x = NULL, y = "PR_AUC")
# 
# 
# # fine tune, based on coarse grid
# fine_grid <- grid_regular(
#   mtry(range = c(2,9)),
#   min_n(range = c(2,12)),
#   levels = 6
# )
# 
# # SLOW ~ 15 min
# set.seed(456)
# regular_res <- tune_grid(
#   rf_tune_wf,
#   resamples = folds,
#   metrics = my_metrics,
#   grid = fine_grid
# )
# 
# regular_res %>%
#   collect_metrics() %>%
#   mutate(min_n = factor(min_n)) %>%
#   ggplot(aes(x=mtry, y=mean, color = min_n)) +
#   geom_linerange(aes(ymax = mean + std_err,
#                      ymin = mean - std_err), alpha = 0.2, linewidth = 4) +
#   geom_line(alpha = 0.5, size = 1.5) +
#   geom_point() +
#   labs(y = "PR AUC")
# 
# regular_res %>% 
#   collect_metrics()
# 
# 
# # finalize the model
# best_auc <- select_best(regular_res, "pr_auc")
# 
# tuned_rf <- finalize_model(
#   rf_tune_spec,
#   best_auc
# )

# tuned_wf <- update_model(rf_tune_wf, tuned_rf)

# mtry = 3, min_n = 2

# shortcut to avoid re-running the tuning:

tuned_rf <- rand_forest(trees=500, mtry = 3, min_n = 2) %>% # use 1000 trees for my final model
  set_mode("classification") %>% 
  set_engine("ranger", probability = TRUE)

tuned_wf <- rf_tune_wf %>% 
  update_model(tuned_rf)








# thresholding-----------------------------
my_beta <- 1 # f measure parameter

control <- control_resamples(save_pred = TRUE)

# SLOW -- 1 min
rf_res <- fit_resamples(tuned_wf, folds, control = control)

# extract predictions and annotate with metadata

# hack to regenerate train set , this time with regino incl
set.seed(123)
split2 <- train_nds %>% 
  select(class, where(is.numeric), -lat, -lon, region) %>%  # set aside ID vars
  mutate(class = fct_relevel(class, "1")) %>% 
  initial_split(strata = class) # since algae are minority class, want to ensure a high number in the test set
train2 <- training(split2)

regions_col <- train2 %>%
  rownames_to_column() %>% 
  select(region, rowname) %>% 
  mutate(rowname=as.numeric(rowname))

preds <-
  rf_res %>% 
  unnest(.predictions) %>% 
  left_join(regions_col, by = c(".row"="rowname")) #%>% 
  # select(region, id, class, .pred_1)

# weighted harmonic mean, where precision is 2x more important then recall
f_beta <- new_class_metric(function(data, truth, estimate, ...){
  f_meas(
    data = data, 
    truth = {{truth}},
    estimate = {{estimate}},
    beta = my_beta # set the beta option here
  )
}, "maximize")

# coarse search
thresholds <- seq(0.1, 0.9, by = 0.05) # test a range of thresholds

global_thresh_coarse <- preds %>%  
  group_by(id) %>% 
  probably::threshold_perf(truth = class, estimate = .pred_1, thresholds = thresholds, metrics = metric_set(f_beta)) %>% 
  group_by(.threshold) %>% 
  summarise(mean = mean(.estimate), ci = 1.96*sd(.estimate))
global_thresh_coarse %>% 
  ggplot(aes(x = .threshold, y = mean)) +
  geom_ribbon(aes(ymin = mean -ci, ymax = mean+ci), fill = "blue", alpha = 0.2) +
  geom_line() +
  labs(y = paste0("F", my_beta))
ggsave(here("figs/map_qc/thresholds.png"))

fine_thresholds <- seq(0.3, 0.6, by = 0.01)

global_thresh_fine <- preds %>%  
  group_by(id) %>% 
  probably::threshold_perf(truth = class, estimate = .pred_1, thresholds = fine_thresholds, metrics = metric_set(f_beta)) %>%
  group_by(.threshold) %>% 
  summarise(mean = mean(.estimate), ci = 1.96*sd(.estimate))  %>% 
  arrange(-mean) 
global_thresh_fine

global_thresh_fine %>% 
  ggplot(aes(x = .threshold, y = mean)) +
  geom_ribbon(aes(ymin = mean -ci, ymax = mean+ci), fill = "blue", alpha = 0.2) +
  geom_line() +
  labs(y = paste("F", my_beta), x = "Threshold")

global_opt_thresh <- global_thresh_fine %>% 
  filter(mean == max(mean)) %>% 
  pull(.threshold)

regions3 <- regions2 %>% 
  add_column(threshGlobF1 = round(global_opt_thresh,2))


## F1 regional thresholding---------------------

regional_thresholds <- seq(0.2, 0.8, by = 0.01)
regional_cv_threshold_res <- preds %>%  
  group_by(id, region) %>% 
  probably::threshold_perf(truth = class, estimate = .pred_1, thresholds = regional_thresholds, metrics = metric_set(f_beta)) %>%
  group_by(.threshold, region) %>% 
  summarise(mean = mean(.estimate), ci = 1.96*sd(.estimate)) %>%
  ungroup() 

opt_thresh_regional <- regional_cv_threshold_res %>% 
  group_by(region) %>% 
  filter(mean == max(mean, na.rm = TRUE)) %>% 
  slice_tail(n=1) %>% # choose the higher of the two thresholds if tie
  ungroup()
opt_thresh_regional 

ll <- opt_thresh_regional %>% 
  select(region, .threshold) %>% 
  mutate(threshRegF1 = round(.threshold,2), .keep="unused") 

regions4 <- regions3 %>% left_join(ll, by = c("thresholdID" = "region"))
  

regional_cv_threshold_res %>% 
  # drop_na() %>% 
  ggplot(aes(x = .threshold, y = mean)) +
  geom_ribbon(aes(ymin = mean -ci, ymax = mean+ci), fill = "blue", alpha = 0.2) +
  geom_line() +
  labs(y = paste("F", my_beta)) +
  facet_wrap(vars(region)) +
  geom_point(data = opt_thresh_regional, aes(x = .threshold, y = mean))


## f0.5 regional thresholding -------------------------------

f05_beta <- new_class_metric(function(data, truth, estimate, ...){
  f_meas(
    data = data, 
    truth = {{truth}},
    estimate = {{estimate}},
    beta = 0.5 # set the beta option here
  )
}, "maximize")

f05_regional_cv_threshold_res <- preds %>%  
  group_by(id, region) %>% 
  probably::threshold_perf(truth = class, estimate = .pred_1, thresholds = regional_thresholds, metrics = metric_set(f05_beta)) %>%
  group_by(.threshold, region) %>% 
  summarise(mean = mean(.estimate), ci = 1.96*sd(.estimate)) %>%
  ungroup() 

f05_opt_thresh_regional <- f05_regional_cv_threshold_res %>% 
  group_by(region) %>% 
  filter(mean == max(mean, na.rm = TRUE)) %>% 
  slice_tail(n=1) %>% # choose the higher of the two thresholds if tie
  ungroup()
f05_opt_thresh_regional 

ll2 <- f05_opt_thresh_regional %>% 
  select(region, .threshold) %>% 
  mutate(threshRegF5 = round(.threshold,2), .keep="unused") 

regions5 <- regions4 %>% left_join(ll2, by = c("thresholdID" = "region"))

# fill in missing with default 0.5
regions6 <- regions5 %>% 
  replace_na(list(threshRegF1 = 0.5, threshRegF5 = 0.5))

st_write(regions6, here("data/s2_classifier_map/regions_plus_thresholds.shp"))
regions6 %>% 
  as_tibble() %>% 
  select(-geometry) %>% 
  write_csv(here("data/s2_classifier_map/regions_plus_thresholds.csv"))


# generate final metrics ------------------

lf <- last_fit(tuned_wf, split, control = control)

best_threshs <- opt_thresh_regional %>% 
  select(region, .threshold) %>% 
  bind_rows(
    default_threshs <- tribble(
      ~region, ~.threshold,
      "an",0.5,
      "hi",0.5,
      "nz", 0.5
    )
  ) %>% 
  arrange(.threshold)


# region meta
# (lf uses the .row from the full set)
regions_col_full_set <- train_nds %>% 
  rownames_to_column() %>% 
  select(region, rowname) %>% 
  mutate(rowname = as.numeric(rowname))


preds <- lf %>%  
  unnest(.predictions) %>% 
  select(.row, class, .pred_1) %>% 
  left_join(regions_col_full_set, by = c(".row"= "rowname")) %>% 
  select(-.row) %>% 
  left_join(best_threshs) %>% 
  mutate(.pred_class = if_else(.pred_1>.threshold, as.factor("1"), as.factor("0")))  # replace default pred class with our custom threshold
  


final_metrics <- metric_set(precision, recall)
final_metrics(preds, truth = class, estimate = .pred_class)

# 
# # A tibble: 2 × 3
# .metric   .estimator .estimate
# <chr>     <chr>          <dbl>
# 1 precision binary         0.932
# 2 recall    binary         0.928

# of course this metric does not take into account the potential for unequal training representation per region... 
# but that;s for the validation set ...




# clustering pts, band values -----------------------------

# density plots of band values
train_distinct %>% 
  pivot_longer(cols = b1:b9) %>% 
  ggplot(aes(value, color = class_name)) +
  geom_density() +
  facet_wrap(vars(name), scales = "free") +
  scale_color_manual(values=classname_colors)


# density plot of ND data
train_nds %>% 
  select(class_name, where(is.numeric),-lat,-lon, -b1:-b9) %>% 
  pivot_longer(-class_name) %>% 
  ggplot(aes(value, color = class_name)) +
  geom_density() +
  facet_wrap(vars(name), scales = "free") +
  scale_color_manual(values=classname_colors)

umap_fitted <- train_nds %>% 
  select(where(is.numeric)) %>% 
  scale() %>% 
  umap()

umap_tbl <- train_nds %>% 
  cbind(umap_fitted$layout) %>% # old school augment
  rename(umap1 = `1`, umap2 = `2`) %>% 
  as_tibble()


umap_tbl %>% 
  ggplot(aes(x = umap1, y = umap2, color = class_name)) +
  geom_point() +
  scale_color_manual(values=classname_colors)

umap_tbl %>% 
  group_by(random) %>% 
  slice_sample(n=10) %>% 
  ungroup() %>% 
  ggplot(aes(x = umap1, y = umap2, color = class_name)) +
  geom_text(aes(label = random)) +
  scale_color_manual(values=classname_colors)



# export ---------------

# export_id <- "trainclean3"
# train_nds %>% 
#   select(class, where(is.numeric), -lat, -lon) %>% 
#   rename_with(toupper) %>% # format for GEE
#   write_csv(paste0(here("data/s2_classifier_map/training/"), export_id, ".csv"))
# 
# # globalSnowAlgae/2_trainDat/
# 
# 
# train_nds %>% 
#   select(basename, class, region, hrs, class_name, lon, lat) %>% 
#   write_csv(paste0(here("data/s2_classifier_map/training/"), "meta_", export_id, ".csv"))
# 


# clean train 2 --removed algae w sus dust comtam

# scrach  ---------------------



# iteratively remove features ---------------------------

# start with the lowest VI
worst_vars <- all_features %>% head(n = 5)
# TO DO all combinations of the bad variables
bad_var_combs <- worst_vars 
# 
set.seed(123)
v <- 3
# rep_folds <- 
vfold_cv(train, v = v, strata = class, repeats = length(bad_var_combs)) %>% 
  add_column(vars  = rep(bad_var_combs, each = v)) %>% 
  mutate(splits = map(splits, ~training(.x) %>% select(-vars)))



set.seed(123)
new_folds <- train %>% 
  select(-worst_vars[1]) %>% 
  
  
  new_res <- fit_resamples(tuned_wf, new_folds)

# tuning... ok to use full splits to tune, including thresholding, since this data is not used for error estimation 
# feature selection


# generate probability predictions




# compare with region-specific thresholds----------
# generate probability predictions
set.seed(123) # try with different seeds (too lazy to set up resampling)
folds <- train_nds %>% 
  filter(region=="na") %>% 
  mutate(class = fct_relevel(class, "1")) %>% # ensure "1" is the "event" level
  vfold_cv()


rf_mod <- rand_forest(trees=500, mtry = 2, min_n = 3) %>% 
  set_mode("classification") %>% 
  set_engine("ranger", probability = TRUE)

rf_rec <- recipe(class~., data = train_nds %>% filter(region=="na")) %>% 
  update_role(basename,new_role = "id")  %>%
  update_role(region, new_role = "id") %>% 
  update_role(hrs, new_role = "id") %>% 
  update_role(class_name, new_role = "id") %>%  
  update_role(lon, new_role = "id") %>% 
  update_role(lat, new_role = "id") 

rf_wf <- workflow() %>%
  add_recipe(rf_rec) %>%
  add_model(rf_mod)

control <- control_resamples(save_pred = TRUE)

# SLOW
rf_res2 <- fit_resamples(rf_wf, folds, control = control)

# extract predictions
preds2 <- rf_res2 %>% 
  unnest(.predictions) %>% 
  select(id, class, .pred_1)



# thresholding
library(probably)

# coarse search
thresholds <- seq(0.1, 0.9, by = 0.05) # test a range of thresholds
preds2 %>%  
  probably::threshold_perf(truth = class, estimate = .pred_1, thresholds = thresholds, metrics = metric_set(f_meas)) %>% 
  ggplot(aes(x = .threshold, y = .estimate)) +
  geom_line()


fine_thresholds <- seq(0.3, 0.6, by = 0.01)
preds2 %>%  
  probably::threshold_perf(truth = class, estimate = .pred_1, thresholds = fine_thresholds, metrics = metric_set(f_meas)) %>%
  arrange(-.estimate) %>% #0.41
  ggplot(aes(x = .threshold, y = .estimate)) +
  geom_line()


