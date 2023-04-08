#*
#*
#* This script processes the algal cover area and validation data output from GEE
#*
#* Includes QC graphs and stats relating to validation data coverage
#*
#* 
#* Following the approach of Murray et al 2022
#*
#* to combine regions of error, group by decile and summarise decile area, then treat as single dataset
#* but the sampling effort is not even between regions. would need to additionally weight as proportion of TOTAL region area, 
#* eg caucasus+HMA TOTAL region area, weight caucasus decile area out of TOTAL 
#* 
#* 
#* TO DO:
#* 
#*

library(tidyverse)
library(here)
library(janitor)
library(fs)
library(tidymodels)
library(sf) 
library(officer)
library(flextable)
library(snakecase)
library(kableExtra)

theme_set(theme_bw())
# read data into R ---------------------


validation_raw <- dir_ls(here("data/s2_classifier_map/validation/"), regexp = "csv") %>% 
  map_df(read_csv, id = "path") 

# decile areas
decile_area_raw <- dir_ls( here("data/s2_classifier_map/area/per_region_F1"), regexp = "T.kml$") %>% 
  map_df(st_read) %>% 
  distinct(region, .keep_all=TRUE) # remove duplicates

# algae areas
algae_areas_raw <- dir_ls(here("data/s2_classifier_map/area/"), regexp="thrs") %>% 
  map_df(read_csv, id = "path")

thresholds <- read_csv(here("data/s2_classifier_map/regions_plus_thresholds.csv")) %>% 
  clean_names()


# wrangle  ------------------------------

# check for duplicate points
nrow(validation_raw) # 1541
length(unique(validation_raw$.geo)) # 1541


threshs <- thresholds %>% 
  rename(region = name) %>% 
  select(-threshold_id) %>% 
  mutate(region = if_else(region=="antarcticaNoSheet", "antarctica",region),
         region = if_else(str_detect(region,"Greenland"), "greenland",region),
         region = if_else(str_detect(region, "norway"), "norway", region) # since we didnt collect norwayNorth or South spec error
  ) %>% 
  bind_rows(tribble(
    ~region, ~subregion_of,~thresh_glob_f1,~thresh_reg_f1, ~thresh_reg_f5,
    "arctic", "arctic",    0.39,             0.57,          0.57,
    "northAmerica", "northAmerica", 0.39,   0.2,          0.35,
    "europe", "europe",    0.39,            0.27,            0.55)) %>% 
  distinct() 

# get area per decile-region for weighting
decile_areas <- decile_area_raw %>% 
  janitor::clean_names() %>% 
  select(region, starts_with("x")) %>% 
  mutate(across(c(matches("d")), as.numeric )) %>% 
  pivot_longer(cols = starts_with("x")) %>% 
  mutate(decile = paste0(str_sub(name, 2,2)) %>% as.factor(), .keep="unused") %>% 
  rename(area_km2 = value) %>% 
  filter(region!="antarcticaNoSheet") %>%  # not used
  as_tibble() %>% 
  select(-geometry) %>% 
  left_join(select(threshs, region, subregion_of))
count(decile_areas, subregion_of, region) %>% print(n=35)

# clean the validation data
validation1 <- validation_raw %>% 
  janitor::clean_names() %>% 
  mutate(truth = a_true_class, 
         region = str_split_fixed(basename(path), "_", 3)[,1], 
         coords = str_extract(geo, "(?<=\\[)[:graph:]+(?=\\])"),#[:digit:]+ (?=\\,)
         lon = str_split_fixed(coords, ",", 2)[,1] %>% parse_number(),
         lat = str_split_fixed(coords, ",", 2)[,2] %>% parse_number(),
         truth = as.factor(truth),
         decile = as.factor(decile),
         .keep="unused") %>% 
  filter(is.na(c_training_image) ) %>% # data used for training should not be considered in validation
  select(region, decile, truth, class_prob, lon, lat, b_date_of_true_class) %>% 
  rownames_to_column("pt_id") 



## algae area data -----------------------

algae_areas1 <- algae_areas_raw %>% 
  mutate(threshname = basename(path) %>% str_extract("(?<=thrs)[:alnum:]{3}")) %>% 
  mutate(name = if_else(str_detect(name, "norway"), "norway", name),
         name = if_else(str_detect(name, "Greenland"), "greenland", name)) %>% # case sensitive, not incl greenlandNoSheet 
  rename(region = name) %>% 
  select(threshname, region, algae, glacier, redmask) %>% 
  left_join(select(threshs, region, subregion_of), by = c("region")) %>% 
  filter(region!="antarcticaNoSheet") %>% 
  group_by(threshname, region, subregion_of) %>% 
  summarise(algae= sum(algae), glacier = sum(glacier), redmask = sum(redmask)) %>% 
  ungroup()

# need areas of arctic, europe, northAmerica

algae_areas_select_continents <- algae_areas1 %>% 
  filter(subregion_of %in% c("arctic", "europe", "northAmerica", "highMtnAsia")) %>% 
  group_by(threshname, subregion_of) %>% 
  summarise(algae= sum(algae), glacier = sum(glacier), redmask = sum(redmask)) %>% 
  ungroup()  %>% 
  mutate(region = subregion_of)
  
world_area <- algae_areas1 %>% 
  filter(region!= "greenlandNoSheet") %>% #  dont double count greenland
  group_by(threshname) %>% 
  summarise(algae= sum(algae), glacier = sum(glacier), redmask = sum(redmask)) %>% 
  add_column(region = "world", subregion_of = "world") %>% 
  ungroup()

algae_areas <- algae_areas1 %>% 
  filter( region != "highMtnAsia") %>% # report highMtnAsia *including* the Altai
  bind_rows(algae_areas_select_continents) %>% 
  bind_rows(world_area)

count(algae_areas, region, subregion_of) %>% 
  arrange(n)

## assign binary CLASS using various thresholds -----------
validation <- validation1 %>% 
  left_join(threshs) %>% 
  mutate(class_gf1 = if_else(class_prob>thresh_glob_f1, 1, 0),
         class_rf1 = if_else(class_prob>thresh_reg_f1, 1, 0),
         class_rf5= if_else(class_prob>thresh_reg_f5, 1, 0)) %>% 
  mutate(subregion_of = if_else(is.na(subregion_of), region, subregion_of))

validation %>% distinct(region, subregion_of) %>% arrange(subregion_of)
nrow(validation) # 1541
validation %>% distinct(decile, subregion_of, region, pt_id, truth, class_gf1, class_rf1, class_rf5) %>% nrow() # 1541--no dups


## export val data GEE validation app ------------
validation %>% 
  arrange(region, -class_prob) %>% 
  rename(date = b_date_of_true_class) %>% 
  drop_na(truth) %>% # removes 2 pts
  as_tibble() %>% 
  mutate(truth_chr = truth %>%
           fct_recode( "Algae predominant" = "1",
                       "No algae"="0",
                       "Algae maybe present"="2",
                       "Algae maybe present"="3",
                       "Unsure"="4")) %>% 
  select(pt_id, region, decile, truth, truth_chr, class_prob, class_rf1, class_gf1, class_rf5, lon, lat, date, subregion_of) %>% 
  group_by(region) %>%
  nest() %>% 
  mutate(data = map(data, ~rownames_to_column(.x, "region_pt_id") )) %>% 
  unnest(cols = c(data)) %>%
  ungroup() %>%
  write_csv(here("data/s2_classifier_map/validation/final_combined/validation_final.csv"))




# plot validation barplot by decile-region -----------------------------

# exclude caucasus
my_ord <- c("northAmerica","greenland","europe","arctic","highMtnAsia","kamchatka", "andes","antarctica","newZealand")

# calculate percentage breakdown of labels per decile-region
continent_decile_plot_data <- validation %>% 
  drop_na(truth) %>% 
  mutate(truth_chr = truth %>%
           fct_collapse("Maybe present" = c("2", "3")) %>% # collapse 2 and 3 codes 
           fct_recode( "Predominant" = "1",
                       "Not present"="0",
                       "Unsure"="4") ) %>% 
  # because subregions are sampled unevenly, weight by decile-region area
  left_join(decile_areas) %>% 
  group_by(subregion_of, decile, truth_chr) %>%
  summarise(n_wtd = sum(area_km2)) %>% 
  group_by(subregion_of, decile) %>%
  mutate(percent = 100* (n_wtd/sum(n_wtd) ) ) %>% 
  ungroup() %>% 
  mutate(subregion_of = fct_relevel(subregion_of, my_ord),
         truth_chr = fct_relevel(truth_chr, c("Predominant", "Maybe present","Not present", "Unsure"))  )  
    
# annnotate sample size and area in plot
cont_decile_areas <- decile_areas %>% 
  group_by(subregion_of, decile) %>% 
  summarise(area_km2 = sum(area_km2)) %>% 
  ungroup() %>% 
  mutate(area_short = case_when(area_km2>999 ~ formatC(area_km2, format = "e", digits = 0) %>% str_remove("\\+0"),
                                area_km2>1 ~ formatC(area_km2, format = "d", digits = 0),
                                area_km2<1 ~ "<1")) %>% 
  left_join( count(validation, subregion_of, decile) ) %>% 
  mutate(subregion_of = fct_relevel(subregion_of, my_ord))
count(cont_decile_areas, subregion_of) # should be 10 deciles per continent

#              red         orange       teal       grey
err_colors <- c("#fb8072", "#fdb462", "#8dd3c7", "grey50")


ggplot(continent_decile_plot_data, aes(x=decile, y=percent)) +
    geom_bar(stat="identity", aes(fill = truth_chr)) +
    facet_wrap(vars(subregion_of), ncol=2) +
    scale_fill_manual(values = err_colors) +
    # scale_x_discrete(breaks = seq(0,10,1), limits = c(-5,100)) +
    scale_y_continuous(breaks = c(0,50,100), minor_breaks = seq(0,100,10), limits = c(-5,120)) +
    labs( x = "Classifier probability decile", y = "% Label", fill = "Validation label:\nSnow algae...") +
    theme_minimal() +
    geom_text(data=cont_decile_areas, aes(x = decile, y = 100, label = paste0(area_short,"\n",n)), vjust = 0, size = 3) #+
  # annotate("text", label = "A=       \nn=       ", x = as.factor("0"), y = 115, size = 3)

# ggsave(here("figs/distribution/2_validation_error/validation_barplot.pdf"), width = 8, height = 10, units = "in", dpi = 300)



# estimate error ------------------------------
# this is the error of algae predominant vs not predominant-- green and red on the Fig S2. not including the orange and grey data (partial coverage, or )
my_metrics <- metric_set(precision, recall)

# discard intermediate or ambiguous labelled points
val_01 <- validation %>% 
  filter(truth %in% as.factor(0:1)) %>% # discard 2,3,4
  mutate(truth = as.character(truth) %>% fct_relevel("1"),# must relevel to remove higher order factor levels
         class_rf1 = as.factor(class_rf1) %>% fct_relevel("1"),
         class_gf1 = as.factor(class_gf1) %>% fct_relevel("1"),
         class_rf5 = as.factor(class_rf5) %>% fct_relevel("1")) %>% 
  left_join(decile_areas) 
nrow(val_01)

count(val_01, subregion_of, region) %>% arrange(-n)

my_boot_error <- function(df, grp = region, class = class_rf1, nboots = 10, weight = area_km2, strat = decile){
  df %>% 
    group_by({{grp}}) %>% 
    nest() %>% 
    mutate(boots = map(data, bootstraps, times = nboots, strata = {{strat}}), .keep="unused") %>% 
    unnest(boots) %>% 
    mutate(training = map(splits, training),
           metrics = map(training, ~my_metrics(.x,  truth = truth, estimate = {{class}}, case_weights = {{weight}})),
           .keep="unused") %>% 
    unnest(metrics) %>% 
    group_by({{grp}}, .metric) %>% 
    summarise(mean = mean(.estimate), 
              ci = 1.96*sd(.estimate),
              min = if_else(mean-ci<0, 0, mean-ci)) %>% 
    ungroup()
}




ord <- c("northAmerica",    "alaskaPeninsula", "alaskaRange",     "coastNorth" ,     "coastSouth",
         "interiorNorth",   "interiorSouth",   "cascades",        "greenland",       "greenlandNoSheet",
         "europe",
         "iceland",         "norway",          "alps",            "arctic",          "svalbard",
         "highMtnAsia",     "altai",           "caucasus",        "kamchatka",       "andes",    "newZealand",
         "antarctica", "world" )

# uncomment this section to re-run the bootstrap analysis
# runtime: ~ 30 min


# ##  regional F1 thresh error ----------------------------
# bb <- 1000
# 
# # runtime: ~ 5 min for 500 replicates
# set.seed(123)
# start <- Sys.time()
# subregion_err <- val_01 %>% 
#   filter(!region %in% c("europe", "northAmerica")) %>% # calculate these seperately, because they are composed of the aggregate of multiple 
#   my_boot_error(nboots = bb)
# print( Sys.time() - start )
# 
# set.seed(123)
# cont_err <- val_01 %>% 
#   filter(subregion_of %in% c("europe","northAmerica")) %>% 
#   my_boot_error(grp = subregion_of, nboots = bb) %>% 
#   rename(region = subregion_of)
# 
# set.seed(123)
# world_err <- val_01 %>% 
#   my_boot_error(grp=NULL, nboots = bb) %>% 
#   add_column(region="world")
# 
# f1_err <- subregion_err %>% 
#   bind_rows(cont_err) %>% 
#   bind_rows(world_err)%>% 
#   add_column(threshname = "RF1")
# 
# ## table 
# 
# f1_err_tt <- f1_err %>% 
#   select(-min, -threshname) %>% 
#   mutate(across(where(is.numeric), ~round(.x, digits = 2))) %>% 
#   unite(estimate, mean, ci, sep = " ± ") %>% 
#   pivot_wider(names_from = ".metric", values_from = "estimate") %>% 
#   mutate(region = fct_relevel(region, ord)) %>% 
#   arrange(region)
# 
# f1_err_tt %>% 
#   mutate(region = to_title_case(as.character(region))) %>% 
#   flextable() %>% 
#   save_as_docx( path = here("figs/distribution/rf1_error.docx"))
# 
# ## global F1 thresh err  -----------------------------------
# set.seed(123)
# gf1_subregion_err <- val_01 %>% 
#   filter(!region %in% c("europe", "northAmerica")) %>% 
#   my_boot_error(class = class_gf1, nboots = bb)
# 
# set.seed(123)
# gf1_cont_err <- val_01 %>% 
#   filter(subregion_of %in% c("europe","northAmerica")) %>% 
#   my_boot_error(class = class_gf1, grp = subregion_of, nboots = bb) %>% 
#   rename(region = subregion_of)
# 
# set.seed(123)
# gf1_world_err <- val_01 %>% 
#   my_boot_error(class = class_gf1, grp=NULL, nboots = bb) %>% 
#   add_column(region="world")
# 
# gf1_err <- gf1_subregion_err %>% 
#   bind_rows(gf1_cont_err) %>% 
#   bind_rows(gf1_world_err)%>% 
#   add_column(threshname = "GF1")
# 
# ## regional F0.5 thresh err  -----------------------------------
# set.seed(123)
# rf5_subregion_err <- val_01 %>% 
#   filter(!region %in% c("europe", "northAmerica")) %>% 
#   my_boot_error(class = class_rf5, nboots = bb)
# 
# set.seed(123)
# rf5_cont_err <- val_01 %>% 
#   filter(subregion_of %in% c("europe","northAmerica")) %>% 
#   my_boot_error(class = class_rf5, grp = subregion_of, nboots = bb) %>% 
#   rename(region = subregion_of)
# 
# set.seed(123)
# rf5_world_err <- val_01 %>% 
#   my_boot_error(class = class_rf5, grp=NULL, nboots = bb) %>% 
#   add_column(region="world")
# 
# rf5_err <- rf5_subregion_err %>% 
#   bind_rows(rf5_cont_err) %>% 
#   bind_rows(rf5_world_err) %>% 
#   add_column(threshname = "RF5")
# 
# 
# ## calculate min and max area ---------------
# 
# error_full <- f1_err %>% 
#   bind_rows(gf1_err) %>% 
#   bind_rows(rf5_err) %>% 
#   left_join(algae_areas) 
# 
# error <- error_full %>% 
#   select(threshname, region, algae, glacier, .metric, min) %>% 
#   pivot_wider(names_from = .metric, values_from = min) %>% 
#   rename(min_precision = precision, min_recall = recall) %>% 
#   mutate(min_algae = algae * min_precision,
#          max_algae = if_else(is.na(min_recall), NA, 
#                              algae/min_recall ) ) %>% 
#   select(-min_precision, -min_recall) %>% 
#   mutate(percent = algae/glacier,
#          min_p = min_algae/glacier,
#          max_p = max_algae/glacier) %>% 
#   relocate(glacier, .after = max_algae)
# 
# ## compare estimates across thresholds --------------
# 
# error %>% write_csv(here("data/s2_classifier_map/error.csv"))

# import the pre-computed bootstrap stats
error <- read_csv(here("data/s2_classifier_map/error.csv"))

my_continents <-  c("northAmerica", "greenland","iceland","norway", "alps","arctic","highMtnAsia","kamchatka","andes","newZealand","antarctica","world")
error %>% 
  mutate(threshname = case_when(threshname=="RF1"~"Region-spec F1",
                                threshname=="GF1"~"Global F1",
                                threshname=="RF5"~"Region-spec F0.5"),
         percent = percent*100, min_p = min_p*100, max_p = max_p*100,
         region = fct_relevel(region, my_continents)) %>% 
  filter(region %in% my_continents) %>% 
  ggplot(aes(x = threshname, y= percent, color = threshname)) +
  geom_point(position = position_dodge(width = 1)) +
  geom_errorbar(position = "dodge", stat = "identity", aes(ymin = min_p, ymax= max_p), width =0.2) +
  labs(y = "Percent Algal Cover", color = "Threshold")+
  scale_y_continuous(limits = c(0,NA))+ #breaks = seq(0,12,1), minor_breaks = seq(0,12,1), 
  facet_wrap(vars(region), scales="free") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())#,
# legend.position = "bottom")
ggsave(here("figs/distribution/area_tables/threshold_compare.png"), width = 8, height = 5, units="in")
# manually add in the percents for low area cover regions

## table of F1 area ----------------

round100 <- function(x){
  case_when(x>10000 ~ round(x,-2),
            x>1000 ~ round(x,-1),
            # x>100 ~ round(x,-1),
            .default = round(x, 0))
}

table_s3_data <- error %>% 
  filter(threshname=="RF1") %>% 
  select(-threshname) %>% 
  # set indentation level
  left_join(select(threshs, region, subregion_of )) %>% 
  mutate(type = case_when(region==subregion_of~ "continent",
                          region=="greenlandNoSheet" ~ "range",
                          region %in% c("world", "northAmerica","europe") ~"continent",
                          .default = "range") ) %>% 
  # round
  mutate(across(algae:max_algae, ~round100(.x))) %>% 
  mutate(glacier = round(glacier, -1)) %>% 
  mutate(across(percent:max_p, ~round(.x*100, 1) %>% as.character())) %>% # set as character to avoid putting commas in 
  replace_na(list(max_p = "-")) %>% 
  # sort in desired order `ord`
  mutate(region = fct_relevel(region, ord)) %>% 
  arrange(region) %>% 
  # prettify for table
  mutate(region =  to_title_case(as.character(region))) %>% 
  mutate(region = if_else(type=="range", paste("- ",region), region),
         region = if_else(region %>% str_detect("No Sheet"), "-  excluding sheet", region)) %>% 
  select(-type, -subregion_of)  %>% 
  mutate(across(c(percent, min_p, max_p), ~as.character(.x) ))


colnames(table_s3_data) <- c("Region", "Algae km2", "Min algae km2","Max algae km2", "Glacier km2", "Algae %", "Min %", "Max %")


table_s3_data %>% 
  flextable() %>% 
  colformat_double(big.mark = ",", digits = 0, na_str = "-") %>% 
  save_as_docx(path = here("figs/distribution/3_area/rf1_area.docx"))


# Fig 2, panel D -------------

svalbard_algal_area <- algae_areas %>% 
  filter(threshname=="RF1", region=="svalbard") %>% 
  select(-threshname, -subregion_of, -redmask)

svalbard_row <- error_full %>% 
  filter(threshname=="RF1", region=="arctic") %>% 
  select(.metric, min) %>% 
  pivot_wider(names_from = .metric, values_from = min) %>% 
  rename(min_precision = precision, min_recall = recall) %>% 
  add_column(svalbard_algal_area) %>% 
  mutate(min_algae = algae * min_precision,
         max_algae = if_else(is.na(min_recall), NA, 
                             algae/min_recall ) ) %>% 
  select(-min_precision, -min_recall) %>% 
  mutate(percent = algae/glacier,
         min_p = min_algae/glacier,
         max_p = max_algae/glacier) %>% 
  relocate(glacier, .after = max_algae)

fig2d <- error %>% 
  filter(threshname=="RF1", region %in% c("northAmerica", "greenland", "norway", "kamchatka")) %>% 
  select(-threshname) %>% 
  bind_rows(svalbard_row) %>% 
  mutate(region = fct_relevel(region, ord)) %>% 
  arrange(region) %>% 
  mutate(across(algae:max_algae, ~round100(.x))) %>%
  mutate(glacier = round(glacier, -1)) %>% 
  mutate(across(percent:max_p, ~round(.x*100, 1))) %>%
  mutate(region = if_else(region=="greenlandNoSheet", "greenland", region)) %>% 
  mutate(region = fct_relevel(region, "northAmerica", "greenland", "svalbard","norway", "kamchatka")) %>% 
  arrange(region) %>% 
  select(-glacier) %>% 
  unite(ci, min_algae, max_algae, sep = " to ") %>% 
  unite(ci_p, min_p, max_p, sep = " to ") %>% 
  mutate(ci  = paste0("(",ci,")")) %>% 
  mutate(ci_p  = paste0("(",ci_p,")")) %>% 
  mutate(region =  to_title_case(as.character(region))) %>% 
  unite(algae, algae, ci, sep = " ") %>% 
  unite(percent, percent, ci_p, sep = " ")
fig2d

flextable(fig2d) %>% 
  save_as_docx(path = here("figs/distribution/fig2d.docx"))
# format in docx for consistancy







