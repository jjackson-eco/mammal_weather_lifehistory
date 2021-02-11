####################################################
##                                                ##
##     Global climate and population dynamics     ##
##                                                ##
##             Naive weather effects              ##
##                                                ##
##                Dec 9th 2020                    ##
##                                                ##
####################################################

# Record-wise regressions linking weather to population growth rates, completely naive i.e.
# just the weather effect

rm(list = ls())
options(width = 100)

library(tidyverse)
library(viridis)
library(patchwork)


##__________________________________________________________________________________________________
#### 1. Load data ####

# mammal data
load("../rawdata/mammal.RData")
glimpse(mammal)

# annual weather anomaly - focus on just the mean anomaly in this script at a 5km range
mam_chelsa_annual <- readRDS("data/mam_chelsa_annual.RDS") %>% 
  filter(scale == "scale_5km") %>% 
  dplyr::select(ID,year, weather_scale = scale, mean_temp_anomaly, mean_precip_anomaly)
glimpse(mam_chelsa_annual)

# Species names to merge
load("../rawdata/GBIF_species_names_mamUPDATE.RData", verbose = TRUE)

##__________________________________________________________________________________________________
#### 2. Joining data ####

# linking to weather data and species names 
mammal_weather <- mammal %>% 
  left_join(., y = mam_chelsa_annual, by = c("ID", "year")) %>% 
  left_join(., y = dplyr::select(lpd_gbif, Binomial, gbif_species = gbif.species.tree),
            by = "Binomial") %>% 
  mutate(year_s = as.numeric(scale(year)))

glimpse(mammal_weather)


##__________________________________________________________________________________________________
#### 3. Models for each record ####

pgr_weather_naive <- mammal_weather %>% 
  group_by(ID_block) %>% 
  group_modify(~{
    
    # Temperature
    mod_temp = lm(pop_growth_rate ~ mean_temp_anomaly, data = .)
    coef_tempmod = coef(mod_temp)
    
    # Precipitation + dealing with NA values
    if(length(which(is.na(.$mean_precip_anomaly) == T)) == 0){
      mod_precip = lm(pop_growth_rate ~ mean_precip_anomaly, data = .)
      coef_precipmod = coef(mod_precip)}
    else{coef_precipmod = rep(NA,4)}     # Arbitrary long NA vector
    
    tibble(.[1,],
           coef_temp = unname(coef_tempmod[2]),
           coef_precip = unname(coef_precipmod[2]),
           n_obs = nrow(.))
  }) 

glimpse(pgr_weather_naive)

##__________________________________________________________________________________________________
#### 4. Saving data from GAM ####

mnanom_5km_naive <- pgr_weather_naive
save(mnanom_5km_naive, file = "data/pgr_weather/mnanom_5km_naive.RData")


