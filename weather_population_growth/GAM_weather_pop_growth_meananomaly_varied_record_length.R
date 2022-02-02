###########################################################
##                                                       ##
##        Global climate and population dynamics         ##
##                                                       ##
##  GAMM timeseries method - Varying timeseries lengths  ##
##                                                       ##
##                    Feb 2nd 2022                       ##
##                                                       ##
###########################################################

# Record-wise regressions linking weather to population growth rates, 
# accounting for autocorrelation with GAMM and formal AR(1) time-series analysis
# Models just for 5km buffer radius and the mean weather anomaly.

rm(list = ls())
options(width = 100)

library(tidyverse)
library(ggridges)
library(viridis)
library(patchwork)
library(gridExtra)
library(mgcv)         # Simon Wood to the rescue again. All Hail

temp_colour <- "#d45371"
precip_colour <- "#30738e"

##__________________________________________________________________________________________________
#### 1. Load data ####

# mammal data
load("../rawdata/mammal_alternate_lengths.RData", verbose = TRUE)
glimpse(mammal5)
glimpse(mammal20)

# annual weather anomaly - focus on just the mean anomaly in this script at a 5km range
mam_chelsa_annual <- readRDS("data/mam_chelsa_annual.RDS") %>% 
  filter(scale == "scale_5km") %>% 
  dplyr::select(ID,year, weather_scale = scale, mean_temp_anomaly, mean_precip_anomaly)
glimpse(mam_chelsa_annual)

# Species names to merge
load("../rawdata/GBIF_species_names_mamUPDATE.RData", verbose = TRUE)

##__________________________________________________________________________________________________
#### 2. Short data ####

# linking to weather data and species names 
mammal_weather_5 <- mammal5 %>% 
  left_join(., y = mam_chelsa_annual, by = c("ID", "year")) %>% 
  left_join(., y = dplyr::select(lpd_gbif, Binomial, gbif_species = gbif.species.tree),
            by = "Binomial") %>% 
  mutate(year_s = as.numeric(scale(year)))

glimpse(mammal_weather_5)

# model for each record
mnanom_5km_GAM_5yr <- mammal_weather_5 %>% 
  group_by(ID_block) %>% 
  group_modify(~{
    
    # Temperature
    mod_temp = gamm(pop_growth_rate ~ s(year, bs = "tp", k = 2) + mean_temp_anomaly,
                    data = ., family = gaussian,
                    correlation = corARMA(form = ~ year, p = 1),
                    method = "REML")
    coef_tempmod = coef(mod_temp$gam)
    
    # Precipitation + dealing with NA values
    if(length(which(is.na(.$mean_precip_anomaly) == T)) == 0){
      mod_precip = gamm(pop_growth_rate ~ s(year, bs = "tp", k = 2) + mean_precip_anomaly,
                        data = ., family = gaussian,
                        correlation = corARMA(form = ~ year, p = 1),
                        method = "REML")
      coef_precipmod = coef(mod_precip$gam)}
    else{coef_precipmod = rep(NA,15)}     # Arbitrary long NA vector
    
    tibble(.[1,],
           coef_temp = unname(coef_tempmod[2]),
           coef_precip = unname(coef_precipmod[2]),
           n_obs = nrow(.))
  }) 

glimpse(mnanom_5km_GAM_5yr)

##__________________________________________________________________________________________________
#### 3. Long data ####

# linking to weather data and species names 
mammal_weather_20 <- mammal20 %>% 
  left_join(., y = mam_chelsa_annual, by = c("ID", "year")) %>% 
  left_join(., y = dplyr::select(lpd_gbif, Binomial, gbif_species = gbif.species.tree),
            by = "Binomial") %>% 
  mutate(year_s = as.numeric(scale(year)))

glimpse(mammal_weather_20)

# model for each record
mnanom_5km_GAM_20yr <- mammal_weather_20 %>% 
  group_by(ID_block) %>% 
  group_modify(~{
    
    # Temperature
    mod_temp = gamm(pop_growth_rate ~ s(year, bs = "tp", k = 5) + mean_temp_anomaly,
                    data = ., family = gaussian,
                    correlation = corARMA(form = ~ year, p = 1),
                    method = "REML")
    coef_tempmod = coef(mod_temp$gam)
    
    # Precipitation + dealing with NA values
    if(length(which(is.na(.$mean_precip_anomaly) == T)) == 0){
      mod_precip = gamm(pop_growth_rate ~ s(year, bs = "tp", k = 5) + mean_precip_anomaly,
                        data = ., family = gaussian,
                        correlation = corARMA(form = ~ year, p = 1),
                        method = "REML")
      coef_precipmod = coef(mod_precip$gam)}
    else{coef_precipmod = rep(NA,15)}     # Arbitrary long NA vector
    
    tibble(.[1,],
           coef_temp = unname(coef_tempmod[2]),
           coef_precip = unname(coef_precipmod[2]),
           n_obs = nrow(.))
  }) 

glimpse(mnanom_5km_GAM_20yr)

##__________________________________________________________________________________________________
#### 5. Save data ####

save(mnanom_5km_GAM_5yr, mnanom_5km_GAM_20yr,
     file = "data/pgr_weather/mnanom_5km_GAM_varied_record_length.RData")


