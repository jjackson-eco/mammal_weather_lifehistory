####################################################
##                                                ##
##     Global climate and population dynamics     ##
##                                                ##
##    GAMM timeseries method - Weather variance   ##
##                                                ##
##                Aug 9th 2021                    ##
##                                                ##
####################################################

# Record-wise regressions linking the variance in weather to population growth rates, 
# accounting for autocorrelation with GAMM and formal AR(1) time-series analysis

rm(list = ls())
options(width = 100)

library(tidyverse)
library(viridis)
library(patchwork)
library(mgcv) 

##__________________________________________________________________________________________________
#### 1. Load data ####

# mammal data
load("../rawdata/mammal.RData")
glimpse(mammal)

# annual weather anomaly - focus on just the mean anomaly in this script at a 5km range
mam_chelsa_annual <- readRDS("data/mam_chelsa_annual.RDS") %>% 
  filter(scale == "scale_5km") %>% 
  dplyr::select(ID,year, weather_scale = scale, temp_variance, precip_variance) %>% 
  # Scaling variance values
  mutate(temp_variance = as.numeric(scale(temp_variance)),
         precip_variance = as.numeric(scale(precip_variance)))
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

pgr_weathervar_gam <- mammal_weather %>% 
  group_by(ID_block) %>% 
  group_modify(~{
    
    # Temperature
    mod_temp = gamm(pop_growth_rate ~ s(year, bs = "tp", k = 5) + temp_variance,
                    data = ., family = gaussian,
                    correlation = corARMA(form = ~ year, p = 1),
                    method = "REML")
    coef_tempmod = coef(mod_temp$gam)
    
    # Precipitation + dealing with NA values
    if(length(which(is.na(.$precip_variance) == T)) == 0){
      mod_precip = gamm(pop_growth_rate ~ s(year, bs = "tp", k = 5) + precip_variance,
                        data = ., family = gaussian,
                        correlation = corARMA(form = ~ year, p = 1),
                        method = "REML")
      coef_precipmod = coef(mod_precip$gam)}
    else{coef_precipmod = rep(NA,15)}     # Arbitrary long NA vector
    
    tibble(.[1,],
           coef_tempvar = unname(coef_tempmod[2]),
           coef_precipvar = unname(coef_precipmod[2]),
           n_obs = nrow(.))
  }) 

glimpse(pgr_weathervar_gam)

hist(pgr_weathervar_gam$coef_tempvar)
hist(pgr_weathervar_gam$coef_precipvar)

##__________________________________________________________________________________________________
#### 4. Saving data from GAM ####

var_5km_GAM <- pgr_weathervar_gam
save(var_5km_GAM, file = "data/pgr_weather/var_5km_GAM.RData")







