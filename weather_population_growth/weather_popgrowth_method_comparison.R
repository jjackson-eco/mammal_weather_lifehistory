####################################################
##                                                ##
##     Global climate and population dynamics     ##
##                                                ##
##  Comparing methods of fitting weather effects  ##
##                                                ##
##                Dec 11th 2020                   ##
##                                                ##
####################################################

# Record-wise regressions linking weather to population growth rates, 
# comparing 5 models that account (or not) for temporal autocorrelation 

rm(list = ls())
options(width = 100)

library(tidyverse)
library(viridis)
library(patchwork)
library(glmmTMB)  
library(mgcv)
library(psych)

temp_colour <- viridis(20, option = "C")[13]
precip_colour <- viridis(20, option = "D")[10]

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
  mutate(year_s = as.numeric(scale(year)),
         year_f = as.factor(year))

glimpse(mammal_weather)

##__________________________________________________________________________________________________
#### 3. Models for each record ####

pgr_weather <- mammal_weather %>% 
  group_by(ID_block) %>% 
  group_modify(~{
    
    #_______________________________________________
    # Temperature
    
    # Model 1. Fully naive model - simple linear regression
    temp_naive = lm(pop_growth_rate ~ mean_temp_anomaly, data = .)
    coef_temp_naive = coef(temp_naive)[2]
    
    # Model 2. Linear regression accounting for trend 
    temp_lintr = lm(pop_growth_rate ~ mean_temp_anomaly + year, data = .)
    coef_temp_lintr = coef(temp_lintr)[2]
    
    # Model 3. Linear regression with trend and past abundance
    temp_linear = lm(pop_growth_rate ~ mean_temp_anomaly + ln_abundance + year, data = .)
    coef_temp_linear = coef(temp_linear)[2]
    
    # Model 4. glmmTMB with autoregression(1) for year
    temp_TMB = glmmTMB(pop_growth_rate ~ mean_temp_anomaly + ar1(as.factor(year_f) + 0 | ID), 
                       family = 'gaussian', data = .)
    coef_temp_TMB = as.numeric(coef(temp_TMB)$cond$ID[length(coef(temp_TMB)$cond$ID)])
    
    # Model 5. GAMM with a coarse year smoothing term and an autoregression (1) correlation structure 
    temp_gamm = gamm(pop_growth_rate ~ s(year, bs = "tp", k = 3) + mean_temp_anomaly,
                    data = ., family = gaussian,
                    correlation = corARMA(form = ~ year, p = 1),
                    method = "REML")
    coef_temp_gamm = coef(temp_gamm$gam)[2]
    
    #_______________________________________________
    # Precipitation
    if(length(which(is.na(.$mean_precip_anomaly) == T)) == 0){
      
      # Model 1. Fully naive model - simple linear regression
      precip_naive = lm(pop_growth_rate ~ mean_precip_anomaly, data = .)
      coef_precip_naive = coef(precip_naive)[2]
      
      # Model 2. Linear regression accounting for trend 
      precip_lintr = lm(pop_growth_rate ~ mean_precip_anomaly + year, data = .)
      coef_precip_lintr = coef(precip_lintr)[2]
      
      # Model 3. Linear regression with trend and past abundance
      precip_linear = lm(pop_growth_rate ~ mean_precip_anomaly + ln_abundance + year, data = .)
      coef_precip_linear = coef(precip_linear)[2]
      
      # Model 4. glmmTMB with autoregression(1) for year
      precip_TMB = glmmTMB(pop_growth_rate ~ mean_precip_anomaly + ar1(as.factor(year_f) + 0 | ID), 
                         family = 'gaussian', data = .)
      coef_precip_TMB = as.numeric(coef(precip_TMB)$cond$ID[length(coef(precip_TMB)$cond$ID)])
      
      # Model 5. GAMM with a coarse year smoothing term and an autoregression (1) correlation structure 
      precip_gamm = gamm(pop_growth_rate ~ s(year, bs = "tp", k = 3) + mean_precip_anomaly,
                       data = ., family = gaussian,
                       correlation = corARMA(form = ~ year, p = 1),
                       method = "REML")
      coef_precip_gamm = coef(precip_gamm$gam)[2]}
    else{coef_precip_naive = NA
         coef_precip_lintr = NA
         coef_precip_linear = NA
         coef_precip_TMB = NA
         coef_precip_gamm = NA}    
    
    tibble(.[1,c(1,4,6,13)],
           n_obs = nrow(.),
           coef_temp_naive = coef_temp_naive,
           coef_temp_lintr = coef_temp_lintr,
           coef_temp_linear = coef_temp_linear,
           coef_temp_TMB = coef_temp_TMB,
           coef_temp_gamm = coef_temp_gamm,
           coef_precip_naive = coef_precip_naive,
           coef_precip_lintr = coef_precip_lintr,
           coef_precip_linear = coef_precip_linear,
           coef_precip_TMB = coef_precip_TMB,
           coef_precip_gamm = coef_precip_gamm)
  }) 


##__________________________________________________________________________________________________
#### 4. Comparing the coefficients ####

## Temperature
jpeg("plots/weather_pop_growth/temp_effect_comparison.jpeg",
     width = 20, height = 20, units = "cm",res = 400)
pairs.panels(as.data.frame(pgr_weather %>% ungroup() %>% dplyr::select(contains("temp"))),
             smooth = FALSE, lm = TRUE,  ellipses = FALSE)
dev.off()

## Precipitation
jpeg("plots/weather_pop_growth/precip_effect_comparison.jpeg",
     width = 20, height = 20, units = "cm",res = 400)
pairs.panels(as.data.frame(pgr_weather %>% ungroup() %>% dplyr::select(contains("precip"))),
             smooth = FALSE, lm = TRUE,  ellipses = FALSE)
dev.off()



