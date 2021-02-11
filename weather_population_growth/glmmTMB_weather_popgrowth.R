####################################################
##                                                ##
##     Global climate and population dynamics     ##
##                                                ##
##          glmmTMB timeseries method             ##
##                                                ##
##                Dec 11th 2020                   ##
##                                                ##
####################################################

# Record-wise regressions linking weather to population growth rates, 
# accounting for autocorrelation with glmmTMB and a formal AR(1) time-series correlation structure

rm(list = ls())
options(width = 100)

library(tidyverse)
library(viridis)
library(patchwork)
library(glmmTMB)  
library(mgcv)

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
#### 3. White tailed-deer example ####

wtd <- filter(mammal_weather, ID_block == "22039_1")

#____________________________________
## Temperature

wtd_glmmTMB <- glmmTMB(pop_growth_rate ~ mean_temp_anomaly + ar1(as.factor(year_f) + 0 | ID_block),
                       family = 'gaussian', data = wtd)

wtd_gamarma <- gamm(pop_growth_rate ~ s(year, bs = "tp") + mean_temp_anomaly,
                    data = wtd, family = gaussian,
                    correlation = corARMA(form = ~ year, p = 1),
                    method = "REML")

wtd_linear <- lm(pop_growth_rate ~ mean_temp_anomaly + ln_abundance + year, data = wtd)
wtd_naive <- lm(pop_growth_rate ~ mean_temp_anomaly, data = wtd)

## The GAMM looks to be doing a similar job, whilst formally accounting for temporal autocorrelation


