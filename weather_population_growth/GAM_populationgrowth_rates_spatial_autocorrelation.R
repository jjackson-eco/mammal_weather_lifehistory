####################################################
##                                                ##
##     Global climate and population dynamics     ##
##                                                ##
##     GAMM Spatial Autocorrelation model full    ##
##                                                ##
##                Feb 12th 2021                   ##
##                                                ##
####################################################

# One large GAM linking weather to population growth rates, split for each record. Also
# accounting for autocorrelation with GAMM with a formal AR(1) time-series mixed effect.

rm(list = ls())
options(width = 100)

library(tidyverse)
library(viridis)
library(patchwork)
library(mgcv)      # Simon Wood to the rescue again. All Hail

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
  mutate(year_s = as.numeric(scale(year)))

glimpse(mammal_weather)

##__________________________________________________________________________________________________
#### 3. Modelling population change with respect to weather ####

# taking a random sample of data as a test
set.seed(100)
mamdat <- mammal_weather %>% 
  filter(ID_block %in% )

mod_temp <- gam(pop_growth_rate ~ s(year, by = ID_block, bs = "tp", k = 5) +
                mean_temp_anomaly:ID_block, data = mammal_weather, 
                family = gaussian, method = "REML")


## simulate example as a test
set.seed(10)
simdat <- expand_grid(pop = as.factor(LETTERS[1:10]), year = 1:10) %>% 
  mutate(pgr = rnorm(n(), 0, 0.5),
         weather = rnorm(n(), 0, 0.5))


ggplot(simdat, aes(x = year, y = pgr, colour = pop)) +
  geom_line()

simmod <- gam(pgr ~ s(year, by = pop) + weather,
              data = simdat, 
              family = gaussian, method = "REML")


