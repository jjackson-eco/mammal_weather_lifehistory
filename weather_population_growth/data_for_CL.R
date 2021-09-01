####################################################
##                                                ##
##     Global climate and population dynamics     ##
##                                                ##
##       Weather Popgrowth timeseries data        ##
##                                                ##
##                Aug 26th 2021                   ##
##                                                ##
####################################################

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
#### 3. Data for Christie ####

# data selection
mw <- mammal_weather %>% 
  group_by(ID_block) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  dplyr::rename(raw_abundance_no0 = raw_abundance2) %>% 
  dplyr::select(ID, Family, Binomial, gbif_species, Biome = biome, Latitude, Longitude,
                year, year_s, n, mean_temp_anomaly, mean_precip_anomaly, raw_abundance, raw_abundance_no0,
                pop_growth_rate) %>% 
  arrange(n, ID)

# selecting long records and short records
uID <- unique(mw$ID)

mammal_weather_CL <- mw %>% 
  filter(ID %in% uID[1:5] | ID %in% uID[(length(uID) - 4): length(uID)])

save(mammal_weather_CL, file = "../rawdata/mammal_weather_CL.RData")


