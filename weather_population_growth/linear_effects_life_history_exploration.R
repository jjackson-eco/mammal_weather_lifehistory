####################################################
##                                                ##
##     Global climate and population dynamics     ##
##                                                ##
##     Life-history weather effect exploration    ##
##                                                ##
##                Oct 28th 2020                   ##
##                                                ##
####################################################

## Exploration of potential links between weather effects and life-history

rm(list = ls())
options(width = 100)

library(tidyverse)

##__________________________________________________________________________________________________
#### 1. Load data ####

# mammal coefficient data - Mean anomaly in a 5km buffer radius
load("data/pgr_weather/mnanom_5km.RData")
glimpse(mnanom_5km)

# life-history data
load("../rawdata/mam_dski.RData")
glimpse(mam_dski)

# species names 
load("../../LPI_rawdata/GBIF_species_names.RData", verbose = T)
rm(dski_gbif, mamMCC_gbif)

##__________________________________________________________________________________________________
#### 2. Merging data ####

## 2a. Species names for mammal data
mamcoef <- mnanom_5km %>% 
  left_join(x = ., y = lpd_gbif, by = "Binomial")

# any not in gbif? NOPE :)
length(which(is.na(mamcoef$gbif.species.id)))














