#####################################################
##                                                 ##
##     Global climate and population dynamics      ##
##                                                 ##
##   Annual Abundance Timeseries Gap Exploration   ##
##                                                 ##
##               April 23rd 2020                   ##
##                                                 ##
#####################################################
rm(list = ls())
options(width = 100)

library(tidyverse)
library(gridExtra)
library(grid)

## Load in the raw mammal data
load("../rawdata/mam.RData", verbose = T)

##__________________________________________________________________________________________________
#### 1. Set up data ####

# Keeping only years from CHELSA and only 5 years of study beyond that.
mam <- mam %>% 
  filter(year >= 1979 & year <= 2013) %>% # important to do this first
  group_by(ID) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  filter(n >= 5) %>% 
  dplyr::select(-n)

##__________________________________________________________________________________________________
#### 2. Time series gaps for each study ####

dupyer <- mam %>% 
  group_by(ID) %>% 
  summarise(Binomial = Binomial[1],
            duplicate_years = n() -n_distinct(year))

mam %>% 
  group_by(ID) %>% 
  group_modify(~{
    cyears = .$year
  })






