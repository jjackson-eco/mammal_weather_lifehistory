######################################################
##                                                  ##
##     Global climate and population dynamics       ##
##                                                  ##
## Linear Detrending of Annual Abundance Timeseries ##
##                                                  ##
##                May 12th 2020                     ##
##                                                  ##
######################################################
rm(list = ls())
options(width = 100)

library(tidyverse)
library(gridExtra)
library(grid)

##__________________________________________________________________________________________________
#### 1. Load and set up data ####

# Loading the raw mammal data with the blocks
load("../rawdata/mammal_data.RData")
glimpse(mammal)

# Checking number of blocks and records
mammal %>% 
  group_by(ID) %>% 
  summarise(n_block = n_distinct(block), n_ID = 1) %>% 
  ungroup() %>% 
  mutate(num = n_ID*n_block) %>% 
  summarise(sum(num))
  
n_distinct(mammal$ID_block) # Seem to match up ok

##__________________________________________________________________________________________________
#### 2. Linear detrend for each block of each study ####

# extracting the residuals after fitting a linear trend to each timeseries
mam_detrend <- mammal %>% 
  group_by(ID_block) %>% 
  group_modify(~{
    mod = lm(scaled_abundance ~ year, data = .) 
    
    mutate(., residual_abundance = mod$residuals,
           coef = mod$coefficients[2])
  }) %>% 
  ungroup()















