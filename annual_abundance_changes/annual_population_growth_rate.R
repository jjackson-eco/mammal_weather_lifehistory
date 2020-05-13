########################################################
##                                                    ##
##      Global climate and population dynamics        ##
##                                                    ##
## Annual per-capita residual population growth rates ##
##                                                    ##
##                 May 13th 2020                      ##
##                                                    ##
########################################################
rm(list = ls())
options(width = 100)

library(tidyverse)
library(gridExtra)
library(grid)

##__________________________________________________________________________________________________
#### 1. Load data ####

load("data/mam_detrend.RData")
glimpse(mam_detrend)

##__________________________________________________________________________________________________
#### 2. Calculate per-capita residual growth rate between time t and t+1 ####

# This will remove one year from each ID_block - Does this matter?










