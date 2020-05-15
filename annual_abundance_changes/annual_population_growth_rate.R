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
mammal <- mam_detrend %>% 
  group_by(ID_block) %>% 
  group_modify(~{
    resid_t0 <- .$residual_abundance[-(length(.$residual_abundance))]
    resid_t1 <- .$residual_abundance[-1]
    
    mutate(., pop_growth_rate = c(resid_t1/resid_t0,NA))
  }) %>% 
  ungroup() %>% 
  filter(is.na(pop_growth_rate) == F)


# Seems to be some huge values
ggplot(mammal, aes(x = pop_growth_rate)) +
  geom_histogram(bins = 15)

filter(mammal, pop_growth_rate >100)

filter(mammal, ID_block == "10372_1") %>% 
  dplyr::select(ID_block, year, residual_abundance, pop_growth_rate)





