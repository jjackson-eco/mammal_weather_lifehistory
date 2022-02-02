###################################################################
##                                                               ##
##            Global climate and population dynamics             ##
##                                                               ##
##   Obtaining population growth data for short + long records   ##
##                                                               ##
##                      Feb 2nd 2022                             ##
##                                                               ##
###################################################################
rm(list = ls())
options(width = 100)

library(tidyverse)
library(patchwork)
library(flextable)

## Load in the raw mammal data
load("../rawdata/mam.RData", verbose = T)

mam_raw <- mam # for summary later
summarise(mam_raw, n = n(), n_records = n_distinct(ID))

##__________________________________________________________________________________________________
#### 1. Basic data operations to extract blocks and gaps ####

# Keeping only years from CHELSA and only 5 years of study beyond that.
mam <- mam %>% 
  filter(year >= 1979 & year <= 2013) %>% # important to do this first
  group_by(ID) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  filter(n >= 5) %>% 
  dplyr::select(-n)

# gaps in abundance
mam_gaps <- mam %>% 
  group_by(ID) %>% 
  group_modify(~{
    cyears = .$year
    diff_cyears = diff(cyears)
    cumsum_blocks = cumsum(c(1, diff_cyears != 1))
    
    summarise(., Binomial = Binomial[1],
              record_length = length(cyears),
              no_consecutive_blocks = n_distinct(cumsum_blocks),
              prop_1year_transitions = sum(diff_cyears == 1)/ length(diff_cyears),
              longest_block = max(table(cumsum_blocks)))
  }) %>% 
  ungroup()

# consecutive blocks
mam_blocks <- mam %>% 
  group_by(ID) %>%
  mutate(block = cumsum(c(1, diff(year) != 1)),
         max_block = max(block)) %>% 
  ungroup() %>% 
  dplyr::select(ID, Binomial, Order, ln_abundance, 
                year, block, max_block) %>% 
  left_join(x = ., y = dplyr::select(mam_gaps, -c(Binomial, no_consecutive_blocks)),
            by = "ID") %>% 
  arrange(desc(longest_block)) %>% 
  mutate(ID = factor(ID, levels = unique(.$ID)))

##__________________________________________________________________________________________________
#### 2. Short blocks of 5 years or more ####

# IDs and blocks that we want to keep
ID_block_keep_short <- mam_blocks %>% 
  mutate(ID = as.numeric(as.character(ID))) %>% 
  group_by(ID, block) %>% 
  summarise(ID_block = paste0(ID[1],"_",block[1]),
            block_keep = if_else(n() >= 5, 1, 0)) %>% ## <- specify here that we want consecutive blocks of 5 years or more
  ungroup() %>% 
  filter(block_keep == 1)

# Restricting the dataset
mam_IDblocks_5yr <- mam %>% 
  group_by(ID) %>%
  mutate(block = cumsum(c(1, diff(year) != 1)),
         ID_block = paste0(ID[1],"_",block)) %>% 
  ungroup() %>% 
  filter(ID_block %in% ID_block_keep_short$ID_block == T) %>% 
  dplyr::select(1,21,22,2:9,11:20)

# Summarising this in a nice table
mam_datasum5 <- data.frame(Dataset = c("Raw data", "Study data >= 5 years"),
                          Observations  = c(nrow(mam_raw), 
                                            nrow(mam_IDblocks_5yr)),
                          Records = c(n_distinct(mam_raw$ID), 
                                      n_distinct(mam_IDblocks_5yr$ID)),
                          Species = c(n_distinct(mam_raw$Binomial),
                                      n_distinct(mam_IDblocks_5yr$Binomial)))
flextable(mam_datasum5, cwidth = 2)

##__________________________________________________________________________________________________
#### 3. Long blocks of 20 years or more ####

# IDs and blocks that we want to keep
ID_block_keep_long <- mam_blocks %>% 
  mutate(ID = as.numeric(as.character(ID))) %>% 
  group_by(ID, block) %>% 
  summarise(ID_block = paste0(ID[1],"_",block[1]),
            block_keep = if_else(n() >= 20, 1, 0)) %>% ## <- specify here that we want consecutive blocks of 20 years or more
  ungroup() %>% 
  filter(block_keep == 1)

# Restricting the dataset
mam_IDblocks_20yr <- mam %>% 
  group_by(ID) %>%
  mutate(block = cumsum(c(1, diff(year) != 1)),
         ID_block = paste0(ID[1],"_",block)) %>% 
  ungroup() %>% 
  filter(ID_block %in% ID_block_keep_long$ID_block == T) %>% 
  dplyr::select(1,21,22,2:9,11:20)

# Summarising this in a nice table
mam_datasum20 <- data.frame(Dataset = c("Raw data", "Study data >= 20 years"),
                          Observations  = c(nrow(mam_raw), 
                                            nrow(mam_IDblocks_20yr)),
                          Records = c(n_distinct(mam_raw$ID), 
                                      n_distinct(mam_IDblocks_20yr$ID)),
                          Species = c(n_distinct(mam_raw$Binomial),
                                      n_distinct(mam_IDblocks_20yr$Binomial)))
flextable(mam_datasum20, cwidth = 2)

##__________________________________________________________________________________________________
#### 4. Population growth rates for each study ####

# Not restricting based on the 0 observations for this round - this is to explore the implications of the 10 year choice
# This will remove one year from each ID_block

# short
mammal5 <- mam_IDblocks_5yr %>% 
  mutate(raw_abundance2 = raw_abundance + 1) %>% # For the 0 observations, adjusting all to non-zeros
  group_by(ID_block) %>% 
  group_modify(~{
    t0 <- .$raw_abundance2[-(length(.$raw_abundance2))] # get rid the last obs
    t1 <- .$raw_abundance2[-1]                        # get rid of the first obs
    
    mutate(., pop_growth_rate = c(log(t1/t0),NA))
  }) %>% 
  ungroup() %>% 
  filter(is.na(pop_growth_rate) == F)

# long
mammal20 <- mam_IDblocks_20yr %>% 
  mutate(raw_abundance2 = raw_abundance + 1) %>% # For the 0 observations, adjusting all to non-zeros
  group_by(ID_block) %>% 
  group_modify(~{
    t0 <- .$raw_abundance2[-(length(.$raw_abundance2))] # get rid the last obs
    t1 <- .$raw_abundance2[-1]                        # get rid of the first obs
    
    mutate(., pop_growth_rate = c(log(t1/t0),NA))
  }) %>% 
  ungroup() %>% 
  filter(is.na(pop_growth_rate) == F)

##__________________________________________________________________________________________________
#### 5. Save ####

save(mammal5, mammal20, file = "../rawdata/mammal_varied_record_length.RData")



