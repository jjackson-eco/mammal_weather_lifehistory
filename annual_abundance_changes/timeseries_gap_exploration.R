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

# Check for duplicated years
# dupyer <- mam %>% 
#   group_by(ID) %>% 
#   summarise(Binomial = Binomial[1],
#             duplicate_years = n() -n_distinct(year))

##__________________________________________________________________________________________________
#### 2. Testing - Fossa  ####

test <- filter(mam, ID == 25396)

ggplot(test, aes(x = year, y = scaled_abundance)) + 
  geom_point()

tyear <- test$year
diff_tyear <- diff(tyear)
cumsum_blocks <- cumsum(c(1, diff_tyear != 1))

no_consecutive_blocks <- n_distinct(cumsum_blocks)
prop_1year_transitions <- sum(diff_tyear == 1)/ length(diff_tyear)
longest_block <- max(table(cumsum_blocks))

# # interesting function
# consecutive_blocks <-  length(split(tyear, cumsum(c(1, diff(tyear) != 1))))

##__________________________________________________________________________________________________
#### 3. Exploring gaps for all studies ####

mam_gaps <- mam %>% 
  group_by(ID) %>% 
  group_modify(~{
    cyears = .$year
    diff_cyears = diff(cyears)
    cumsum_blocks = cumsum(c(1, diff_cyears != 1))
    
    summarise(., Binomial = Binomial[1],
              study_length = length(cyears),
              no_consecutive_blocks = n_distinct(cumsum_blocks),
              prop_1year_transitions = sum(diff_cyears == 1)/ length(diff_cyears),
              longest_block = max(table(cumsum_blocks)))
  })

##__________________________________________________________________________________________________
#### 4. Exploration plots ####

ggplot(mam_gaps, aes(x = prop_1year_transitions)) +
  geom_histogram(bins = 10, fill = "lightblue", 
                 colour = "black", size = 0.1) +
  labs(x = "Proportion of study in 1-year transitions",
       y = "Number of studies") +
  theme_bw(base_size = 11) +
  ggsave("plots/annual_abundance/Proportion_1year_transitions.jpeg",
         width = 4, height = 3, units = "in", dpi = 400)




