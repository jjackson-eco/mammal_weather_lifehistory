#####################################################
##                                                 ##
##     Global climate and population dynamics      ##
##                                                 ##
##         Annual population growth rates          ##
##                                                 ##
##                 June 8th 2020                   ##
##                                                 ##
#####################################################
rm(list = ls())
options(width = 100)

library(tidyverse)
library(viridis)
library(grid)
library(gridExtra)
library(psych)

##__________________________________________________________________________________________________
#### 1. Load data ####

load("../rawdata/mam_IDblocks.RData")
glimpse(mam_IDblocks)

##__________________________________________________________________________________________________
#### 2. Calculate per-capita residual growth rate between time t and t+1 ####

# This will remove one year from each ID_block - Does this matter?
mammal <- mam_IDblocks %>% 
  group_by(ID_block) %>% 
  group_modify(~{
    t0 <- .$ln_abundance[-(length(.$ln_abundance))]
    t1 <- .$ln_abundance[-1]
    
    mutate(., pop_growth_rate = c(t1/t0,NA))
  }) %>% 
  ungroup() %>% 
  filter(is.na(pop_growth_rate) == F)

##__________________________________________________________________________________________________
#### 3. Exploratory plots ####

# 3a. Histogram
ggplot(mammal, aes(x = pop_growth_rate)) +
  geom_histogram(bins = 100, fill = "black") +
  labs(x = "Per-capita population growth rate", 
       y = "Frequency") +
  theme_bw(base_size = 12) +
  ggsave(filename = "plots/annual_abundance/pop_growth_rate_histogram.jpeg",
         width = 6, height = 4, units = "in", dpi = 400)

# 3b. Density dependence
ggplot(mammal, aes(x = ln_abundance, y = pop_growth_rate)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_point(alpha = 0.15, colour = viridis(10)[6], size = 2.5) +
  geom_smooth(method = "lm", se = F, colour = "black") +
  coord_cartesian(ylim = c(0,3)) +
  labs(x = "ln Abundance", y = "Per-capita population growth rate") +
  theme_bw(base_size = 18) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggsave(filename = "plots/annual_abundance/density_dependence_mam.jpeg",
         width = 7, height = 7, units = "in", dpi = 400)


