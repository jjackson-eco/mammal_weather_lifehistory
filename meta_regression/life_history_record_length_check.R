####################################################
##                                                ##
##  Global climate and mammal population dynamics ##
##                                                ##
##        Life-history record length link         ##
##                                                ##
##                Jun 6th 2022                    ##
##                                                ##
####################################################

rm(list = ls())
options(width = 100)

#_______________
## Packages
library(tidyverse)
library(patchwork)

#_______________
## Analysis data
load("data/mammal_analysis_data_GAM.RData", verbose = TRUE)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 1. Wrangle and plot ####

lh_n <- mam_coef %>% 
  group_by(species) %>% 
  summarise(longevity = longevity[1], litter = litter[1], 
            bodymass = bodymass[1], mn_n_obs = mean(n_obs),
            n_records = n())

cor.test(lh_n$longevity, lh_n$mn_n_obs)

long_n <- ggplot(lh_n, aes(x = longevity, y = mn_n_obs, size = n_records)) +
  geom_point(alpha =0.7) +
  annotate("text", label = "r = 0.13, p = 0.11", x = -1, y = 25, size = 3) +
  scale_size_continuous(range = c(1,9)) +
  labs(x = "Standardised maximum longevity", y = "Mean length of record (years)",
       size = "Number of\nrecords") +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())

cor.test(lh_n$litter, lh_n$mn_n_obs)

lit_n <- ggplot(lh_n, aes(x = litter, y = mn_n_obs, size = n_records)) +
  geom_point(alpha =0.7) +
  annotate("text", label = "r = -0.10, p = 0.21", x = 2, y = 25, size = 3) +
  scale_size_continuous(range = c(1,9)) +
  labs(x = "Standardised mean litter size", y = "Mean length of record (years)",
       size = "Number of\nrecords") +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())

ggsave(long_n + lit_n, filename = "plots/manuscript_figures/Supplementary figures/lh_nobs.jpeg",
       width = 34, height = 12.5, units = "cm", dpi = 1000)

  
