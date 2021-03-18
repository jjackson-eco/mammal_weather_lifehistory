####################################################
##                                                ##
##     Global climate and population dynamics     ##
##                                                ##
##          Prior Predictive Simulation           ##
##                                                ##
##                Mar 18th 2021                   ##
##                                                ##
####################################################

## Simple prior simulations for the parameters of our meta-regression on the
## terrestrial mammals

rm(list = ls())
options(width = 100)

library(tidyverse)
library(ape)
library(brms)
library(ggridges)

# colours
temp_colour <- "#990a80"
precip_colour <- "#287f79"

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 1. Load data ####

load("data/mammal_analysis_data_GAM.RData")

glimpse(mam_coef)
glimpse(mamMCC_coef)
glimpse(A_coef)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 2. Intercept PPS ####

set.seed(1000)
tibble(`Normal(0, 10)\nWeak` = rnorm(n = 500, mean = 0, sd = 10),
       `Normal(0, 2)\nMid` = rnorm(n = 500, mean = 0, sd = 2),
       `Normal(0, 0.5)\nRegularising` = rnorm(n = 500, mean = 0, sd = 0.5),
       row = 1:500) %>% 
  pivot_longer(-row) %>% 
  bind_rows(tibble(row = rep(1:nrow(mam_coef), 2), 
                   name = rep(c("Raw\ntemperature\ncoefficient\ndata",
                                "Raw\nprecipitation\ncoefficient\ndata"), 
                              each = nrow(mam_coef)), 
                   value = c(mam_coef$coef_temp, mam_coef$coef_precip))) %>% 
  mutate(name = factor(name, levels = c("Raw\ntemperature\ncoefficient\ndata",
                                        "Raw\nprecipitation\ncoefficient\ndata",
                                        "Normal(0, 10)\nWeak", 
                                        "Normal(0, 2)\nMid",
                                        "Normal(0, 0.5)\nRegularising"))) %>% 
  ggplot(aes(x = value, fill = name, y = name)) +
  geom_density_ridges(show.legend = F) +
  labs(x = expression(paste("Weather coefficient - Global intercept value ", bar(alpha))),
       y = NULL) +
  coord_cartesian(xlim = c(-10,10)) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank()) +
  ggsave("plots/meta_regression/prior_predictive_simulation/weather_coefficient_pps.jpeg",
         width = 15, height = 13, units = "cm", dpi = 500)
  
##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3. Biome effects ####












