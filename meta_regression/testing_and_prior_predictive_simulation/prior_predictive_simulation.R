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

# simulate intercepts for the coefficients in a tible
tibble(`Normal(0, 10)\nWeak` = rnorm(n = 500, mean = 0, sd = 10),
       `Normal(0, 2)\nMid` = rnorm(n = 500, mean = 0, sd = 2),
       `Normal(0, 0.5)\nRegularising` = rnorm(n = 500, mean = 0, sd = 0.5),
       row = 1:500) %>% 
  # pivot longer and put them with raw data
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
  geom_density_ridges(show.legend = F, size = 0.2) +
  labs(x = expression(paste("Weather coefficient - Global intercept value ", bar(alpha))),
       y = NULL) +
  coord_cartesian(xlim = c(-10,10)) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank()) +
  ggsave("plots/meta_regression/prior_predictive_simulation/weather_coefficient_pps.jpeg",
         width = 15, height = 13, units = "cm", dpi = 500)
  
##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3. Beta difference effects ####

# Raw data - pairwise estimates of records as a baseline
raw_dat_diffs <- tibble(name = "Raw\npairwise\ncoefficient\ndifferences",
                        value = apply(combn(mam_coef$coef_temp,2), 2, diff)) 

set.seed(1000)
# Prior simulation - normal beta priors
tibble(`Normal(0, 10)\nWeak` = rnorm(n = 500, mean = 0, sd = 10),
       `Normal(0, 2)\nMid` = rnorm(n = 500, mean = 0, sd = 2),
       `Normal(0, 0.5)\nRegularising` = rnorm(n = 500, mean = 0, sd = 0.5)) %>% 
  pivot_longer(cols = 1:3) %>% 
  bind_rows(raw_dat_diffs) %>% 
    mutate(name = factor(name, levels = c("Raw\npairwise\ncoefficient\ndifferences",
                                          "Normal(0, 10)\nWeak", 
                                          "Normal(0, 2)\nMid",
                                          "Normal(0, 0.5)\nRegularising"))) %>% 
  ggplot(aes(x = name, y = value, fill = name)) +
  geom_violin(show.legend = F, size = 0.2) +
  coord_cartesian(ylim = c(-10,10)) +
  labs(x = NULL, y = "Differences in weather coefficients") +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank()) +
  ggsave("plots/meta_regression/prior_predictive_simulation/coefficient_differences_beta.jpeg",
         width = 15, height = 13, units = "cm", dpi = 500)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3. Simulated life-history effects ####

# Simple linear lognormal model
set.seed(1000)
tibble(sim = rep(1:100, each = 50),
       life_history = rep(seq(-2,2,length.out = 50), times = 100),
       # simulate life-history slopes on linear predictor scale using normal
       beta_weak = rep(rnorm(100, 0, 10), each = 50),
       beta_mid = rep(rnorm(100, 0, 2), each = 50),
       beta_reg = rep(rnorm(100, 0, 0.5), each = 50)) %>% 
  # predictions with a back-transformation
  mutate(`Normal(0, 10)\nWeak` = exp(beta_weak*life_history),
         `Normal(0, 2)\nMid` = exp(beta_mid*life_history),
         `Normal(0, 0.5)\nRegularising` = exp(beta_reg*life_history)) %>%
  dplyr::select(-c(beta_weak, beta_mid, beta_reg)) %>% 
  pivot_longer(-c(life_history, sim)) %>% 
  filter(value < 100) %>% 
  mutate(name = factor(name, levels = c("Normal(0, 10)\nWeak", 
                                        "Normal(0, 2)\nMid",
                                        "Normal(0, 0.5)\nRegularising"))) %>% 
  ggplot(aes(x = life_history, y = value, group = sim, colour = name)) +
  geom_line(alpha = 0.5, show.legend = F) +
  geom_hline(yintercept = max(mam_coef$abs_temp)) +
  geom_hline(yintercept = max(mam_coef$abs_precip, na.rm = TRUE), linetype = "dashed") +
  labs(x = "Standardised life-history value", y = "|Weather coefficient|") +
  facet_wrap(~name, scales = "free") +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank()) +
  ggsave(filename = "plots/meta_regression/prior_predictive_simulation/life_history_effect_pps.jpeg",
         width = 26, height = 13, units = "cm", dpi = 500)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 4. Error terms - Exponential priors ####

# Simulate exponentials with different rates
rates <- c(0.5,1,2,5,10,20)
exp_dat <- tibble(rate = rep(rates, each = 500),
                  value = 0)
for(i in rates){
  exp_dat[exp_dat$rate == i,]$value <- rexp(500, rate = i)
}

ggplot(exp_dat, aes(x = value, y = as.factor(rate),
                    fill = as.factor(rate))) +
  geom_density_ridges(show.legend = F, rel_min_height=0.005, size = 0.2) +
  geom_vline(xintercept = 1) +
  labs(x = "Random effect variance term value", y = "Exponential rate") +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank()) +
  ggsave("plots/meta_regression/prior_predictive_simulation/random_effect_variance_pps.jpeg",
         width = 17, height = 9, units = "cm", dpi = 500)




