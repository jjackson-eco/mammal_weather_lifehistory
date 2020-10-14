####################################################
##                                                ##
##     Global climate and population dynamics     ##
##                                                ##
##    Simple brms Phylogenetic Meta-regression    ##
##                                                ##
##                Oct 14th 2020                   ##
##                                                ##
####################################################

## Simplest single observation for each species for mean temp/precip anomalies at 5km resolution
## Testing brms

rm(list = ls())
options(width = 100)

library(tidyverse)
library(ape)
library(brms)

##__________________________________________________________________________________________________
#### 1. Load data ####

## Weather effects
load("data/pgr_weather/mnanom_5km.RData", verbose = TRUE)

## Phylogenetic
load("../Phlogenetics/mam_lpd_pruned.RData", verbose = TRUE)

##__________________________________________________________________________________________________
#### 2. Mammal coefficient data ####

mam_coef <- mnanom_5km %>% 
  ungroup() %>% 
  dplyr::select(ID, Order, coef_temp, coef_precip) %>% 
  left_join(x = ., 
            y = dplyr::select(LPD_tree_update, ID, checked_speciesname),
            by = "ID") %>% 
  mutate(checked_speciesname = gsub(" ", "_", checked_speciesname)) %>% 
  # Average coefficient for each species
  group_by(checked_speciesname) %>% 
  summarise(Order = Order[1],
            coef_temp = mean(coef_temp), 
            coef_precip = mean(coef_precip),
            .groups = "drop") %>% 
  dplyr::rename(spp = checked_speciesname)

##__________________________________________________________________________________________________
#### 3. Phylogenetic covariance matrix ####

# Trim tree to right data
mamMCC_weatheff <- keep.tip(mamMCC_pruned, mam_coef$spp)

# Covariance matrix - Brownian motion model
A <- ape::vcv.phylo(mamMCC_weatheff)

##__________________________________________________________________________________________________
#### 4. Sampling from the posterior ####

temp_simple <- brm(
  coef_temp ~ 1 + (1|gr(spp, cov = A)), 
  data = mam_coef, family = gaussian(), 
  data2 = list(A = A),
  prior = c(
    prior(normal(0, 0.5), "Intercept"),
    prior(student_t(3, 0, 20), "sd")),
  control = list(adapt_delta = 0.95),
  chains = 4, cores = 4, iter = 2000, warmup = 500
)

precip_simple <- brm(
  coef_precip ~ 1 + (1|gr(spp, cov = A)), 
  data = mam_coef, family = gaussian(), 
  data2 = list(A = A),
  prior = c(
    prior(normal(0, 0.5), "Intercept"),
    prior(student_t(3, 0, 20), "sd")),
  control = list(adapt_delta = 0.95),
  chains = 4, cores = 4, iter = 2000, warmup = 500
)


