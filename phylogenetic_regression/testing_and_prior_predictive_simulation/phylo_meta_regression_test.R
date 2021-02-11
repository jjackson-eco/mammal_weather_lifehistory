####################################################
##                                                ##
##     Global climate and population dynamics     ##
##                                                ##
##       brms Phylogenetic meta-regression        ##
##                                                ##
##                 Nov 9th 2020                   ##
##                                                ##
####################################################

## Investigating the base model with a meta-regression framework
## Intercept and phylogenetic regression

rm(list = ls())
options(width = 100)

library(tidyverse)
library(ape)
library(brms)
library(nlme)
library(patchwork)
library(ggridges)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 1. Load data ####

## Weather effects
load("data/pgr_weather/mnanom_5km.RData", verbose = TRUE)
glimpse(mnanom_5km)

## Phylogenetic
load("../Phlogenetics/mam_lpd_pruned.RData", verbose = TRUE)
str(mamMCC_pruned)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 2. Mammal coefficient data ####

mam_temp <- mnanom_5km %>% 
  dplyr::select(ID, ID_block, species = gbif_species, 
                Order, biome, Latitude, coef_temp, 
                coef_precip, n_obs) %>% 
  mutate(species = gsub(" ", "_", species),
         phylo = species,
         ## z transform the coefficients
         coef_temp = as.numeric(scale(coef_temp)),
         coef_precip = as.numeric(scale(coef_precip)),
         #z transform absolute latitude for models
         lat = as.numeric(scale(abs(Latitude))),
         OLRE = 1:n()) 

## restrict for precipitation - dropping NAs
mam_precip <- drop_na(mam_temp, coef_precip)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3. Phylogenetic covariance matrix ####

# Trim tree to right data
mamMCC_temp <- keep.tip(mamMCC_pruned, mam_temp$phylo)
mamMCC_precip <- keep.tip(mamMCC_pruned, mam_precip$phylo)

# Covariance matrix - Brownian motion model
A_temp <- ape::vcv.phylo(mamMCC_temp)
A_precip <- ape::vcv.phylo(mamMCC_precip)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 4. Single chain tests ####

# Simple phylogenetic model

# Tested with OLRE, but actually we have meaningful group level variance with species repeats.
# Lots of species have only one observation, so having both of these terms is essentially redundant

set.seed(666)
temp_test <- brm(
  coef_temp | se(sqrt(1 / (n_obs - 3))) ~ 1 + (1|gr(phylo, cov = A_temp)) + (1| species),  
  data = mam_temp, family = gaussian(), 
  data2 = list(A_temp = A_temp),
  prior = c(
    prior(normal(0, 0.5), class =  Intercept),
    prior(exponential(3), class = sd, group = "phylo"),
    prior(exponential(1), class = sd, group = "species")),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 1, iter = 4000, warmup = 500
)

summary(temp_test)
plot(temp_test)

precip_test <- brm(
  coef_precip | se(sqrt(1 / (n_obs - 3))) ~ 1 + (1|gr(phylo, cov = A_precip)) + (1| species), 
  data = mam_precip, family = gaussian(), 
  data2 = list(A_precip = A_precip),
  prior = c(
    prior(normal(0, 0.5), class =  Intercept),
    prior(exponential(3), class = sd, group = "phylo"),
    prior(exponential(3), class = sd, group = "species")),
  control = list(adapt_delta = 0.99,
                 max_treedepth = 15),
  chains = 1, iter = 4000, warmup = 500
)

summary(precip_test)
plot(precip_test)


### Question that remains is whether we need to have the residual variance in full with an OLRE,
### or whether our species grouping level term is enough. Need to know




