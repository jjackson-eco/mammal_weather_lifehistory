##############################################################
##                                                          ##
##          Global climate and population dynamics          ##
##                                                          ##
## Spatially autocorrelated meta-regression for Temperature ##
##                                                          ##
##                     Feb 25th 2021                        ##
##                                                          ##
##############################################################

## Using a nearest neighbours spatial weighting to fit a spatially autocorrelated
## meta-regression model on the GAM weather coefficients.

rm(list = ls())
options(width = 100)

library(tidyverse)
library(sf)
library(spdep)
library(brms)
library(tidybayes)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 1. Load data ####

load("data/mammal_analysis_data_GAM.RData", verbose = TRUE)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 2. Nearest neighbours spatial weights matrix for analysis ####

# specify the coordinates of the data
coordinates(mam_coef) <- ~ Longitude + Latitude

# return k nearest neighbours for each coordinate point
knea <- knearneigh(coordinates(mam_coef), longlat = TRUE)

# convert to a neighbours list
neighbours <- knn2nb(knea)

# Spatial weighting matrix
Wmat <- nb2mat(neighbours, style = "W")

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3. Temperature model ####

## Base model
set.seed(666)
temp_base <- brm(
  coef_temp ~ 1 + sample_size + biome + (1| species),  
  data = mam_coef, family = gaussian(),
  prior = c( # lagsar gets a flat prior bound between 0 and 1
    prior(normal(0, 0.3), class =  Intercept),
    prior(normal(0, 0.3), class = b, coef = "sample_size"),
    prior(exponential(8), class = sd, group = "species")),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 4000, warmup = 2000
)

## spatial model
set.seed(666)
temp_sp <- brm(
  coef_temp ~ 1 + sample_size + biome + sar(Wmat, type = "lag") + (1| species),  
  data = mam_coef, family = gaussian(),
  data2 = list(Wmat = Wmat),
  prior = c( # lagsar gets a flat prior bound between 0 and 1
    prior(normal(0, 0.3), class =  Intercept),
    prior(normal(0, 0.3), class = b, coef = "sample_size"),
    prior(exponential(8), class = sd, group = "species")),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 4000, warmup = 2000
)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 4 Model comparisons ####

## Model comparisons
temp_base <- add_criterion(temp_base, criterion = c("loo","waic"))
temp_biome <- add_criterion(temp_biome, criterion = c("loo","waic"))

mod_comp_temp <- as.data.frame(loo_compare(temp_base, temp_biome, criterion = "loo"))

