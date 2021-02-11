####################################################
##                                                ##
##     Global climate and population dynamics     ##
##                                                ##
##      Absolute temperature effect with GAM      ##
##                                                ##
##        brms Phylogenetic meta-regression       ##
##        life-history and spatial effects        ##
##                                                ##
##                 Dec 16th 2020                  ##
##                                                ##
####################################################

## Investigating the effects of biome and life-history on absolute population responses
## from GAM ARMA models. Implementing the meta-regression framework for temperature

rm(list = ls())
options(width = 100)

library(tidyverse)
library(tidybayes)
library(ape)
library(brms)
library(rethinking)
library(nlme)
library(patchwork)
library(ggridges)
library(ggdist)
library(viridis)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 1. Load data ####

## Weather effects
load("data/pgr_weather/mnanom_5km_GAM.RData", verbose = TRUE)
glimpse(mnanom_5km_GAM)

## Phylogenetic
load("../Phlogenetics/mam_lpd_pruned.RData", verbose = TRUE)
str(mamMCC_pruned)

## Life history data
load("data/lifehistory.RData", verbose = TRUE)

## Species names to link data
load("../rawdata/GBIF_species_names_mamUPDATE.RData", verbose = TRUE)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 2. Mammal coefficient data ####

mam_temp <- mnanom_5km_GAM %>% 
  ungroup() %>%  ## <- The groups ruin the z-transformations
  left_join(x = ., y = dplyr::select(lpd_gbif, Binomial, gbif.species.id), 
            by = "Binomial") %>% 
  left_join(x = ., y = dplyr::select(lifehistory, c(3,6:9)),
            by = "gbif.species.id") %>% 
  mutate(species = gsub(" ", "_", gbif_species),
         phylo = species,
         # absolute values of z transformed coefficients
         abs_temp = abs(as.numeric(scale(coef_temp))),
         abs_precip = abs(as.numeric(scale(coef_precip))),
         # z transform absolute latitude for models
         lat = as.numeric(scale(abs(Latitude))),
         # observation-level term for residual variance (not sure if needed)
         OLRE = 1:n(),
         # iucn change  
         IUCNstatus = if_else(is.na(IUCNstatus) == T, "NotAss", IUCNstatus),
         sample_size = as.numeric(scale(log(n_obs)))) %>% 
  dplyr::select(id = ID, id_block = ID_block, n_obs, sample_size, 
                order = Order, species, phylo, 
                biome, lat, iucn = IUCNstatus, litter,
                longevity, bodymass, abs_temp, 
                abs_precip) %>% 
  drop_na(litter, longevity, bodymass)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3. Phylogenetic covariance matrix ####

# Trim tree to right data
mamMCC_temp <- keep.tip(mamMCC_pruned, mam_temp$phylo) 

# Covariance matrix - Brownian motion model
A_temp <- ape::vcv.phylo(mamMCC_temp)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 4. Temperature models ####

# Base model chain test
plot(density(rgamma(1000,shape = 2, scale = 1)))
plot(density(rexp(1000, rate = 10)))

set.seed(666)
temp_base_test <- brm(
  abs_temp ~ 1 + sample_size + (1|gr(phylo, cov = A_temp)) + (1| species),
  data = mam_temp,
  family = Gamma(link = "log"),
  data2 = list(A_temp = A_temp),
  prior = c(
    prior(normal(0, 1), class =  Intercept),
    prior(normal(0, 1), class = b, coef = "sample_size"),
    prior(exponential(11), class = sd, group = "phylo"),
    prior(exponential(3), class = sd, group = "species"),
    prior(gamma(2,1), class = shape)),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 1, iter = 2000, warmup = 400
)

## Base model
set.seed(666)
temp_base <- brm(
  abs_temp ~ 1 + sample_size + (1|gr(phylo, cov = A_temp)) + (1| species),  
  data = mam_temp, 
  family = Gamma(link = "log"), 
  data2 = list(A_temp = A_temp),
  prior = c(
    prior(normal(0, 1), class =  Intercept),
    prior(normal(0, 1), class = b, coef = "sample_size"),
    prior(exponential(10), class = sd, group = "phylo"),
    prior(exponential(2), class = sd, group = "species"),
    prior(gamma(2,1), class = shape)),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 2000, warmup = 400
)

#_______________________________________________________________________________
### 4a. Univariate predictor models

## Longevity
set.seed(666)
temp_longevity <- brm(
  abs_temp ~ 1 + longevity + sample_size + (1|gr(phylo, cov = A_temp)) + (1| species),  
  data = mam_temp, 
  family = Gamma(link = "log"), 
  data2 = list(A_temp = A_temp),
  prior = c(
    prior(normal(0, 1), class =  Intercept),
    prior(normal(0, 1), class = b, coef = "longevity"),
    prior(normal(0, 1), class = b, coef = "sample_size"),
    prior(exponential(10), class = sd, group = "phylo"),
    prior(exponential(2), class = sd, group = "species"),
    prior(gamma(2,1), class = shape)),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 2000, warmup = 400
)

## Litter size
set.seed(666)
temp_litter <- brm(
  abs_temp ~ 1 + litter + sample_size + (1|gr(phylo, cov = A_temp)) + (1| species),  
  data = mam_temp, 
  family = Gamma(link = "log"), 
  data2 = list(A_temp = A_temp),
  prior = c(
    prior(normal(0, 1), class =  Intercept),
    prior(normal(0, 1), class = b, coef = "litter"),
    prior(normal(0, 1), class = b, coef = "sample_size"),
    prior(exponential(10), class = sd, group = "phylo"),
    prior(exponential(2), class = sd, group = "species"),
    prior(gamma(2,1), class = shape)),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 2000, warmup = 400
)

## Body mass
set.seed(666)
temp_bodymass <- brm(
  abs_temp ~ 1 + bodymass + sample_size + (1|gr(phylo, cov = A_temp)) + (1| species),  
  data = mam_temp, 
  family = Gamma(link = "log"), 
  data2 = list(A_temp = A_temp),
  prior = c(
    prior(normal(0, 1), class =  Intercept),
    prior(normal(0, 1), class = b, coef = "bodymass"),
    prior(normal(0, 1), class = b, coef = "sample_size"),
    prior(exponential(11), class = sd, group = "phylo"),
    prior(exponential(3), class = sd, group = "species"),
    prior(gamma(2,1), class = shape)),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 2000, warmup = 400
)

## Biome
set.seed(666)
temp_biome <- brm(
  abs_temp ~ 1 + biome + sample_size + (1|gr(phylo, cov = A_temp)) + (1| species),  
  data = mam_temp, 
  family = Gamma(link = "log"), 
  data2 = list(A_temp = A_temp),
  prior = c(
    prior(normal(0, 1), class =  Intercept),
    prior(normal(0, 1), class = b),
    prior(normal(0, 1), class = b, coef = "sample_size"),
    prior(exponential(11), class = sd, group = "phylo"),
    prior(exponential(3), class = sd, group = "species"),
    prior(gamma(2,1), class = shape)),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 20),
  chains = 3, cores = 3, iter = 2000, warmup = 400
)

#_______________________________________________________________________________
### 4b. Univariate model comparisons

## Model comparisons
temp_base <- add_criterion(temp_base, criterion = c("loo","waic"))
temp_longevity <- add_criterion(temp_longevity, criterion = c("loo","waic"))
temp_litter <- add_criterion(temp_litter, criterion = c("loo","waic"))
temp_bodymass <- add_criterion(temp_bodymass, criterion = c("loo","waic"))
temp_biome <- add_criterion(temp_biome, criterion = c("loo","waic"))

loo_compare(temp_base, temp_longevity, temp_bodymass, 
            temp_litter, temp_biome, criterion = "loo")


