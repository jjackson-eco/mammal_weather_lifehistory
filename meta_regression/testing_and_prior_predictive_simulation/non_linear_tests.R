####################################################
##                                                ##
##     Global climate and population dynamics     ##
##                                                ##
##        brms Phylogenetic meta-regression       ##
##            non-linear model testing            ##
##                                                ##
##                Nov 13th 2020                   ##
##                                                ##
####################################################

## Investigating potential non-linear effects between life-history and weather effects

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
load("data/pgr_weather/mnanom_5km.RData", verbose = TRUE)
glimpse(mnanom_5km)

## Phylogenetic
load("../Phlogenetics/mam_lpd_pruned.RData", verbose = TRUE)
str(mamMCC_pruned)

## Life history data
load("data/lifehistory.RData", verbose = TRUE)

## Species names to link data
load("../rawdata/GBIF_species_names_mamUPDATE.RData", verbose = TRUE)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 2. Mammal coefficient data ####

mam_temp <- mnanom_5km %>% 
  left_join(x = ., y = dplyr::select(lpd_gbif, Binomial, gbif.species.id), 
            by = "Binomial") %>% 
  left_join(x = ., y = dplyr::select(lifehistory, c(3,6:9)),
            by = "gbif.species.id") %>% 
  mutate(species = gsub(" ", "_", gbif_species),
         phylo = species,
         # z transform the coefficients
         coef_temp = as.numeric(scale(coef_temp)),
         coef_precip = as.numeric(scale(coef_precip)),
         # z transform absolute latitude for models
         lat = as.numeric(scale(abs(Latitude))),
         # observation-level term for residual variance (not sure if needed)
         OLRE = 1:n(),
         # iucn change  
         IUCNstatus = if_else(is.na(IUCNstatus) == T, "NotAss", IUCNstatus),
         sample_size = as.numeric(scale(log(n_obs))),
         # longevity as a factor
         longevity_factor = if_else(longevity > mean(longevity,na.rm = TRUE), "high", "low")) %>% 
  dplyr::select(id = ID, id_block = ID_block, n_obs, sample_size, 
                order = Order, species, phylo, 
                biome, lat, iucn = IUCNstatus, litter,
                longevity, longevity_factor, bodymass, coef_temp, 
                coef_precip) %>% 
  drop_na(litter, longevity, bodymass)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3. Phylogenetic covariance matrix ####

# Trim tree to right data
mamMCC_temp <- keep.tip(mamMCC_pruned, mam_temp$phylo) 

# Covariance matrix - Brownian motion model
A_temp <- ape::vcv.phylo(mamMCC_temp)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 4. Non-linear model test ####

#_______________________________________________________________________________
### 4a. Simplest models with no random effects

# base
set.seed(666)
base_lon_simple <- brm(
  coef_temp ~ 1,  
  data = mam_temp, family = gaussian(),
  prior = c(
    prior(normal(0, 0.2), class =  Intercept)),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 1, iter = 2000, warmup = 500
)

# linear effect
set.seed(666)
linear_lon_simple <- brm(
  coef_temp ~ 1 + longevity,  
  data = mam_temp, family = gaussian(),
  prior = c(
    prior(normal(0, 0.2), class =  Intercept),
    prior(normal(0, 0.2), class = b, coef = "longevity")),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 1, iter = 2000, warmup = 500
)

# threshold model
set.seed(666)
thresh_lon_simple <- brm(
  coef_temp ~ 1 + longevity*longevity_factor,  
  data = mam_temp, family = gaussian(),
  prior = c(
    prior(normal(0, 0.2), class =  Intercept),
    prior(normal(0, 0.2), class = b)),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 1, iter = 2000, warmup = 500
)

# simple exponential model
set.seed(666)
nonlinear_lon_simple <- brm(
  bf(coef_temp ~ V0 * exp(r * longevity), V0 + r ~ 1, nl = TRUE),  
  data = mam_temp, family = gaussian(),
  prior = c(
    prior(normal(0, 0.1), nlpar = "V0"),
    prior(normal(-1.3, 0.1), nlpar = "r")),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 1, iter = 2000, warmup = 500
)

## Model comparisons
base_lon_simple <- add_criterion(base_lon_simple, criterion = c("loo","waic"))
linear_lon_simple <- add_criterion(linear_lon_simple, criterion = c("loo","waic"))
thresh_lon_simple <- add_criterion(thresh_lon_simple, criterion = c("loo","waic"))
nonlinear_lon_simple <- add_criterion(nonlinear_lon_simple, criterion = c("loo","waic"))

loo_compare(base_lon_simple, linear_lon_simple, 
            thresh_lon_simple, nonlinear_lon_simple,
            criterion = "loo")

## Potential support for non-linear diminishing returns exponential model, but lets see

# Quick plot of the results
tibble(longevity = seq(-2,2, length.out = 1000)) %>%
  mutate(longevity_factor = if_else(longevity > mean(longevity,na.rm = TRUE), 
                                    "high", "low")) %>% 
  mutate(post_lin = colMeans(posterior_predict(linear_lon_simple, newdata = .)),
         post_nlin = colMeans(posterior_predict(nonlinear_lon_simple, newdata = .)),
         post_thresh = colMeans(posterior_predict(thresh_lon_simple, newdata = .))) %>% 
  ggplot(aes(x = longevity, y = post_lin)) +
  geom_point(data = mam_temp, aes(x = longevity, y = coef_temp), alpha = 0.3, size = 3) +
  geom_line() +
  geom_line(aes(y = post_thresh, group = longevity_factor), colour = "red") +
  geom_line(aes(y = post_nlin), colour = "green") +
  coord_cartesian(ylim = c(-1,1))

#_______________________________________________________________________________
### 4b. Mixed effects

# Linear
set.seed(666)
temp_lon <- brm(
  coef_temp ~ 1 + longevity + sample_size + (1|gr(phylo, cov = A_temp)) + (1| species),  
  data = mam_temp, family = gaussian(),
  data2 = list(A_temp = A_temp),
  prior = c(
    prior(normal(0, 0.2), class =  Intercept),
    prior(normal(0, 0.2), class = b, coef = "longevity"),
    prior(normal(0, 0.5), class = b, coef = "sample_size"),
    prior(exponential(5), class = sd, group = "phylo"),
    prior(exponential(5), class = sd, group = "species")),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 4, cores = 4, iter = 3000, warmup = 500
)

# Simple exponential
set.seed(666)
nonlinear_lon <- brm(
  bf(coef_temp ~ V0 * exp(r * longevity), 
     V0 ~ 1 + sample_size + (1|gr(phylo, cov = A_temp)) + (1| species),
     r ~ 1, nl = TRUE),  
  data = mam_temp, family = gaussian(),
  data2 = list(A_temp = A_temp),
  prior = c(
    prior(normal(0, 0.1), nlpar = "V0"),
    prior(normal(-1.3, 0.1), nlpar = "r"),
    prior(normal(0, 0.5), class = b, coef = "sample_size", nlpar = "V0"),
    prior(exponential(5), class = sd, group = "phylo", nlpar = "V0"),
    prior(exponential(5), class = sd, group = "species", nlpar = "V0")),
  control = list(adapt_delta = 0.99,
                 max_treedepth = 15),
  chains = 4, cores = 4, iter = 3000, warmup = 500
)

# Model comparions
temp_lon <- add_criterion(temp_lon, criterion = c("loo","waic"))
nonlinear_lon <- add_criterion(nonlinear_lon, criterion = c("loo","waic"))


loo_compare(base_lon_simple, linear_lon_simple, 
            thresh_lon_simple, nonlinear_lon_simple,
            temp_lon, nonlinear_lon,
            criterion = "loo")


