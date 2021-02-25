##########################################################
##                                                      ##
##       Global climate and population dynamics         ##
##                                                      ##
##  Absolute precipitation effect with GAM - log normal ##
##                                                      ##
##          brms Phylogenetic meta-regression           ##
##          life-history and spatial effects            ##
##                                                      ##
##                   Dec 16th 2020                      ##
##                                                      ##
##########################################################

## Investigating the effects of biome and life-history on absolute population responses
## from GAM ARMA models. Implementing the meta-regression framework for Precipitation. Using log-normal models

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

mam_precip <- mnanom_5km_GAM %>% 
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
         # logged absolute coefficients
         log_abs_temp = log(abs_temp),
         log_abs_precip = log(abs_precip),
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
                abs_precip, log_abs_temp, log_abs_precip) %>% 
  drop_na(litter, longevity, bodymass, log_abs_precip)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3. Phylogenetic covariance matrix ####

# Trim tree to right data
mamMCC_precip <- keep.tip(mamMCC_pruned, mam_precip$phylo) 

# Covariance matrix - Brownian motion model
A_precip <- ape::vcv.phylo(mamMCC_precip)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 4. Precipitation models ####

# Base model chain test
set.seed(666)
precip_base_test <- brm(
  log_abs_precip ~ 1 + sample_size + (1|gr(phylo, cov = A_precip)) + (1| species),  
  data = mam_precip, family = gaussian(),
  data2 = list(A_precip = A_precip),
  prior = c(
    prior(normal(0, 1), class =  Intercept),
    prior(normal(0, 1), class = b, coef = "sample_size"),
    prior(exponential(6), class = sd, group = "phylo"),
    prior(exponential(5), class = sd, group = "species")),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 1, iter = 4000, warmup = 1000
)

# Base model 
set.seed(666)
precip_base <- brm(
  log_abs_precip ~ 1 + sample_size + (1|gr(phylo, cov = A_precip)) + (1| species),  
  data = mam_precip, family = gaussian(),
  data2 = list(A_precip = A_precip),
  prior = c(
    prior(normal(0, 1), class =  Intercept),
    prior(normal(0, 1), class = b, coef = "sample_size"),
    prior(exponential(6), class = sd, group = "phylo"),
    prior(exponential(5), class = sd, group = "species")),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

#_______________________________________________________________________________
### 4a. Univariate predictor models

# Longevity
set.seed(666)
precip_longevity <- brm(
  log_abs_precip ~ 1 + sample_size + longevity + (1|gr(phylo, cov = A_precip)) + (1| species),  
  data = mam_precip, family = gaussian(),
  data2 = list(A_precip = A_precip),
  prior = c(
    prior(normal(0, 0.8), class =  Intercept),
    prior(normal(0, 0.6), class = b, coef = "sample_size"),
    prior(normal(0, 0.5), class = b, coef = "longevity"),
    prior(exponential(6), class = sd, group = "phylo"),
    prior(exponential(5), class = sd, group = "species")),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

# Litter size
set.seed(666)
precip_litter <- brm(
  log_abs_precip ~ 1 + sample_size + litter + (1|gr(phylo, cov = A_precip)) + (1| species),  
  data = mam_precip, family = gaussian(),
  data2 = list(A_precip = A_precip),
  prior = c(
    prior(normal(0, 1), class =  Intercept),
    prior(normal(0, 0.8), class = b, coef = "sample_size"),
    prior(normal(0, 0.8), class = b, coef = "litter"),
    prior(exponential(6), class = sd, group = "phylo"),
    prior(exponential(5), class = sd, group = "species")),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

# Adult body mass
set.seed(666)
precip_bodymass <- brm(
  log_abs_precip ~ 1 + sample_size + bodymass + (1|gr(phylo, cov = A_precip)) + (1| species),  
  data = mam_precip, family = gaussian(),
  data2 = list(A_precip = A_precip),
  prior = c(
    prior(normal(0, 1), class =  Intercept),
    prior(normal(0, 0.8), class = b, coef = "sample_size"),
    prior(normal(0, 0.8), class = b, coef = "bodymass"),
    prior(exponential(6), class = sd, group = "phylo"),
    prior(exponential(5), class = sd, group = "species")),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

# Biome
set.seed(666)
precip_biome <- brm(
  log_abs_precip ~ 1 + sample_size + biome + (1|gr(phylo, cov = A_precip)) + (1| species),  
  data = mam_precip, family = gaussian(),
  data2 = list(A_precip = A_precip),
  prior = c(
    prior(normal(0, 1), class =  Intercept),
    prior(normal(0, 0.8), class = b, coef = "sample_size"),
    prior(normal(0, 0.3), class = b), # more conservative as lots of category levels
    prior(exponential(6), class = sd, group = "phylo"),
    prior(exponential(5), class = sd, group = "species")),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

#_______________________________________________________________________________
### 4b. Univariate model comparisons

## Model comparisons
precip_base <- add_criterion(precip_base, criterion = c("loo","waic"))
precip_longevity <- add_criterion(precip_longevity, criterion = c("loo","waic"))
precip_litter <- add_criterion(precip_litter, criterion = c("loo","waic"))
precip_bodymass <- add_criterion(precip_bodymass, criterion = c("loo","waic"))
precip_biome <- add_criterion(precip_biome, criterion = c("loo","waic"))

as.data.frame(loo_compare(precip_base, precip_longevity, precip_bodymass, 
                          precip_litter, precip_biome, criterion = "loo"))

#_______________________________________________________________________________
### 4c. Life-history models full

# More conservative priors due to lots of parameters
lh_priors <- c(
  prior(normal(0, 0.5), class =  Intercept),
  prior(normal(0, 0.3), class = b),
  prior(normal(0, 0.3), class = b, coef = "sample_size"),
  prior(exponential(8), class = sd, group = "phylo"),
  prior(exponential(5), class = sd, group = "species"))

## Life-history full
set.seed(666)
precip_lh <- brm(
  log_abs_precip ~ 1 + longevity + bodymass + litter +
    longevity:bodymass + litter:bodymass + litter:longevity + 
    sample_size + (1|gr(phylo, cov = A_precip)) + (1| species),  
  data = mam_precip, 
  family = gaussian(), 
  data2 = list(A_precip = A_precip),
  prior = lh_priors,
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## Life-history univariate only
set.seed(666)
precip_lh_uni <- brm(
  log_abs_precip ~ 1 + longevity + bodymass + litter +
    sample_size + (1|gr(phylo, cov = A_precip)) + (1| species),  
  data = mam_precip, 
  family = gaussian(), 
  data2 = list(A_precip = A_precip),
  prior = lh_priors,
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## Longevity and Body Mass full
set.seed(666)
precip_lonbod <- brm(
  log_abs_precip ~ 1 + longevity*bodymass + sample_size + (1|gr(phylo, cov = A_precip)) + (1| species),  
  data = mam_precip, 
  family = gaussian(), 
  data2 = list(A_precip = A_precip),
  prior = lh_priors,
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## Litter and Body Mass full
set.seed(666)
precip_litbod <- brm(
  log_abs_precip ~ 1 + litter*bodymass + sample_size + (1|gr(phylo, cov = A_precip)) + (1| species),  
  data = mam_precip, 
  family = gaussian(),
  data2 = list(A_precip = A_precip),
  prior = lh_priors,
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## Body mass correction - additive only
set.seed(666)
precip_lonbod_simple <- brm(
  log_abs_precip ~ 1 + longevity + bodymass + sample_size + (1|gr(phylo, cov = A_precip)) + (1| species),  
  data = mam_precip, 
  family = gaussian(), 
  data2 = list(A_precip = A_precip),
  prior = lh_priors,
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

set.seed(666)
precip_litbod_simple <- brm(
  log_abs_precip ~ 1 + litter + bodymass + sample_size + (1|gr(phylo, cov = A_precip)) + (1| species),  
  data = mam_precip, 
  family = gaussian(), 
  data2 = list(A_precip = A_precip),
  prior = lh_priors,
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

#_______________________________________________________________________________
### 4b. Full model comparisons

precip_lh <- add_criterion(precip_lh, criterion = c("loo","waic"))
precip_lh_uni <- add_criterion(precip_lh_uni, criterion = c("loo","waic"))
precip_lonbod <- add_criterion(precip_lonbod, criterion = c("loo","waic"))
precip_litbod<- add_criterion(precip_litbod, criterion = c("loo","waic"))
precip_lonbod_simple <- add_criterion(precip_lonbod_simple, criterion = c("loo","waic"))
precip_litbod_simple <- add_criterion(precip_litbod_simple, criterion = c("loo","waic"))

precip_modcomp_lognormal <- as.data.frame(loo_compare(precip_base, precip_longevity, precip_bodymass, 
                                                precip_litter, precip_biome, precip_lonbod_simple,
                                                precip_litbod_simple, precip_lh_uni, precip_lonbod, precip_litbod,
                                                precip_lh, criterion = "loo"))


# Looks like best is litter and body mass, but also strong evidence for longevity separately relative to base model.
# Not good evidence for interactions though. Therefore we pick the univariate only model.
precip_modcomp_lognormal
precip_lh_uni

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 4. Save ####

save(precip_lh_uni, file = "results/lognormal_models/precip_lh_uni.RData")
save(precip_modcomp_lognormal, file = "results/lognormal_models/precip_modcomp_lognormal.RData")


