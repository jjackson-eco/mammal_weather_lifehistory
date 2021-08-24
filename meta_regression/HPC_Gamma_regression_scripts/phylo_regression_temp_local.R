####################################################
##                                                ##
##     Global climate and population dynamics     ##
##                                                ##
##      Absolute temperature effect with GAM      ##
##                                                ##
##    Local - brms Phylogenetic meta-regression  ##
##        life-history and spatial effects        ##
##                                                ##
##                 Aug 23rd 2021                  ##
##                                                ##
####################################################

## Models regressing biome and life-history on absolute population responses
## from GAM ARMA models. Implementing the meta-regression framework for temperature

rm(list = ls())
options(width = 100)

# Packages
library(tidyverse)
library(tidybayes)
library(ape)
library(brms)
library(nlme)
library(patchwork)

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
plot(density(rgamma(1000,shape = 2, scale = 0.6)))
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
    prior(exponential(2), class = sd, group = "species"),
    prior(gamma(2,0.5), class = shape)),
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
    prior(exponential(11), class = sd, group = "phylo"),
    prior(exponential(2), class = sd, group = "species"),
    prior(gamma(2,0.5), class = shape)),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 4000, warmup = 2000
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
    prior(normal(0, 0.5), class =  Intercept),
    prior(normal(0, 0.3), class = b, coef = "longevity"),
    prior(normal(0, 0.5), class = b, coef = "sample_size"),
    prior(exponential(11), class = sd, group = "phylo"),
    prior(exponential(2), class = sd, group = "species"),
    prior(gamma(2,0.7), class = shape)),
  control = list(adapt_delta = 0.99,
                 max_treedepth = 17),
  chains = 3, cores = 3, iter = 4000, warmup = 2000
)

## Litter size
set.seed(666)
temp_litter <- brm(
  abs_temp ~ 1 + litter + sample_size + (1|gr(phylo, cov = A_temp)) + (1| species),  
  data = mam_temp, 
  family = Gamma(link = "log"), 
  data2 = list(A_temp = A_temp),
  prior = c(
    prior(normal(0, 0.8), class =  Intercept),
    prior(normal(0, 0.8), class = b, coef = "litter"),
    prior(normal(0, 0.8), class = b, coef = "sample_size"),
    prior(exponential(11), class = sd, group = "phylo"),
    prior(exponential(2), class = sd, group = "species"),
    prior(gamma(2,0.5), class = shape)),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 4000, warmup = 2000
)

## Body mass
set.seed(666)
temp_bodymass <- brm(
  abs_temp ~ 1 + bodymass + sample_size + (1|gr(phylo, cov = A_temp)) + (1| species),  
  data = mam_temp, 
  family = Gamma(link = "log"), 
  data2 = list(A_temp = A_temp),
  prior = c(
    prior(normal(0, 0.8), class =  Intercept),
    prior(normal(0, 0.5), class = b, coef = "bodymass"),
    prior(normal(0, 1), class = b, coef = "sample_size"),
    prior(exponential(11), class = sd, group = "phylo"),
    prior(exponential(2), class = sd, group = "species"),
    prior(gamma(2,0.5), class = shape)),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 4000, warmup = 2000
)

## Biome
set.seed(666)
temp_biome <- brm(
  abs_temp ~ 1 + biome + sample_size + (1|gr(phylo, cov = A_temp)) + (1| species),  
  data = mam_temp, 
  family = Gamma(link = "log"), 
  data2 = list(A_temp = A_temp),
  prior = c(
    prior(normal(0, 0.8), class =  Intercept),
    prior(normal(0, 0.3), class = b),
    prior(normal(0, 1), class = b, coef = "sample_size"),
    prior(exponential(11), class = sd, group = "phylo"),
    prior(exponential(2), class = sd, group = "species"),
    prior(gamma(2,0.5), class = shape)),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 4000, warmup = 2000
)

#_______________________________________________________________________________
### 4b. Univariate model comparisons

## Model comparisons
temp_base <- add_criterion(temp_base, criterion = c("loo","waic"))
temp_longevity <- add_criterion(temp_longevity, criterion = c("loo","waic"))
temp_litter <- add_criterion(temp_litter, criterion = c("loo","waic"))
temp_bodymass <- add_criterion(temp_bodymass, criterion = c("loo","waic"))
temp_biome <- add_criterion(temp_biome, criterion = c("loo","waic"))

univar_modcomp <- as.data.frame(loo_compare(temp_base, temp_longevity, temp_bodymass, 
            temp_litter, temp_biome, criterion = "loo"))

univar_modcomp
save(univar_modcomp, file = "results/local_gamma_models/temperature_univariate_model_comparisons.RData")

#_______________________________________________________________________________
### 4c. Life-history models full

lh_priors <- c(
  prior(normal(0, 0.3), class =  Intercept),
  prior(normal(0, 0.2), class = b),
  prior(normal(0, 0.3), class = b, coef = "sample_size"),
  prior(exponential(11), class = sd, group = "phylo"),
  prior(exponential(2), class = sd, group = "species"),
  prior(gamma(2,0.6), class = shape))

## Life-history full
set.seed(666)
temp_lh <- brm(
  abs_temp ~ 1 + longevity + bodymass + litter +
    longevity:bodymass + litter:bodymass + litter:longevity + 
    sample_size + (1|gr(phylo, cov = A_temp)) + (1| species),  
  data = mam_temp, 
  family = Gamma(link = "log"), 
  data2 = list(A_temp = A_temp),
  prior = lh_priors,
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 4000, warmup = 2000
)

## Life-history univariate only
set.seed(666)
temp_lh_uni <- brm(
  abs_temp ~ 1 + longevity + bodymass + litter +
    sample_size + (1|gr(phylo, cov = A_temp)) + (1| species),  
  data = mam_temp, 
  family = Gamma(link = "log"), 
  data2 = list(A_temp = A_temp),
  prior = lh_priors,
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 4000, warmup = 2000
)

## Longevity and Body Mass full
set.seed(666)
temp_lonbod <- brm(
  abs_temp ~ 1 + longevity*bodymass + sample_size + (1|gr(phylo, cov = A_temp)) + (1| species),  
  data = mam_temp, 
  family = Gamma(link = "log"), 
  data2 = list(A_temp = A_temp),
  prior = lh_priors,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 4000, warmup = 2000
)

## Litter and Body Mass full
set.seed(666)
temp_litbod <- brm(
  abs_temp ~ 1 + litter*bodymass + sample_size + (1|gr(phylo, cov = A_temp)) + (1| species),  
  data = mam_temp, 
  family = Gamma(link = "log"), 
  data2 = list(A_temp = A_temp),
  prior = lh_priors,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 4000, warmup = 2000
)

## Body mass correction - additive only
set.seed(666)
temp_lonbod_simple <- brm(
  abs_temp ~ 1 + longevity + bodymass + sample_size + (1|gr(phylo, cov = A_temp)) + (1| species),  
  data = mam_temp, 
  family = Gamma(link = "log"), 
  data2 = list(A_temp = A_temp),
  prior = lh_priors,
  control = list(adapt_delta = 0.98,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 4000, warmup = 2000
)

set.seed(666)
temp_litbod_simple <- brm(
  abs_temp ~ 1 + litter + bodymass + sample_size + (1|gr(phylo, cov = A_temp)) + (1| species),  
  data = mam_temp, 
  family = Gamma(link = "log"), 
  data2 = list(A_temp = A_temp),
  prior = lh_priors,
  control = list(adapt_delta = 0.98,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 4000, warmup = 2000
)

#_______________________________________________________________________________
### 4b. Full model comparisons

temp_lh <- add_criterion(temp_lh, criterion = c("loo","waic"))
temp_lh_uni <- add_criterion(temp_lh_uni, criterion = c("loo","waic"))
temp_lonbod <- add_criterion(temp_lonbod, criterion = c("loo","waic"))
temp_litbod<- add_criterion(temp_litbod, criterion = c("loo","waic"))
temp_lonbod_simple <- add_criterion(temp_lonbod_simple, criterion = c("loo","waic"))
temp_litbod_simple <- add_criterion(temp_litbod_simple, criterion = c("loo","waic"))

mod_comp <- as.data.frame(loo_compare(temp_base, temp_longevity, temp_bodymass, 
            temp_litter, temp_biome, temp_lonbod_simple,
            temp_litbod_simple, temp_lh_uni, temp_lonbod, temp_litbod,
            temp_lh, criterion = "loo"))

save(mod_comp, file = "results/local_gamma_models/temperature_model_comparisons.RData")

## Best model - temp_lh_uni
saveRDS(temp_lh_uni, "results/local_gamma_models/temp_lh_uni.RDS")


