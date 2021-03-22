####################################################
##                                                ##
##     Global climate and population dynamics     ##
##                                                ##
##        Temperature with GAM coefficients       ##
##                                                ##
##        brms Phylogenetic meta-regression       ##
##              and spatial effects               ##
##                                                ##
##                 Nov 9th 2020                   ##
##                                                ##
####################################################

## Investigating the effects of biome on population responses
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
library(flextable)

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
         # z transform the coefficients
         coef_temp = as.numeric(scale(coef_temp)),
         coef_precip = as.numeric(scale(coef_precip)),
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
                longevity, bodymass, coef_temp, 
                coef_precip) %>% 
  drop_na(litter, longevity, bodymass)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3. Phylogenetic covariance matrix ####

# Trim tree to right data
mamMCC_temp <- keep.tip(mamMCC_pruned, mam_temp$phylo) 

# Covariance matrix - Brownian motion model
A_temp <- ape::vcv.phylo(mamMCC_temp)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 4. Temperature models ####

## Base model
set.seed(666)
temp_base <- brm(
  coef_temp ~ 1 + sample_size + (1|gr(phylo, cov = A_temp)) + (1| species),  
  data = mam_temp, family = gaussian(),
  data2 = list(A_temp = A_temp),
  prior = c(
    prior(normal(0, 0.3), class =  Intercept),
    prior(normal(0, 0.3), class = b, coef = "sample_size"),
    prior(exponential(8), class = sd, group = "phylo"),
    prior(exponential(8), class = sd, group = "species")),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 4000, warmup = 2000
)

## Biome
set.seed(666)
temp_biome <- brm(
  coef_temp ~ 1 + biome + sample_size + (1|gr(phylo, cov = A_temp)) + (1| species),  
  data = mam_temp, family = gaussian(),
  data2 = list(A_temp = A_temp),
  prior = c(
    prior(normal(0, 0.3), class =  Intercept),
    prior(normal(0, 0.15), class = b),
    prior(normal(0, 0.3), class = b, coef = "sample_size"),
    prior(exponential(8), class = sd, group = "phylo"),
    prior(exponential(8), class = sd, group = "species")),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 4000, warmup = 2000
)

#_______________________________________________________________________________
### 4b. Model comparisons

## Model comparisons
temp_base <- add_criterion(temp_base, criterion = c("loo","waic"))
temp_biome <- add_criterion(temp_biome, criterion = c("loo","waic"))

mod_comp_temp <- as.data.frame(loo_compare(temp_base, temp_biome, criterion = "loo"))

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 5. Save ####

save(mod_comp_temp, file = "results/gaussian_models/model_comparison_temp_rawcoef.RData")
save(temp_biome, file = "results/gaussian_models/temp_biome_rawcoef.RData")

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 6. Table and model plots ####

load("results/gaussian_models/model_comparison_temp_rawcoef.RData")
load("results/gaussian_models/temp_biome_rawcoef.RData")

# Model comparison table
mod_comp_temp %>% 
  mutate(model = c("Base model", "Biome effect")) %>% 
  dplyr::select(model, elpd_loo, se_elpd_loo, elpd_diff, se_diff, looic) %>% 
  flextable(cwidth = 1.2) %>% 
  set_header_labels(model = "Model",
                    elpd_loo = "LOO elpd",
                    se_elpd_loo = "LOO elpd error",
                    elpd_diff = "elpd difference",
                    se_diff = "elpd error difference",
                    looic = "LOO information criterion") %>% 
  colformat_double(digits = 2) %>% 
  save_as_image("plots/meta_regression/temperature_gaussian_model_comparison.png")

# Model plot
temp_biome_pars <- parnames(temp_biome)

jpeg("plots/meta_regression/temp_biome_mod_parms.jpeg", 
     width = 13, height = 15, res = 500, units = "cm")
plot(temp_biome, pars = c("b_Intercept", "b_sample_size","sd_species__Intercept", 
                          "sd_phylo__Intercept", "sigma"))
dev.off()
