####################################################
##                                                ##
##     Global climate and population dynamics     ##
##                                                ##
##     Weather Variance with GAM coefficients     ##
##                                                ##
##        brms Phylogenetic meta-regression       ##
##              and spatial effects               ##
##                                                ##
##                 Aug 23rd 2021                  ##
##                                                ##
####################################################

## Investigating the effects of biome on population responses to weather variance
## from GAM ARMA models. Implementing the meta-regression framework for both temp and precip

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
library(tidybayes)
library(viridis)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 1. Load data ####

## Weather effects
load("data/pgr_weather/var_5km_GAM.RData", verbose = TRUE)
glimpse(var_5km_GAM)

## Phylogenetic
load("../Phlogenetics/mam_lpd_pruned.RData", verbose = TRUE)
str(mamMCC_pruned)

## Life history data
load("data/lifehistory.RData", verbose = TRUE)

## Species names to link data
load("../rawdata/GBIF_species_names_mamUPDATE.RData", verbose = TRUE)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 2. Mammal coefficient data ####

mam_weathervar <- var_5km_GAM %>% 
  ungroup() %>%  ## <- The groups ruin the z-transformations
  left_join(x = ., y = dplyr::select(lpd_gbif, Binomial, gbif.species.id), 
            by = "Binomial") %>% 
  left_join(x = ., y = dplyr::select(lifehistory, c(3,6:9)),
            by = "gbif.species.id") %>% 
  mutate(species = gsub(" ", "_", gbif_species),
         phylo = species,
         # z transform the coefficients
         coef_temp = as.numeric(scale(coef_tempvar)),
         coef_precip = as.numeric(scale(coef_precipvar)),
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
  drop_na(litter, longevity, bodymass, coef_precip)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3. Phylogenetic covariance matrix ####

# Trim tree to right data
mamMCC_weathervar <- keep.tip(mamMCC_pruned, mam_weathervar$phylo) 

# Covariance matrix - Brownian motion model
A_weathervar <- ape::vcv.phylo(mamMCC_weathervar)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 4. Temperature Variance models ####

## Base model
set.seed(666)
temp_base <- brm(
  coef_temp ~ 1 + sample_size + (1|gr(phylo, cov = A_weathervar)) + (1| species),  
  data = mam_weathervar, family = gaussian(),
  data2 = list(A_weathervar = A_weathervar),
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
  coef_temp ~ 1 + biome + sample_size + (1|gr(phylo, cov = A_weathervar)) + (1| species),  
  data = mam_weathervar, family = gaussian(),
  data2 = list(A_weathervar = A_weathervar),
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
#### 5. Precipitation Variance models ####

## Base model
set.seed(666)
precip_base <- brm(
  coef_precip ~ 1 + sample_size + (1|gr(phylo, cov = A_weathervar)) + (1| species),  
  data = mam_weathervar, family = gaussian(),
  data2 = list(A_weathervar = A_weathervar),
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
precip_biome <- brm(
  coef_precip ~ 1 + biome + sample_size + (1|gr(phylo, cov = A_weathervar)) + (1| species),  
  data = mam_weathervar, family = gaussian(),
  data2 = list(A_weathervar = A_weathervar),
  prior = c(
    prior(normal(0, 0.3), class =  Intercept),
    prior(normal(0, 0.15), class = b),
    prior(normal(0, 0.3), class = b, coef = "sample_size"),
    prior(exponential(8), class = sd, group = "phylo"),
    prior(exponential(8), class = sd, group = "species")),
  control = list(adapt_delta = 0.98,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 4000, warmup = 2000)

#_______________________________________________________________________________
### 5b. Model comparisons

## Model comparisons
precip_base <- add_criterion(precip_base, criterion = c("loo","waic"))
precip_biome <- add_criterion(precip_biome, criterion = c("loo","waic"))

mod_comp_precip <- as.data.frame(loo_compare(precip_base, precip_biome, criterion = "loo"))

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 5. Save ####

save(mod_comp_temp, file = "results/gaussian_models/model_comparison_tempvar_rawcoef.RData")
save(temp_biome, file = "results/gaussian_models/tempvar_biome_rawcoef.RData")

save(mod_comp_precip, file = "results/gaussian_models/model_comparison_precipvar_rawcoef.RData")
save(precip_biome, file = "results/gaussian_models/precipvar_biome_rawcoef.RData")

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 6. Model plots ####

temp_colour <- "#990a80"
precip_colour <- "#287f79"

## Temperature
temp_plot <- temp_biome %>%
  gather_draws(`b_Intercept|sd_.*|b_sample_size|sigma`, regex = TRUE) %>% #tidybayes
  ungroup() %>% 
  ggplot(aes(y = .variable, x = .value)) + 
  stat_halfeye(show.legend = FALSE, fill = temp_colour) +
  geom_vline(xintercept = 0, size = 0.8) +
  scale_y_discrete(labels = c(expression(paste("Global intercept ", bar(alpha))),
                              expression(paste("Sample size ", beta[N])),
                              expression(paste("Phylogenetic covariance ", sigma[PHYLO])),
                              expression(paste("Species level variance ", sigma[SPECIES])),
                              "Population-level variance")) +
  labs(x = "Posterior estimate", y = NULL, tag = "a)") +
  theme_ridges(center_axis_labels = TRUE, grid = T, line_size = 0.3) +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 10))

## Precipitation
precip_plot <- precip_biome %>%
  gather_draws(`b_Intercept|sd_.*|b_sample_size|sigma`, regex = TRUE) %>% #tidybayes
  ungroup() %>% 
  ggplot(aes(y = .variable, x = .value)) + 
  stat_halfeye(show.legend = FALSE, fill = precip_colour) +
  geom_vline(xintercept = 0, size = 0.8) +
  scale_y_discrete(labels = c(expression(paste("Global intercept ", bar(alpha))),
                              expression(paste("Sample size ", beta[N])),
                              expression(paste("Phylogenetic covariance ", sigma[PHYLO])),
                              expression(paste("Species level variance ", sigma[SPECIES])),
                              "Population-level variance")) +
  labs(x = "Posterior estimate", y = NULL, tag = "b)") +
  theme_ridges(center_axis_labels = TRUE, grid = T, line_size = 0.3) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

temp_plot + precip_plot

ggsave(temp_plot + precip_plot,
       filename = "plots/meta_regression/weathervar_posterior.jpeg",
         width = 35, height = 14, units = "cm", dpi = 500)



