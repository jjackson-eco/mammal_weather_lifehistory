####################################################
##                                                ##
##     Global climate and population dynamics     ##
##                                                ##
##        Abundance measure and IUCN status       ##
##                                                ##
##        brms Phylogenetic meta-regression       ##
##              and spatial effects               ##
##                                                ##
##                 Aug 31st 2021                  ##
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
library(gghalves)
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

mam_coef<- mnanom_5km_GAM %>% 
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
                biome, lat, iucn = IUCNstatus, abundance_measure, litter,
                longevity, bodymass, coef_temp, 
                coef_precip) %>% 
  drop_na(litter, longevity, bodymass, coef_precip, abundance_measure) ## Getting rid of the NA rows just for this verification

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3. Phylogenetic covariance matrix ####

# Trim tree to right data
mamMCC <- keep.tip(mamMCC_pruned, mam_coef$phylo) 

# Covariance matrix - Brownian motion model
A_mam <- ape::vcv.phylo(mamMCC)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 4.Base models ####

set.seed(666)
temp_base <- brm(
  coef_temp ~ 1 + sample_size + (1|gr(phylo, cov = A_mam)) + (1| species),  
  data = mam_coef, family = gaussian(),
  data2 = list(A_mam = A_mam),
  prior = c(
    prior(normal(0, 0.2), class =  Intercept),
    prior(normal(0, 0.3), class = b, coef = "sample_size"),
    prior(exponential(6), class = sd, group = "phylo"),
    prior(exponential(9), class = sd, group = "species")),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 4000, warmup = 2000
)

set.seed(666)
precip_base <- brm(
  coef_precip ~ 1 + sample_size + (1|gr(phylo, cov = A_mam)) + (1| species),  
  data = mam_coef, family = gaussian(),
  data2 = list(A_mam = A_mam),
  prior = c(
    prior(normal(0, 0.3), class =  Intercept),
    prior(normal(0, 0.3), class = b, coef = "sample_size"),
    prior(exponential(8), class = sd, group = "phylo"),
    prior(exponential(8), class = sd, group = "species")),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 4000, warmup = 2000
)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 5. Temperature models ####

set.seed(666)
temp_measure <- brm(
  coef_temp ~ 1 + sample_size + abundance_measure + (1|gr(phylo, cov = A_mam)) + (1| species),  
  data = mam_coef, family = gaussian(),
  data2 = list(A_mam = A_mam),
  prior = c(
    prior(normal(0, 0.2), class =  Intercept),
    prior(normal(0, 0.3), class = b, coef = "sample_size"),
    prior(normal(0, 0.3), class = b),
    prior(exponential(6), class = sd, group = "phylo"),
    prior(exponential(9), class = sd, group = "species")),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 4000, warmup = 2000
)

set.seed(666)
temp_iucn <- brm(
  coef_temp ~ 1 + sample_size + iucn + (1|gr(phylo, cov = A_mam)) + (1| species),  
  data = mam_coef, family = gaussian(),
  data2 = list(A_mam = A_mam),
  prior = c(
    prior(normal(0, 0.2), class =  Intercept),
    prior(normal(0, 0.3), class = b, coef = "sample_size"),
    prior(normal(0, 0.3), class = b),
    prior(exponential(6), class = sd, group = "phylo"),
    prior(exponential(9), class = sd, group = "species")),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 4000, warmup = 2000
)

#_______________________________________________________________________________
### 5b. Model comparisons

## Model comparisons
temp_base <- add_criterion(temp_base, criterion = c("loo","waic"))
temp_measure <- add_criterion(temp_measure, criterion = c("loo","waic"))#
temp_iucn <- add_criterion(temp_iucn, criterion = c("loo","waic"))

mod_comp_temp <- as.data.frame(loo_compare(temp_base, temp_measure, temp_iucn, criterion = "loo"))

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 6. Precipitation models ####

set.seed(666)
precip_measure <- brm(
  coef_precip ~ 1 + sample_size + abundance_measure + (1|gr(phylo, cov = A_mam)) + (1| species),  
  data = mam_coef, family = gaussian(),
  data2 = list(A_mam = A_mam),
  prior = c(
    prior(normal(0, 0.3), class =  Intercept),
    prior(normal(0, 0.3), class = b, coef = "sample_size"),
    prior(normal(0, 0.3), class = b),
    prior(exponential(8), class = sd, group = "phylo"),
    prior(exponential(8), class = sd, group = "species")),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 4000, warmup = 2000
)

set.seed(666)
precip_iucn <- brm(
  coef_precip ~ 1 + sample_size + iucn + (1|gr(phylo, cov = A_mam)) + (1| species),  
  data = mam_coef, family = gaussian(),
  data2 = list(A_mam = A_mam),
  prior = c(
    prior(normal(0, 0.3), class =  Intercept),
    prior(normal(0, 0.3), class = b, coef = "sample_size"),
    prior(normal(0, 0.3), class = b),
    prior(exponential(8), class = sd, group = "phylo"),
    prior(exponential(8), class = sd, group = "species")),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 4000, warmup = 2000
)

#_______________________________________________________________________________
### 5b. Model comparisons

## Model comparisons
precip_base <- add_criterion(precip_base, criterion = c("loo","waic"))
precip_measure <- add_criterion(precip_measure, criterion = c("loo","waic"))#
precip_iucn <- add_criterion(precip_iucn, criterion = c("loo","waic"))

mod_comp_precip <- as.data.frame(loo_compare(precip_base, precip_measure, precip_iucn, criterion = "loo"))

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 7. Plots of the posterior ####

temp_colour <- "#d45371"
precip_colour <- "#30738e"

temp_measure_plot <- as.data.frame(brms::posterior_predict(temp_measure, nsamples = 2000)) %>% # reducing to make data more manageable 
  mutate(sim = 1:2000) %>% 
  pivot_longer(-sim) %>% 
  bind_cols(., slice(mam_coef, rep(1:nrow(mam_coef), 2000))) %>% 
  ggplot(aes(x = abundance_measure, y = value)) +
  geom_hline(yintercept = 0) +
  stat_halfeye(width = 0.4, alpha = 0.7, fill = temp_colour) +
  geom_half_point(data = mam_coef, 
                  aes(y = coef_temp, x = abundance_measure), 
                      colour = temp_colour,
                  side = "l", ## draw jitter on the left
                  range_scale = 0.4, ## control range of jitter
                  alpha = 0.7, size = 1) +
  scale_x_discrete(labels = c("Density", "Estimate", "Full\npopulation\ncount",
                              "Index", "Monitoring\nper unit\neffort", "Proxy")) +
  labs(x = "Abundance measure", y = "Temperature effect") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank())

precip_iucn_plot <- as.data.frame(brms::posterior_predict(precip_iucn, nsamples = 2000)) %>% # reducing to make data more manageable 
  mutate(sim = 1:2000) %>% 
  pivot_longer(-sim) %>% 
  bind_cols(., slice(mam_coef, rep(1:nrow(mam_coef), 2000))) %>% 
  ggplot(aes(x = iucn, y = value)) +
  geom_hline(yintercept = 0) +
  stat_halfeye(width = 0.4, alpha = 0.7, fill = precip_colour) +
  geom_half_point(data = mam_coef, 
                  aes(y = coef_temp, x = iucn), 
                  colour = precip_colour,
                  side = "l", ## draw jitter on the left
                  range_scale = 0.4, ## control range of jitter
                  alpha = 0.7, size = 1) +
  scale_x_discrete(labels = c("Critically\nendangered", "Endangered", "Least\nconcern",
                              "Not Assessed", "Near\nthreatened", "Vulnerable")) +
  labs(x = "IUCN status", y = "Precipitation effect") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank())


ggsave(temp_measure_plot/ precip_iucn_plot,
       filename = "plots/meta_regression/abundance_measure_iucn_exploration.jpeg",
       width = 22, height = 20, units = "cm", dpi = 600)




