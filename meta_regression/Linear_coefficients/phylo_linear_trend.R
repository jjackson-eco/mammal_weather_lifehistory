
####################################################
##                                                ##
##     Global climate and population dynamics     ##
##                                                ##
##        Life history with trend coefficients    ##
##                                                ##
##                Dec 10th 2020                   ##
##                                                ##
####################################################

## Investigating the effects of biome and life-history on linear trend
## of population abundance. Naive linear trend term.

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

mam_trend <- mnanom_5km %>% 
  ungroup() %>%  ## <- The groups ruin the z-transformations
  left_join(x = ., y = dplyr::select(lpd_gbif, Binomial, gbif.species.id), 
            by = "Binomial") %>% 
  left_join(x = ., y = dplyr::select(lifehistory, c(3,6:9)),
            by = "gbif.species.id") %>% 
  mutate(species = gsub(" ", "_", gbif_species),
         phylo = species,
         # z transform the coefficient
         coef_trend = as.numeric(scale(coef_trend)),
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
                longevity, bodymass, coef_trend) %>% 
  drop_na(litter, longevity, bodymass)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3. Phylogenetic covariance matrix ####

# Trim tree to right data
mamMCC_trend <- keep.tip(mamMCC_pruned, mam_trend$phylo) 

# Covariance matrix - Brownian motion model
A_trend <- ape::vcv.phylo(mamMCC_trend)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 4. Exploring Trend relationship with lifehistory ####

mam_trend %>% 
  dplyr::select(id_block, n_obs, species, longevity, 
                litter, bodymass, coef_trend) %>% 
  pivot_longer(c(longevity, litter, bodymass), 
               names_to = "demovar") %>% 
  ggplot(aes(x = value, y = coef_trend, size = n_obs, colour = demovar)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", colour = "black") +
  facet_wrap(~ demovar, ncol = 3) +
  coord_cartesian(ylim = c(-3,5)) +
  scale_colour_viridis_d(option = "C", begin = 0.2, end = 0.8, guide = F) +
  scale_size_continuous(range = c(1,5), guide = F) +
  labs(x = "Scaled demographic value", y = "Linear population trend") +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank())

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 5. Trend models ####

## Base model
set.seed(666)
trend_base <- brm(
  coef_trend ~ 1 + sample_size + (1|gr(phylo, cov = A_trend)) + (1| species),  
  data = mam_trend, family = gaussian(),
  data2 = list(A_trend = A_trend),
  prior = c(
    prior(normal(0, 0.2), class =  Intercept),
    prior(normal(0, 0.5), class = b, coef = "sample_size"),
    prior(exponential(5), class = sd, group = "phylo"),
    prior(exponential(5), class = sd, group = "species")),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 4, cores = 4, iter = 3000, warmup = 500
)

#_______________________________________________________________________________
### 4a. Univariate predictor models

## Longevity
set.seed(666)
trend_lon <- brm(
  coef_trend ~ 1 + longevity + bodymass + sample_size + (1|gr(phylo, cov = A_trend)) + (1| species),  
  data = mam_trend, family = gaussian(),
  data2 = list(A_trend = A_trend),
  prior = c(
    prior(normal(0, 0.1), class =  Intercept),
    prior(normal(0, 0.2), class = b),
    prior(normal(0, 0.5), class = b, coef = "sample_size"),
    prior(exponential(6), class = sd, group = "phylo"),
    prior(exponential(5), class = sd, group = "species")),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 4, cores = 4, iter = 3000, warmup = 500
)

## Litter size
set.seed(666)
trend_lit <- brm(
  coef_trend ~ 1 + litter + bodymass + sample_size + (1|gr(phylo, cov = A_trend)) + (1| species),  
  data = mam_trend, family = gaussian(),
  data2 = list(A_trend = A_trend),
  prior = c(
    prior(normal(0, 0.1), class =  Intercept),
    prior(normal(0, 0.2), class = b),
    prior(normal(0, 0.5), class = b, coef = "sample_size"),
    prior(exponential(6), class = sd, group = "phylo"),
    prior(exponential(5), class = sd, group = "species")),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 4, cores = 4, iter = 3000, warmup = 500
)


## Model comparisons
trend_base <- add_criterion(trend_base, criterion = c("loo","waic"))
trend_lon <- add_criterion(trend_lon, criterion = c("loo","waic"))
trend_lit <- add_criterion(trend_lit, criterion = c("loo","waic"))

loo_compare(trend_base, trend_lon, trend_lit, criterion = "loo")


