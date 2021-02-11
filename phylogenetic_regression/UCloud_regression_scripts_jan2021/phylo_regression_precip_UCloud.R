####################################################
##                                                ##
##     Global climate and population dynamics     ##
##                                                ##
##      Absolute precipitation effect with GAM    ##
##                                                ##
##    UCloud - brms Phylogenetic meta-regression  ##
##        life-history and spatial effects        ##
##                                                ##
##                 Jan 20th 2020                  ##
##                                                ##
####################################################

## Models regressing biome and life-history on absolute population responses
## from GAM ARMA models. Implementing the meta-regression framework for precipitation

setwd("phylogenetic_regression/")

# Loading libraries
libraries <- c('tidyverse', 'ape', 'brms', 'viridis', 'patchwork')
for(i in libraries){
  if(! i %in% installed.packages()) lapply(i, install.packages)
  lapply(libraries, require, character.only = TRUE)
}

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 1. Load data ####

## Weather effects
load("data/mnanom_5km_GAM.RData", verbose = TRUE)
glimpse(mnanom_5km_GAM)

## Phylogenetic
load("data/mam_lpd_pruned.RData", verbose = TRUE)
str(mamMCC_pruned)

## Life history data
load("data/lifehistory.RData", verbose = TRUE)

## Species names to link data
load("data/GBIF_species_names_mamUPDATE.RData", verbose = TRUE)

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
  drop_na(litter, longevity, bodymass, abs_precip)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3. Phylogenetic covariance matrix ####

# Trim tree to right data
mamMCC_precip <- keep.tip(mamMCC_pruned, mam_precip$phylo) 

# Covariance matrix - Brownian motion model
A_precip <- ape::vcv.phylo(mamMCC_precip)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 4. Precipitation models ####

# Base model chain test
plot(density(rgamma(1000,shape = 2, scale = 0.6)))
plot(density(rexp(1000, rate = 10)))

set.seed(666)
precip_base_test <- brm(
  abs_precip ~ 1 + sample_size + (1|gr(phylo, cov = A_precip)) + (1| species),
  data = mam_precip,
  family = Gamma(link = "log"),
  data2 = list(A_precip = A_precip),
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
precip_base <- brm(
  abs_precip ~ 1 + sample_size + (1|gr(phylo, cov = A_precip)) + (1| species),  
  data = mam_precip, 
  family = Gamma(link = "log"), 
  data2 = list(A_precip = A_precip),
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
precip_longevity <- brm(
  abs_precip ~ 1 + longevity + sample_size + (1|gr(phylo, cov = A_precip)) + (1| species),  
  data = mam_precip, 
  family = Gamma(link = "log"), 
  data2 = list(A_precip = A_precip),
  prior = c(
    prior(normal(0, 0.3), class =  Intercept),
    prior(normal(0, 0.3), class = b, coef = "longevity"),
    prior(normal(0, 0.3), class = b, coef = "sample_size"),
    prior(exponential(11), class = sd, group = "phylo"),
    prior(exponential(2), class = sd, group = "species"),
    prior(gamma(2,0.7), class = shape)),
  control = list(adapt_delta = 0.99,
                 max_treedepth = 17),
  chains = 3, cores = 3, iter = 4000, warmup = 2000
)

## Litter size
set.seed(666)
precip_litter <- brm(
  abs_precip ~ 1 + litter + sample_size + (1|gr(phylo, cov = A_precip)) + (1| species),  
  data = mam_precip, 
  family = Gamma(link = "log"), 
  data2 = list(A_precip = A_precip),
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
precip_bodymass <- brm(
  abs_precip ~ 1 + bodymass + sample_size + (1|gr(phylo, cov = A_precip)) + (1| species),  
  data = mam_precip, 
  family = Gamma(link = "log"), 
  data2 = list(A_precip = A_precip),
  prior = c(
    prior(normal(0, 0.8), class =  Intercept),
    prior(normal(0, 0.5), class = b, coef = "bodymass"),
    prior(normal(0, 0.8), class = b, coef = "sample_size"),
    prior(exponential(11), class = sd, group = "phylo"),
    prior(exponential(2), class = sd, group = "species"),
    prior(gamma(2,0.5), class = shape)),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 4000, warmup = 2000
)

## Biome
set.seed(666)
precip_biome <- brm(
  abs_precip ~ 1 + biome + sample_size + (1|gr(phylo, cov = A_precip)) + (1| species),  
  data = mam_precip, 
  family = Gamma(link = "log"), 
  data2 = list(A_precip = A_precip),
  prior = c(
    prior(normal(0, 0.8), class =  Intercept),
    prior(normal(0, 0.3), class = b),
    prior(normal(0, 0.8), class = b, coef = "sample_size"),
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
precip_base <- add_criterion(precip_base, criterion = c("loo","waic"))
precip_longevity <- add_criterion(precip_longevity, criterion = c("loo","waic"))
precip_litter <- add_criterion(precip_litter, criterion = c("loo","waic"))
precip_bodymass <- add_criterion(precip_bodymass, criterion = c("loo","waic"))
precip_biome <- add_criterion(precip_biome, criterion = c("loo","waic"))

univar_modcomp <- as.data.frame(loo_compare(precip_base, precip_longevity, precip_bodymass, 
                                            precip_litter, precip_biome, criterion = "loo"))

univar_modcomp
save(univar_modcomp, file = "output/precipitation_univariate_model_comparisons.RData")

#_______________________________________________________________________________
### 4c. Life-history models full

lh_priors <- c(
  prior(normal(0, 0.3), class =  Intercept),
  prior(normal(0, 0.3), class = b),
  prior(normal(0, 0.3), class = b, coef = "sample_size"),
  prior(exponential(11), class = sd, group = "phylo"),
  prior(exponential(2), class = sd, group = "species"),
  prior(gamma(2,0.6), class = shape))

## Life-history full
set.seed(666)
precip_lh <- brm(
  abs_precip ~ 1 + longevity + bodymass + litter +
    longevity:bodymass + litter:bodymass + litter:longevity + 
    sample_size + (1|gr(phylo, cov = A_precip)) + (1| species),  
  data = mam_precip, 
  family = Gamma(link = "log"), 
  data2 = list(A_precip = A_precip),
  prior = lh_priors,
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 4000, warmup = 2000
)

## Life-history univariate only
set.seed(666)
precip_lh_uni <- brm(
  abs_precip ~ 1 + longevity + bodymass + litter +
    sample_size + (1|gr(phylo, cov = A_precip)) + (1| species),  
  data = mam_precip, 
  family = Gamma(link = "log"), 
  data2 = list(A_precip = A_precip),
  prior = lh_priors,
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 4000, warmup = 2000
)

## Longevity and Body Mass full
set.seed(666)
precip_lonbod <- brm(
  abs_precip ~ 1 + longevity*bodymass + sample_size + (1|gr(phylo, cov = A_precip)) + (1| species),  
  data = mam_precip, 
  family = Gamma(link = "log"), 
  data2 = list(A_precip = A_precip),
  prior = lh_priors,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 4000, warmup = 2000
)

## Litter and Body Mass full
set.seed(666)
precip_litbod <- brm(
  abs_precip ~ 1 + litter*bodymass + sample_size + (1|gr(phylo, cov = A_precip)) + (1| species),  
  data = mam_precip, 
  family = Gamma(link = "log"), 
  data2 = list(A_precip = A_precip),
  prior = lh_priors,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 4000, warmup = 2000
)

## Body mass correction - additive only
set.seed(666)
precip_lonbod_simple <- brm(
  abs_precip ~ 1 + longevity + bodymass + sample_size + (1|gr(phylo, cov = A_precip)) + (1| species),  
  data = mam_precip, 
  family = Gamma(link = "log"), 
  data2 = list(A_precip = A_precip),
  prior = lh_priors,
  control = list(adapt_delta = 0.98,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 4000, warmup = 2000
)

set.seed(666)
precip_litbod_simple <- brm(
  abs_precip ~ 1 + litter + bodymass + sample_size + (1|gr(phylo, cov = A_precip)) + (1| species),  
  data = mam_precip, 
  family = Gamma(link = "log"), 
  data2 = list(A_precip = A_precip),
  prior = lh_priors,
  control = list(adapt_delta = 0.98,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 4000, warmup = 2000
)

#_______________________________________________________________________________
### 4b. Full model comparisons

precip_lh <- add_criterion(precip_lh, criterion = c("loo","waic"))
precip_lh_uni <- add_criterion(precip_lh_uni, criterion = c("loo","waic"))
precip_lonbod <- add_criterion(precip_lonbod, criterion = c("loo","waic"))
precip_litbod<- add_criterion(precip_litbod, criterion = c("loo","waic"))
precip_lonbod_simple <- add_criterion(precip_lonbod_simple, criterion = c("loo","waic"))
precip_litbod_simple <- add_criterion(precip_litbod_simple, criterion = c("loo","waic"))

mod_comp <- as.data.frame(loo_compare(precip_base, precip_longevity, precip_bodymass, 
                                      precip_litter, precip_biome, precip_lonbod_simple,
                                      precip_litbod_simple, precip_lh_uni, precip_lonbod, precip_litbod,
                                      precip_lh, criterion = "loo"))

save(mod_comp, file = "output/precipitation_model_comparisons.RData")



## Best model - litter alone is the most effective, but also evidence for longevity


