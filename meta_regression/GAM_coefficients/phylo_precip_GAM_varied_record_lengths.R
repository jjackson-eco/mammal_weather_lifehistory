#################################################################
##                                                             ##
##           Global climate and population dynamics            ##
##                                                             ##
## Precipitation with GAM coefficients - Varied Record Lengths ##
##                                                             ##
##             brms Phylogenetic meta-regression               ##
##                   and spatial effects                       ##
##                                                             ##
##                      Feb 2nd 2022                           ##
##                                                             ##
#################################################################

## Investigating the effects of biome on population responses
## from GAM ARMA models. Implementing the meta-regression framework for precipitation

rm(list = ls())
options(width = 100)

## Change the .libPaths and R_LIBS_USER to the right thing if you're on a uni computer
if(Sys.info()["nodename"] == "BIO-W-LT-000083" |
   Sys.info()["nodename"] == "BIO-W-DT-02108") {
  .libPaths("C:/Users/zool2541/R-4.1.1/library/")
  .libPaths("!\\\\zoo-suitcase/home$/zool2541/My Documents/R/win-library/4.1")}

## Packages
library(tidyverse)
library(tidybayes)
library(ape)
library(brms)
library(nlme)
library(patchwork)
library(ggridges)
library(ggdist)
library(viridis)
library(flextable)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 1. Load data ####

## Weather effects
load("data/pgr_weather/mnanom_5km_GAM_varied_record_length.RData", verbose = TRUE)

## Phylogenetic
load("../Phlogenetics/mam_lpd_pruned.RData", verbose = TRUE)
str(mamMCC_pruned)

## Life history data
load("data/lifehistory.RData", verbose = TRUE)

## Species names to link data
load("../rawdata/GBIF_species_names_mamUPDATE.RData", verbose = TRUE)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 2. Mammal coefficient data ####

# Short records
mam_precip_5yr <- mnanom_5km_GAM_5yr %>% 
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
  drop_na(litter, longevity, bodymass, coef_precip) %>% 
  filter(phylo != "Damaliscus_korrigum")

# Long records
mam_precip_20yr <- mnanom_5km_GAM_20yr %>% 
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
  drop_na(litter, longevity, bodymass, coef_precip)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3. Phylogenetic covariance matrices ####

## Short Records
# Trim tree to right data
mamMCC_precip_5yr <- keep.tip(mamMCC_pruned, mam_precip_5yr$phylo) 

# Covariance matrix - Brownian motion model
A_precip_5yr <- ape::vcv.phylo(mamMCC_precip_5yr)

## Long Records
# Trim tree to right data
mamMCC_precip_20yr <- keep.tip(mamMCC_pruned, mam_precip_20yr$phylo) 

# Covariance matrix - Brownian motion model
A_precip_20yr <- ape::vcv.phylo(mamMCC_precip_20yr)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 4. Precipitation models - Short records ####

## Base model
set.seed(666)
precip_base_5yr <- brm(
  coef_precip ~ 1 + sample_size + (1|gr(phylo, cov = A_precip_5yr)) + (1| species),  
  data = mam_precip_5yr, family = gaussian(),
  data2 = list(A_precip_5yr = A_precip_5yr),
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
precip_biome_5yr <- brm(
  coef_precip ~ 1 + biome + sample_size + (1|gr(phylo, cov = A_precip_5yr)) + (1| species),  
  data = mam_precip_5yr, family = gaussian(),
  data2 = list(A_precip_5yr = A_precip_5yr),
  prior = c(
    prior(normal(0, 0.2), class =  Intercept),
    prior(normal(0, 0.1), class = b),
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
precip_base_5yr <- add_criterion(precip_base_5yr, criterion = c("loo","waic"))
precip_biome_5yr <- add_criterion(precip_biome_5yr, criterion = c("loo","waic"))

mod_comp_precip_5yr <- as.data.frame(loo_compare(precip_base_5yr, 
                                                 precip_biome_5yr, criterion = "loo"))

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 5. Precipitation models - Long records ####

## Base model
set.seed(666)
precip_base_20yr <- brm(
  coef_precip ~ 1 + sample_size + (1|gr(phylo, cov = A_precip_20yr)) + (1| species),  
  data = mam_precip_20yr, family = gaussian(),
  data2 = list(A_precip_20yr = A_precip_20yr),
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
precip_biome_20yr <- brm(
  coef_precip ~ 1 + biome + sample_size + (1|gr(phylo, cov = A_precip_20yr)) + (1| species),  
  data = mam_precip_20yr, family = gaussian(),
  data2 = list(A_precip_20yr = A_precip_20yr),
  prior = c(
    prior(normal(0, 0.2), class =  Intercept),
    prior(normal(0, 0.1), class = b),
    prior(normal(0, 0.3), class = b, coef = "sample_size"),
    prior(exponential(8), class = sd, group = "phylo"),
    prior(exponential(8), class = sd, group = "species")),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 4000, warmup = 2000
)

#_______________________________________________________________________________
### 5b. Model comparisons

## Model comparisons
precip_base_20yr <- add_criterion(precip_base_20yr, criterion = c("loo","waic"))
precip_biome_20yr <- add_criterion(precip_biome_20yr, criterion = c("loo","waic"))

mod_comp_precip_20yr <- as.data.frame(loo_compare(precip_base_20yr, 
                                                 precip_biome_20yr, criterion = "loo"))

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 6. Model plots ####

precip_colour <- "#287f79"

## 6a. 5 years 
precip_post_5yr <- precip_base_5yr %>%
  gather_draws(`b_Intercept|sd_.*|b_sample_size|sigma`, regex = TRUE) %>% #tidybayes
  ungroup() %>% 
  ggplot(aes(y = .variable, x = .value)) + 
  stat_halfeye(show.legend = FALSE, fill = precip_colour) +
  geom_vline(xintercept = 0, size = 0.6) +
  scale_y_discrete(labels = c(expression(paste("Global intercept ", bar(alpha))),
                              expression(paste("Sample size ", beta[N])),
                              expression(paste("Phylogenetic covariance ", sigma[PHYLO])),
                              expression(paste("Species level variance ", sigma[SPECIES])),
                              "Population-level variance")) +
  labs(x = "Posterior estimate", y = NULL, title = "Records with > 5 years of data") +
  theme_ridges(center_axis_labels = TRUE, grid = T, line_size = 0.3) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))

precip_post_20yr <- precip_base_20yr %>%
  gather_draws(`b_Intercept|sd_.*|b_sample_size|sigma`, regex = TRUE) %>% #tidybayes
  ungroup() %>% 
  ggplot(aes(y = .variable, x = .value)) + 
  stat_halfeye(show.legend = FALSE, fill = precip_colour) +
  geom_vline(xintercept = 0, size = 0.6) +
  scale_y_discrete(labels = NULL) +
  labs(x = "Posterior estimate", y = NULL, title = "Records with > 20 years of data") +
  theme_ridges(center_axis_labels = TRUE, grid = T, line_size = 0.3) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))

precip_post_5yr + precip_post_20yr
ggsave(precip_post_5yr + precip_post_20yr,
       filename = "plots/manuscript_figures/Supplementary figures/precip_varyinglength_posterior.jpeg",
       width = 30, height = 20, units = "cm", dpi = 1000)


