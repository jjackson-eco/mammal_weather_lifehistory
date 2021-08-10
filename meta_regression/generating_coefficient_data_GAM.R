#############################################################
##                                                         ##
##          Global climate and population dynamics         ##
##                                                         ##
##               GAM coefficient data script               ##
##                                                         ##
##                   Aug 10th 2020                         ##
##                                                         ##
#############################################################

## Generating the data to be used in predictions and plots from GAM ARMA models
## for both temperature and precipitation. Includes phylogenetic distance matrix.
## This step is repeated in model selection scripts, but here as a place holder.

rm(list = ls())
options(width = 100)

library(tidyverse)
library(ape)

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

mam_coef <- mnanom_5km_GAM %>% 
  ungroup() %>%  ## <- The groups ruin the z-transformations
  left_join(x = ., y = dplyr::select(lpd_gbif, Binomial, gbif.species.id), 
            by = "Binomial") %>% 
  left_join(x = ., y = dplyr::select(lifehistory, c(3,6:9)),
            by = "gbif.species.id") %>% 
  mutate(species = gsub(" ", "_", gbif_species),
         phylo = species,
         # raw coefficients
         coef_temp_raw = coef_temp,
         coef_precip_raw = coef_precip,
         # z transformed coefficients
         coef_temp = as.numeric(scale(coef_temp)),
         coef_precip = as.numeric(scale(coef_precip)),
         # absolute values of z transformed coefficients
         abs_temp = abs(coef_temp),
         abs_precip = abs(coef_precip),  ## <----- Precipitation studies have some NA values
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
                biome, Latitude, Longitude, lat, iucn = IUCNstatus, litter,
                longevity, bodymass, coef_temp_raw, coef_precip_raw,
                coef_temp, coef_precip, abs_temp, 
                abs_precip, log_abs_temp, log_abs_precip) %>% 
  drop_na(litter, longevity, bodymass)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3. Phylogenetic covariance matrix ####

# Trim tree to right data
mamMCC_coef <- keep.tip(mamMCC_pruned, mam_coef$phylo) 

# Covariance matrix - Brownian motion model
A_coef <- ape::vcv.phylo(mamMCC_coef)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 4. Save the data ####

save(mam_coef, mamMCC_coef, A_coef, file = "data/mammal_analysis_data_GAM.RData")


