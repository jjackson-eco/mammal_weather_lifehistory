####################################################
##                                                ##
##     Global climate and population dynamics     ##
##                                                ##
##        Myotis emarginatus investigation        ##
##                                                ##
##                Oct 20th 2020                   ##
##                                                ##
####################################################

## Investigating the odd abundance patterns from Geoffroy's bat Myotis emarginatus

rm(list = ls())
options(width = 100)

library(tidyverse)
library(ape)
library(brms)
library(nlme)
library(patchwork)
library(ggridges)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 1. Load data ####

## Weather effects
load("data/pgr_weather/mnanom_5km.RData", verbose = TRUE)

## Phylogenetic
load("../Phlogenetics/mam_lpd_pruned.RData", verbose = TRUE)

## Mammal abundance data
load("../rawdata/mam_IDblocks.RData")

## Weather anomaly data 
mam_chelsa_annual <- readRDS("data/mam_chelsa_annual.RDS") %>% 
  filter(scale == "scale_5km") %>% 
  dplyr::select(ID,year, weather_scale = scale, mean_temp_anomaly, mean_precip_anomaly)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 2. Geoffroy's bat data ####

## Raw abundance
My_em <- mam_IDblocks %>% 
  ungroup() %>% 
  dplyr::select(-c(block,Binomial,Class,Genus,Species,Subspecies)) %>% 
  left_join(x = ., 
            y = dplyr::select(LPD_tree_update, ID, spp = checked_speciesname),
            by = "ID") %>% 
  filter(spp == "Myotis emarginatus") %>% 
  left_join(x = ., 
            y = dplyr::select(mam_chelsa_annual, -weather_scale),
            by = c("ID", "year"))

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3. Raw data plots + summaries ####

## Data features
My_em %>% 
  group_by(ID) %>% 
  summarise(ID_block_num = n_distinct(ID_block),
            realm = realm[1], biome = biome[1], Latitude = Latitude[1],
            Longitude = Longitude[1], Specific_location = Specific_location[1],
            abundance_measure = abundance_measure[1],
            .groups = "drop") %>% 
  flextable()

# The crazy study is a full population count

## Abundance timeseries
My_em %>% 
  ggplot(aes(x = year, y = ln_abundance, colour = factor(ID), group = ID)) +
  geom_line() +
  geom_point(size = 2) +
  scale_color_manual(name = "Study ID", values = c("firebrick", "cornflowerblue")) +
  labs(x  = "Year", y = "ln Abundance", title = "Myotis emarginatus - Geoffroy's bat") +
  theme_bw(base_size = 14) +
  ggsave(filename = "plots/meta_regression/repeated_obs_spatial/Myotis_emarginatus.jpeg",
         width = 5, height = 4.5, units = "in", dpi = 400)







