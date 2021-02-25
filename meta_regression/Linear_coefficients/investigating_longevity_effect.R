####################################################
##                                                ##
##     Global climate and population dynamics     ##
##                                                ##
##        Investigating the longevity effect      ##
##                                                ##
##                                                ##
##                Nov 19th 2020                   ##
##                                                ##
####################################################

## Investigating the result that short-lived species have 
## negative temperature responses and positive precipitation responses

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
glimpse(mnanom_5km)

## Phylogenetic
load("../Phlogenetics/mam_lpd_pruned.RData", verbose = TRUE)
str(mamMCC_pruned)

## Life history data
load("data/lifehistory.RData", verbose = TRUE)

## Species names to link data
load("../rawdata/GBIF_species_names_mamUPDATE.RData", verbose = TRUE)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 2. Mammal coefficient data ####

mam_temp <- mnanom_5km %>% 
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
         sample_size = as.numeric(scale(log(n_obs))),
         # longevity as a factor
         longevity_factor = if_else(longevity > mean(longevity,na.rm = TRUE), "high", "low")) %>% 
  dplyr::select(id = ID, id_block = ID_block, n_obs, sample_size, 
                order = Order, species, phylo, 
                biome, lat, iucn = IUCNstatus, litter,
                longevity, longevity_factor, bodymass, coef_temp, 
                coef_precip) %>% 
  drop_na(litter, longevity, bodymass)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3. Phylogenetic covariance matrix ####

# Trim tree to right data
mamMCC_temp <- keep.tip(mamMCC_pruned, mam_temp$phylo) 

# Covariance matrix - Brownian motion model
A_temp <- ape::vcv.phylo(mamMCC_temp)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 4. Binning longevity and summarising data ####

ln_dat <- mam_temp %>% 
  mutate(ln_bin = as.numeric(as.character(cut(longevity, 
                                              breaks = c(seq(-2.1, 2, by = 0.5), 2), 
                                              labels = seq(-2.1, 2, by = 0.5))))) %>% 
  group_by(ln_bin) %>% 
  summarise(mntemp = mean(coef_temp),
            vartemp = var(coef_temp),
            mnprecip = mean(coef_precip, na.rm = T),
            varprecip = var(coef_precip, na.rm = T)) 

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 5. Plots and correlations ####

temp_colour <- viridis(20, option = "C")[13]
precip_colour <- viridis(20, option = "D")[10]

my_theme <- theme_bw(base_size = 11) +
  theme(panel.grid = element_blank())

##_____________________________________________________________________________
## 5a. Correlations
mntemp_cor <- cor.test(ln_dat$ln_bin, ln_dat$mntemp)
vartemp_cor <- cor.test(ln_dat$ln_bin, ln_dat$vartemp)

mnprecip_cor <- cor.test(ln_dat$ln_bin, ln_dat$mnprecip)
varprecip_cor <- cor.test(ln_dat$ln_bin, ln_dat$varprecip)

##_____________________________________________________________________________
## 5b. Temperature plots
mntemp_plot <- ggplot(ln_dat, aes(x = ln_bin, y = mntemp)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 4, colour = temp_colour) +
  geom_smooth(method = "lm", se = F, colour = "black") +
  annotate("text", x = 0.4, y = -0.35, 
           label = paste0("r = ", round(mntemp_cor$estimate, 2),
                          ", p = ", round(mntemp_cor$p.value, 2))) +
  labs(x = "Standardised longevity (binned)", y = "Mean temperature effect") +
  my_theme
  
vartemp_plot <- ggplot(ln_dat, aes(x = ln_bin, y = vartemp)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 4, colour = temp_colour) +
  geom_smooth(method = "lm", se = F, colour = "black") +
  annotate("text", x = 0.2, y = 4, 
           label = paste0("r = ", round(vartemp_cor$estimate, 2),
                          ", p = ", round(vartemp_cor$p.value, 2))) +
  labs(x = "Standardised longevity (binned)", y = "Variance of temperature effect") +
  my_theme

##_____________________________________________________________________________
## 5c. Precipitation plots
mnprecip_plot <- ggplot(ln_dat, aes(x = ln_bin, y = mnprecip)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 4, colour = precip_colour) +
  geom_smooth(method = "lm", se = F, colour = "black") +
  annotate("text", x = -0.2, y = 0.7, 
           label = paste0("r = ", round(mnprecip_cor$estimate, 2),
                          ", p = ", round(mnprecip_cor$p.value, 2))) +
  labs(x = "Standardised longevity (binned)", y = "Mean precipitation effect") +
  my_theme

varprecip_plot <- ggplot(ln_dat, aes(x = ln_bin, y = varprecip)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 4, colour = precip_colour) +
  geom_smooth(method = "lm", se = F, colour = "black") +
  annotate("text", x = -0.2, y = 2.8, 
           label = paste0("r = ", round(varprecip_cor$estimate, 2),
                          ", p = ", round(varprecip_cor$p.value, 2))) +
  labs(x = "Standardised longevity (binned)", y = "Variance of precipitation effect") +
  my_theme

##_____________________________________________________________________________
## 5d. Save

ggsave((mntemp_plot + vartemp_plot) /
       (mnprecip_plot + varprecip_plot),
      filename = "plots/phylogenetic_regression/meta_regression_full/longevity_effect_exploration_bins.jpeg",
      width = 15, height = 15, units = "cm", dpi = 500)





