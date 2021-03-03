#############################################################
##                                                         ##
##          Global climate and population dynamics         ##
##                                                         ##
##       Absolute weather effect with GAM - Gamma          ##
##                                                         ##
##         Posterior plots of Life-history effects         ##
##                                                         ##
##                   Dec 27th 2020                         ##
##                                                         ##
#############################################################

## Predicting and plotting from the posterior of Gamma meta-regression models.
## Absolute temperature and precipitation coefficients from GAM-ARMA models 

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

#_______________
## Analysis data
load("data/mammal_analysis_data_GAM.RData", verbose = TRUE)

#_____________
## Best models
temp_lh_uni <- readRDS("results/UCloud_gamma_models/temp_lh_uni_gamma.rds")
precip_lh_uni <- readRDS("results/UCloud_gamma_models/precip_lh_uni_gamma.RDS")

#______________________
## Model selection data
load("results/UCloud_gamma_models/temperature_model_comparisons.RData")
temp_modcomp <- mod_comp
load("results/UCloud_gamma_models/precipitation_model_comparisons.RData")
precip_modcomp <- mod_comp
rm(mod_comp)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3. Posterior predictions ####

#_______________________________________________________________________________
## 3a. Longevity

# prediction data
preddat_lon <- expand_grid(longevity = seq(-2.1,2, length.out = 50),
                           litter = mean(mam_coef$litter),
                           bodymass = mean(mam_coef$bodymass),
                           sample_size = mean(mam_coef$sample_size),
                           phylo = unique(mam_coef$phylo)) %>% 
  mutate(species = phylo)


# standcode to inspect how to back transform using the link function
brms::stancode(temp_lh_uni) #  mu[n] = shape * exp(-(mu[n]))

# posterior predictions
temp_pred_lon <-  brms::posterior_predict(temp_lh_uni, newdata = preddat_lon, 
                                          type = "response") 
precip_pred_lon <-  brms::posterior_predict(precip_lh_uni, 
                                            newdata = preddat_lon, 
                                            # account for NAs in the precipitation coefficients
                                            allow_new_levels = TRUE,
                                            type = "response") 

# summary data for each longevity value
posterior_summary_longevity <- bind_rows(lapply(unique(preddat_lon$longevity), function(x){
  
  cpos = which(preddat_lon$longevity == x) # vector of positions corresponding to where the values are equal to x
  
  # posterior mean - Exponentiate to get back on response scale
  post_mn_temp = mean(temp_pred_lon[,cpos])
  post_mn_precip = mean(precip_pred_lon[,cpos])
  
  # prediction intervals - 80% from rethinking package
  cPI_temp = PI(temp_pred_lon[,cpos], prob = 0.8)
  cQuant_temp = quantile(temp_pred_lon[,cpos], c(0.025, 0.975))
  
  cPI_precip = PI(precip_pred_lon[,cpos], prob = 0.8)
  cQuant_precip = quantile(precip_pred_lon[,cpos], c(0.025, 0.975))
  
  # return data
  return(tibble(longevity = x, 
                post_mn_temp = post_mn_temp, 
                post_mn_precip = post_mn_precip,
                lwrPI_temp = cPI_temp[1], uprPI_temp = cPI_temp[2], 
                lwr_temp = cQuant_temp[1], upr_temp = cQuant_temp[2],
                lwrPI_precip = cPI_precip[1], uprPI_precip = cPI_precip[2], 
                lwr_precip = cQuant_precip[1], upr_precip = cQuant_precip[2]))
}))

#_______________________________________________________________________________
## 3b. Litter size

## precipitation species names

# prediction data
preddat_lit <- expand_grid(litter = seq(-1.2,3, length.out = 50),
                           longevity = mean(mam_coef$longevity),
                           bodymass = mean(mam_coef$bodymass),
                           sample_size = mean(mam_coef$sample_size),
                           phylo = unique(mam_coef$phylo)) %>% 
  mutate(species = phylo)

# posterior predictions
temp_pred_lit <-  brms::posterior_predict(temp_lh_uni, newdata = preddat_lit, 
                                          type = "response") 
precip_pred_lit <-  brms::posterior_predict(precip_lh_uni, 
                                            newdata = preddat_lit, 
                                            allow_new_levels = TRUE,
                                            type = "response") 

# summary data for each litter size value
posterior_summary_litter <- bind_rows(lapply(unique(preddat_lit$litter), function(x){
  
  cpos = which(preddat_lit$litter == x) # vector of positions corresponding to where the values are equal to x
  
  # posterior mean - Exponentiate to get back on response scale
  post_mn_temp = mean(temp_pred_lit[,cpos])
  post_mn_precip = mean(precip_pred_lit[,cpos])
  
  # prediction intervals - 80% from rethinking package
  cPI_temp = PI(temp_pred_lit[,cpos], prob = 0.8)
  cQuant_temp = quantile(temp_pred_lit[,cpos], c(0.025, 0.975))
  
  cPI_precip = PI(precip_pred_lit[,cpos], prob = 0.8)
  cQuant_precip = quantile(precip_pred_lit[,cpos], c(0.025, 0.975))
  
  # return data
  return(tibble(litter = x, 
                post_mn_temp = post_mn_temp, 
                post_mn_precip = post_mn_precip,
                lwrPI_temp = cPI_temp[1], uprPI_temp = cPI_temp[2], 
                lwr_temp = cQuant_temp[1], upr_temp = cQuant_temp[2],
                lwrPI_precip = cPI_precip[1], uprPI_precip = cPI_precip[2], 
                lwr_precip = cQuant_precip[1], upr_precip = cQuant_precip[2]))
}))

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 4. Life-history posterior plots ####

## Plot colours
temp_colour <- viridis(20, option = "C")[14]
precip_colour <- viridis(20, option = "D")[14]

#_______________________________________________________________________________
## 4a. Creating binned data

longevity_binned_data <- mam_coef %>% 
  mutate(lon_bin = as.numeric(as.character(cut(longevity, 
                                               breaks = c(seq(-2.1, 2, by = 0.2), 2), 
                                               labels = seq(-2.1, 2, by = 0.2))))) %>% 
  group_by(lon_bin) %>% 
  summarise(mntemp = mean(abs_temp),
            setemp = sd(abs_temp)/sqrt(n()),
            mnprecip = mean(abs_precip, na.rm = T),
            seprecip = sd(abs_precip, na.rm = T)/sqrt(n())) 

litter_binned_data <- mam_coef %>% 
  mutate(lit_bin = 
           as.numeric(as.character(cut(litter, 
                                       breaks = seq(-1.2, 2.8, by = 0.2), 
                                       labels = seq(-1.2, 2.6, by = 0.2))))) %>% 
  group_by(lit_bin) %>% 
  summarise(mntemp = mean(abs_temp),
            setemp = sd(abs_temp)/sqrt(n()),
            mnprecip = mean(abs_precip, na.rm = T),
            seprecip = sd(abs_precip, na.rm = T)/sqrt(n())) 

#_______________________________________________________________________________
## 4b. Raw data plots

lon_temp <- ggplot(mam_coef, aes(x = longevity, y = abs_temp)) +
  geom_point(aes(size = n_obs), colour = temp_colour, alpha = 0.5) +
  geom_smooth(data = posterior_summary_longevity,
              aes(y = post_mn_temp, ymax = uprPI_temp, ymin = lwrPI_temp),
              colour = "black", fill = temp_colour, alpha = 0.3,
              stat = "identity") +
  coord_cartesian(ylim = c(0,4)) +
  scale_size_continuous(range = c(2,7), guide = F) +
  labs(x = "Standardised longevity", y = "|Temperature effect on abundance|",
       tag = "a)") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank()) 

lit_temp <- ggplot(mam_coef, aes(x = litter, y = abs_temp)) +
  geom_point(aes(size = n_obs), colour = temp_colour, alpha = 0.5) +
  geom_smooth(data = posterior_summary_litter,
              aes(y = post_mn_temp, ymax = uprPI_temp, ymin = lwrPI_temp),
              colour = "black", fill = temp_colour, alpha = 0.3,
              stat = "identity") +
  coord_cartesian(ylim = c(0,4)) +
  scale_size_continuous(range = c(2,7), guide = F) +
  labs(x = "Standardised litter size", y = "|Temperature effect on abundance|",
       tag = "b)") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank()) 

lon_precip <- ggplot(mam_coef, aes(x = longevity, y = abs_precip)) +
  geom_point(aes(size = n_obs), colour = precip_colour, alpha = 0.5) +
  geom_smooth(data = posterior_summary_longevity,
              aes(y = post_mn_precip, ymax = uprPI_precip, ymin = lwrPI_precip),
              colour = "black", fill = precip_colour, alpha = 0.3,
              stat = "identity") +
  coord_cartesian(ylim = c(0,4)) +
  scale_size_continuous(range = c(2,7), guide = F) +
  labs(x = "Standardised longevity", y = "|Precipitation effect on abundance|",
       tag = "c)") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank()) 

lit_precip <- ggplot(mam_coef, aes(x = litter, y = abs_precip)) +
  geom_point(aes(size = n_obs), colour = precip_colour, alpha = 0.5) +
  geom_smooth(data = posterior_summary_litter,
              aes(y = post_mn_precip, ymax = uprPI_precip, ymin = lwrPI_precip),
              colour = "black", fill = precip_colour, alpha = 0.3,
              stat = "identity") +
  coord_cartesian(ylim = c(0,4)) +
  scale_size_continuous(range = c(2,7), guide = F) +
  labs(x = "Standardised litter size", y = "|Precipitation effect on abundance|",
       tag = "d)") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank()) 

# the final figure
(lon_temp | lit_temp)/
  (lon_precip | lit_precip)

ggsave((lon_temp | lit_temp)/
         (lon_precip | lit_precip),
       file = "plots/meta_regression/meta_regression_full/posterior_lifehistory_Gamma.jpeg",
       width = 20, height = 20, units = "cm", dpi = 600)

#_______________________________________________________________________________
## 4c. Binned data plots

lon_temp_bin <- ggplot(longevity_binned_data, aes(x = lon_bin, y = mntemp)) +
  geom_errorbar(aes(ymax = mntemp + setemp, ymin = mntemp - setemp),
                width = 0.03, colour = temp_colour) +
  geom_point(size = 4, colour = temp_colour) +
  geom_smooth(data = posterior_summary_longevity,
              aes(x = longevity, y = post_mn_temp, 
                  ymax = uprPI_temp, ymin = lwrPI_temp),
              colour = "black", fill = temp_colour, alpha = 0.3,
              stat = "identity") +
  coord_cartesian(ylim = c(0,3)) +
  scale_size_continuous(range = c(2,8), guide = F) +
  labs(x = "Standardised longevity", y = "|Temperature effect on abundance|",
       tag = "a)") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank()) 

lit_temp_bin <- ggplot(litter_binned_data, aes(x = lit_bin, y = mntemp)) +
  geom_errorbar(aes(ymax = mntemp + setemp, ymin = mntemp - setemp),
                width = 0.03, colour = temp_colour) +
  geom_point(size = 4, colour = temp_colour) +
  geom_smooth(data = posterior_summary_litter,
              aes(x = litter, y = post_mn_temp, 
                  ymax = uprPI_temp, ymin = lwrPI_temp),
              colour = "black", fill = temp_colour, alpha = 0.3,
              stat = "identity") +
  scale_size_continuous(range = c(2,8), guide = F) +
  labs(x = "Standardised litter size", y = "|Temperature effect on abundance|",
       tag = "b)") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank()) 

lon_precip_bin <- ggplot(longevity_binned_data, aes(x = lon_bin, y = mnprecip)) +
  geom_errorbar(aes(ymax = mnprecip + seprecip, ymin = mnprecip - seprecip),
                width = 0.03, colour = precip_colour) +
  geom_point(size = 4, colour = precip_colour) +
  geom_smooth(data = posterior_summary_longevity,
              aes(x = longevity, y = post_mn_precip, 
                  ymax = uprPI_precip, ymin = lwrPI_precip),
              colour = "black", fill = precip_colour, alpha = 0.3,
              stat = "identity") +
  coord_cartesian(ylim = c(0,4)) +
  scale_size_continuous(range = c(2,8), guide = F) +
  labs(x = "Standardised longevity", y = "|Precipitation effect on abundance|",
       tag = "a)") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank()) 

lit_precip_bin <- ggplot(litter_binned_data, aes(x = lit_bin, y = mnprecip)) +
  geom_errorbar(aes(ymax = mnprecip + seprecip, ymin = mnprecip - seprecip),
                width = 0.03, colour = precip_colour) +
  geom_point(size = 4, colour = precip_colour) +
  geom_smooth(data = posterior_summary_litter,
              aes(x = litter, y = post_mn_precip, 
                  ymax = uprPI_precip, ymin = lwrPI_precip),
              colour = "black", fill = precip_colour, alpha = 0.3,
              stat = "identity") +
  scale_size_continuous(range = c(2,8), guide = F) +
  labs(x = "Standardised litter size", y = "|Precipitation effect on abundance|",
       tag = "b)") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank()) 

# the final figure
(lon_temp_bin | lit_temp_bin)/
  (lon_precip_bin | lit_precip_bin)

ggsave((lon_temp_bin | lit_temp_bin)/
         (lon_precip_bin | lit_precip_bin),
       file = "plots/meta_regression/meta_regression_full/posterior_lifehistory_binned_Gamma.jpeg",
       width = 20, height = 20, units = "cm", dpi = 600)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 5. Saving data for manuscript figure ####

save(posterior_summary_longevity, posterior_summary_litter,
     longevity_binned_data, litter_binned_data,
     file = "results/UCloud_gamma_models/posterior_summary_Gamma.RData")


