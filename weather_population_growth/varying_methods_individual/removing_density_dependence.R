####################################################
##                                                ##
##     Global climate and population dynamics     ##
##                                                ##
##  Testing the importance of density dependence  ##
##                                                ##
##                Dec 2nd 2020                    ##
##                                                ##
####################################################

# Weather effects with and without the density dependence term to see importance

rm(list = ls())
options(width = 100)

library(tidyverse)
library(viridis)
library(patchwork)

##__________________________________________________________________________________________________
#### 1. Load data ####

# mammal data
load("../rawdata/mammal.RData")
glimpse(mammal)

# annual weather anomaly - focus on just the mean anomaly in this script at a 5km range
mam_chelsa_annual <- readRDS("data/mam_chelsa_annual.RDS") %>% 
  filter(scale == "scale_5km") %>% 
  dplyr::select(ID,year, weather_scale = scale, mean_temp_anomaly, mean_precip_anomaly)
glimpse(mam_chelsa_annual)

# Species names to merge
load("../rawdata/GBIF_species_names_mamUPDATE.RData", verbose = TRUE)

##__________________________________________________________________________________________________
#### 2. Joining data ####

# linking to weather data and species names 
mammal_weather <- mammal %>% 
  left_join(., y = mam_chelsa_annual, by = c("ID", "year")) %>% 
  left_join(., y = dplyr::select(lpd_gbif, Binomial, gbif_species = gbif.species.tree),
            by = "Binomial") %>% 
  mutate(year_s = as.numeric(scale(year)))

glimpse(mammal_weather)

##__________________________________________________________________________________________________
#### 4. Linear models for each record ####

## With and with out the abundance in each year

pgr_weather <- mammal_weather %>% 
  group_by(ID_block) %>% 
  group_modify(~{
    
    # Temperature
    mod_temp = lm(pop_growth_rate ~ mean_temp_anomaly + ln_abundance + year, data = .)
    mod_temp_noDD = lm(pop_growth_rate ~ mean_temp_anomaly + year, data = .)
    
    # Precipitation + dealing with NA values
    if(length(which(is.na(.$mean_precip_anomaly) == T)) == 0){
      mod_precip = lm(pop_growth_rate ~ mean_precip_anomaly + ln_abundance + year, data = .)
      mod_precip_noDD = lm(pop_growth_rate ~ mean_precip_anomaly + year, data = .)
      coef_precipmod = mod_precip$coefficients
      coef_precipmod_noDD = mod_precip_noDD$coefficients}
    else{coef_precipmod = rep(NA,4)
         coef_precipmod_noDD = rep(NA,3)}
    
    tibble(.[1,],
           coef_temp = mod_temp$coefficients[2],
           coef_temp_noDD = mod_temp_noDD$coefficients[2],
           coef_precip = coef_precipmod[2],
           coef_precip_noDD = coef_precipmod_noDD[2],
           coef_abun = mod_temp$coefficients[3], 
           coef_trend = mod_temp$coefficients[4],
           
           coef_abun2 = coef_precipmod[3], 
           coef_trend2 = coef_precipmod[4],
           n_obs = nrow(.))
  }) 

glimpse(pgr_weather)

##__________________________________________________________________________________________________
#### 5. Plotting differences in the coefficients ####

temp_colour <- viridis(20, option = "C")[13]
precip_colour <- viridis(20, option = "D")[10]

#____________________________________
## Temperature

ctest_temp <- cor.test(pgr_weather$coef_temp, pgr_weather$coef_temp_noDD)
temp_rho <-  round(ctest_temp$estimate, 3)
temp_p <- round(ctest_temp$p.value, 3)

temp_plot <- ggplot(pgr_weather, aes(x = coef_temp, y = coef_temp_noDD)) +
  geom_point(size = 3, alpha = 0.8, colour = temp_colour) +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "Temperature coefficient", y = "Temperature coefficient\n(no density dependence)") +
  annotate(geom = "text", label = bquote(rho == .(temp_rho) ~ ", p" < 0.001), x = -5, y = 5) +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank())

#____________________________________
## Precipitation

ctest_precip <- cor.test(pgr_weather$coef_precip, pgr_weather$coef_precip_noDD)
precip_rho <-  round(ctest_precip$estimate, 3)
precip_p <- round(ctest_precip$p.value, 3)

precip_plot <- ggplot(pgr_weather, aes(x = coef_precip, y = coef_precip_noDD)) +
  geom_point(size = 3, alpha = 0.8, colour = precip_colour) +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "Precipitation coefficient", y = "Precipitation coefficient\n(no density dependence)") +
  annotate(geom = "text", label = bquote(rho == .(precip_rho) ~ ", p" < 0.001), x = -0, y = 4) +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank())

ggsave(temp_plot + precip_plot,
       filename = "plots/weather_pop_growth/density_dependence_term_test.jpeg",
       width = 17, height = 11, units = "cm", dpi = 600)


  






