####################################################
##                                                ##
##     Global climate and population dynamics     ##
##                                                ##
##   Annual mean anomalies and population growth  ##
##                                                ##
##                June 9th 2020                   ##
##                                                ##
####################################################
rm(list = ls())
options(width = 100)

library(tidyverse)
library(ggridges)
library(viridis)
library(grid)
library(gridExtra)

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

##__________________________________________________________________________________________________
#### 2. Joining data ####

mammal_weather <- mammal %>% 
  left_join(., y = mam_chelsa_annual, by = c("ID", "year"))

##__________________________________________________________________________________________________
#### 3. Hypothesis plots across regions and taxa ####

# These are to get a sense of how the response and predictors look against one another, 
# not actually accounting for the study ID or other structure variables

# pop growth and anomaly - Species
sp_temp_pgr <- ggplot(mammal_weather, aes(x = mean_temp_anomaly, y = pop_growth_rate)) +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = 0) +
  geom_point(alpha = 0.2, colour = "firebrick", size = 2.5) + 
  scale_x_continuous(breaks = seq(-0.7, 0.7, by = 0.2)) +
  facet_wrap(~Order) +
  labs(x = "Mean annual temperature anomaly", y = "Population growth rate") +
  theme_bw(base_size = 24) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())
  
sp_precip_pgr <- ggplot(mammal_weather, aes(x = mean_precip_anomaly, y = pop_growth_rate)) +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = 0) +
  geom_point(alpha = 0.2, colour = "darkblue", size = 2.5) + 
  facet_wrap(~Order) +
  labs(x = "Mean annual precipitation anomaly", y = "Population growth rate") +
  theme_bw(base_size = 24) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())

ggsave(grid.arrange(sp_temp_pgr, sp_precip_pgr, ncol = 1),
       filename = "plots/weather_pop_growth/species_anomal_pgr.jpeg",
       width = 15, height = 20, units = "in", dpi = 400)

# pop growth and anomaly - Biome
mammal_weather <- mammal_weather %>% 
  mutate(biome_lab = case_when(
    biome == "Mediterranean forests, woodlands and scrub" ~ "Mediterranean forests,\nwoodlands and scrub",
    biome == "Temperate grasslands, savannas and shrublands" ~ "Temperate grasslands,\nsavannas and shrublands",
    biome == "Tropical and subtropical coniferous forests" ~ "Tropical and subtropical\nconiferous forests",
    biome == "Tropical and subtropical dry broadleaf forests" ~ "Tropical and subtropical\ndry broadleaf forests",
    biome == "Tropical and subtropical grasslands, savannas and shrublands" ~ "Tropical and subtropical grasslands,\nsavannas and shrublands",
    biome == "Tropical and subtropical moist broadleaf forests" ~ "Tropical and subtropical\nmoist broadleaf forests",
    TRUE ~ biome
  ))

biome_temp_pgr <- ggplot(mammal_weather, aes(x = mean_temp_anomaly, y = pop_growth_rate)) +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = 0) +
  geom_point(alpha = 0.2, colour = "firebrick", size = 2.5) + 
  scale_x_continuous(breaks = seq(-0.7, 0.7, by = 0.2)) +
  facet_wrap(~biome_lab) +
  labs(x = "Mean annual temperature anomaly", y = "Population growth rate") +
  theme_bw(base_size = 24) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 13))

biome_precip_pgr <- ggplot(mammal_weather, aes(x = mean_precip_anomaly, y = pop_growth_rate)) +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = 0) +
  geom_point(alpha = 0.2, colour = "darkblue", size = 2.5) + 
  facet_wrap(~biome_lab) +
  labs(x = "Mean annual precipitation anomaly", y = "Population growth rate") +
  theme_bw(base_size = 24) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 13))

ggsave(grid.arrange(biome_temp_pgr, biome_precip_pgr, ncol = 1),
       filename = "plots/weather_pop_growth/biome_anomal_pgr.jpeg",
       width = 15, height = 20, units = "in", dpi = 400)

#Not too much going on

##__________________________________________________________________________________________________
#### 4. Linear models for each record ####

pgr_weather <- mammal_weather %>% 
  group_by(ID_block) %>% 
  group_modify(~{
    
    # Temperature
    mod_temp = lm(pop_growth_rate ~ mean_temp_anomaly + ln_abundance + year, data = .)
    
    # Precipitation + dealing with NA values
    if(length(which(is.na(.$mean_precip_anomaly) == T)) == 0){
    mod_precip = lm(pop_growth_rate ~ mean_precip_anomaly + ln_abundance + year, data = .)
    coef_precipmod = mod_precip$coefficients}
    else{coef_precipmod = rep(NA,4)}
    
    tibble(.[1,],
           coef_temp = mod_temp$coefficients[2],
           coef_precip = coef_precipmod[2],
           coef_abun = mod_temp$coefficients[3], 
           coef_trend = mod_temp$coefficients[4],
           
           coef_abun2 = coef_precipmod[3], 
           coef_trend2 = coef_precipmod[4],
           n_obs = nrow(.))
  }) 
  
glimpse(pgr_weather)


##__________________________________________________________________________________________________
#### 5. Density ridge plots for the coefficient distributions ####

ggplot(pgr_weather, aes(x = coef_temp, y = Order, fill = Order)) +
  geom_vline(xintercept = 0) +
  geom_density_ridges(alpha = 0.8) +
  coord_cartesian(xlim = c(-3,3)) +
  scale_fill_viridis_d(guide = F, option = "C") +
  labs(x = "Temperature anomaly coefficient", y = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank())

