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

# Species names to merge
load("../rawdata/GBIF_species_names_mamUPDATE.RData", verbose = TRUE)

##__________________________________________________________________________________________________
#### 2. Joining data ####

# linking to weather data and species names 
mammal_weather <- mammal %>% 
  left_join(., y = mam_chelsa_annual, by = c("ID", "year")) %>% 
  left_join(., y = dplyr::select(lpd_gbif, Binomial, gbif_species = gbif.species.tree),
            by = "Binomial")

glimpse(mammal_weather)

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

# pop growth and anomaly - biome
biome_temp_pgr <- ggplot(mammal_weather, aes(x = mean_temp_anomaly, y = pop_growth_rate)) +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = 0) +
  geom_point(alpha = 0.2, colour = "firebrick", size = 2.5) + 
  scale_x_continuous(breaks = seq(-0.7, 0.7, by = 0.2)) +
  facet_wrap(~ biome) +
  labs(x = "Mean annual temperature anomaly", y = "Population growth rate") +
  theme_bw(base_size = 24) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())

biome_precip_pgr <- ggplot(mammal_weather, aes(x = mean_precip_anomaly, y = pop_growth_rate)) +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = 0) +
  geom_point(alpha = 0.2, colour = "darkblue", size = 2.5) + 
  facet_wrap(~ biome) +
  labs(x = "Mean annual precipitation anomaly", y = "Population growth rate") +
  theme_bw(base_size = 24) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())

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

## 5a. Coefficients overall
coef_ovrdat <- pgr_weather %>% 
  dplyr::select(ID_block, starts_with("coef_")) %>% 
  pivot_longer(-ID_block) %>% 
  mutate(model = if_else(name %in% c("coef_temp", "coef_abun", "coef_trend") == T,
         "Temperature", "Precipitation"),
         label = case_when(
           name == "coef_temp" ~ "Mean temperature anomaly",
           name == "coef_precip" ~ "Mean precipitation anomaly",
           name == "coef_abun" ~ "Abundance - temperature",
           name == "coef_trend" ~ "Trend - temperature",
           name == "coef_abun2" ~ "Abundance - precipitation",
           name == "coef_trend2" ~ "Trend - precipitation",
         )) %>% 
  mutate(label = factor(label, 
                        levels = c("Mean temperature anomaly", 
                                   "Abundance - temperature", 
                                   "Trend - temperature",
                                   "Mean precipitation anomaly",
                                   "Abundance - precipitation",
                                   "Trend - precipitation")))

ggplot(coef_ovrdat, aes(x = value, y = label, fill = model)) +
  geom_vline(xintercept = 0) +
  geom_density_ridges(alpha = 0.7, scale = 1, size = 0.3) +
  coord_cartesian(xlim = c(-1,1)) +
  scale_fill_manual(guide = F, values = c("darkblue", "firebrick")) +
  labs(x = "Model coefficient", y = NULL) +
  theme_ridges(center_axis_labels = TRUE, font_size = 16) +
  ggsave(filename = "plots/weather_pop_growth/overall_coefficients_mnanom_5km.jpeg",
         height = 5, width = 6, units = "in", dpi = 400)

## 5b. Weather coefficients by Order
temp_sp <- ggplot(pgr_weather, aes(x = coef_temp, y = Order, 
                        fill = Order, height = stat(density))) +
  geom_vline(xintercept = 0) +
  geom_density_ridges(alpha = 0.7, scale = 2, stat = "density", size = 0.3) +
  coord_cartesian(xlim = c(-0.5,0.5)) +
  scale_fill_viridis_d(guide = F, option = "C") +
  labs(x = "Temperature anomaly coefficient", y = NULL) +
  theme_ridges(center_axis_labels = TRUE, font_size = 25) 

precip_sp <- ggplot(pgr_weather, aes(x = coef_precip, y = Order, 
                        fill = Order, height = stat(density))) +
  geom_vline(xintercept = 0) +
  geom_density_ridges(alpha = 0.7, scale = 1.5, stat = "density", size = 0.3) +
  coord_cartesian(xlim = c(-0.5,0.5)) +
  scale_fill_viridis_d(guide = F, option = "D") +
  labs(x = "Precipitation anomaly coefficient", y = NULL) +
  theme_ridges(center_axis_labels = TRUE, font_size = 25) +
  theme(axis.text.y = element_blank())

ggsave(grid.arrange(temp_sp, precip_sp, ncol = 2, widths = c(6,4)),
       filename = "plots/weather_pop_growth/coef_order_mnanom_5km.jpeg",
       width = 15, height = 13, units = "in", dpi = 400)

## 5c. Weather coefficients by biome
temp_biome <- ggplot(pgr_weather, aes(x = coef_temp, y = biome, fill = biome,
                                   height = stat(density))) +
  geom_vline(xintercept = 0) +
  geom_density_ridges(alpha = 0.7, scale = 2, stat = "density", size = 0.3) +
  coord_cartesian(xlim = c(-0.5,0.5)) +
  scale_fill_viridis_d(guide = F, option = "C") +
  labs(x = "Temperature anomaly coefficient", y = NULL) +
  theme_ridges(center_axis_labels = TRUE, font_size = 20) 

precip_biome <- ggplot(pgr_weather, aes(x = coef_precip, y = biome, fill = biome,
                                   height = stat(density))) +
  geom_vline(xintercept = 0) +
  geom_density_ridges(alpha = 0.7, scale = 1.5, stat = "density", size = 0.3) +
  coord_cartesian(xlim = c(-0.5,0.5)) +
  scale_fill_viridis_d(guide = F, option = "D") +
  labs(x = "Precipitation anomaly coefficient", y = NULL) +
  theme_ridges(center_axis_labels = TRUE, font_size = 20) +
  theme(axis.text.y = element_blank())

ggsave(grid.arrange(temp_biome, precip_biome, ncol = 2, widths = c(10,4)),
       filename = "plots/weather_pop_growth/coef_biome_mnanom_5km.jpeg",
       width = 45, height = 20, units = "cm", dpi = 400)

## 5d. Weather coefficients by latitude
pgr_lat <- pgr_weather %>% 
  mutate(lat = abs(Latitude) - (abs(Latitude) %% 22.5))

temp_lat <- ggplot(pgr_lat, aes(x = coef_temp, y = factor(lat), 
                                    fill = factor(lat), height = stat(density))) +
  geom_vline(xintercept = 0) +
  stat_density_ridges(quantile_lines = T, quantiles = 2, alpha = 0.7, scale = 1.1) +
  coord_cartesian(xlim = c(-0.5,0.5)) +
  scale_fill_viridis_d(guide = F, option = "C", begin = 0.1, end = 0.9) + 
  labs(x = "Temperature anomaly coefficient", y = "Absolute latitude") +
  theme_ridges(center_axis_labels = TRUE, font_size = 25) 

precip_lat <- ggplot(pgr_lat, aes(x = coef_precip, y = factor(lat), 
                                fill = factor(lat), height = stat(density))) +
  geom_vline(xintercept = 0) +
  stat_density_ridges(quantile_lines = T, quantiles = 2, alpha = 0.7, scale = 1.1) +
  coord_cartesian(xlim = c(-0.5,0.5)) +
  scale_fill_viridis_d(guide = F, option = "D", begin = 0.2, end = 0.8) + 
  labs(x = "Precipitation anomaly coefficient", y = NULL) +
  theme_ridges(center_axis_labels = TRUE, font_size = 25) 

ggsave(grid.arrange(temp_lat, precip_lat, ncol = 2, widths = c(7,6)),
       filename = "plots/weather_pop_growth/coef_lat_mnanom_5km.jpeg",
       width = 15, height = 11, units = "in", dpi = 400)

##__________________________________________________________________________________________________
#### 6. Save data ####

mnanom_5km <- dplyr::select(pgr_weather,
                             -c(year, raw_abundance, 
                                ln_abundance, Subspecies)) %>% 
  ungroup()

save(mnanom_5km, file = "data/pgr_weather/mnanom_5km.RData")



