##############################################################
##                                                          ##
##          Global climate and population dynamics          ##
##                                                          ##
## Range of weather anomalies across life-history variables ##
##                                                          ##
##                      Mar 2nd 2021                        ##
##                                                          ##
##############################################################

## Exploration to see whether there is a good range of raw weather 
## anomalies associated with life-history traits.

rm(list = ls())
options(width = 100)

library(tidyverse)
library(viridis)
library(patchwork)
library(flextable)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 1. Load and merge data ####

# mammal raw data
load("../rawdata/mammal.RData")
glimpse(mammal)

# weather
mam_chelsa_annual <- readRDS("data/mam_chelsa_annual.RDS") %>% 
  filter(scale == "scale_5km") %>% 
  dplyr::select(ID,year, weather_scale = scale, mean_temp_anomaly, mean_precip_anomaly)
glimpse(mam_chelsa_annual)

# life-history data
load("data/lifehistory.RData", verbose = TRUE)

# species names to merge
load("../rawdata/GBIF_species_names_mamUPDATE.RData", verbose = TRUE)

# linking to weather data and species names 
mammal_weather <- mammal %>% 
  left_join(., y = mam_chelsa_annual, by = c("ID", "year")) %>% 
  left_join(., y = dplyr::select(lpd_gbif, Binomial, gbif.species.id),
            by = "Binomial") %>% 
  left_join(x = ., y = dplyr::select(lifehistory, c(3,6:9)),
            by = "gbif.species.id") %>% 
  mutate(year_s = as.numeric(scale(year))) %>% 
  drop_na(litter, longevity, bodymass)

glimpse(mammal_weather)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 2. Explore anomaly range per life-history variable ####

temp_colour <- "#990a80"
precip_colour <- "#287f79"

summary(mammal_weather$longevity)

plotdat <- mammal_weather %>% 
  dplyr::rename(Temperature = mean_temp_anomaly, 
                Precipitation = mean_precip_anomaly) %>% 
  pivot_longer(cols = c(Temperature, Precipitation), 
               names_to = "weather_var", values_to = "anomaly") %>% 
  dplyr::select(ID_block, Binomial, year, longevity, litter, bodymass,
                weather_var, anomaly) %>% 
  # life-history data bins of 0.2 units
  mutate(litter_bin = litter - (litter %% 0.2),
         longevity_bin = longevity - (longevity %% 0.2),
         bodymass_bin = bodymass - (bodymass %% 0.2)) 

# longevity
longplot <- ggplot(plotdat, aes(x = longevity_bin, 
                    group = interaction(longevity_bin, weather_var),
                    y = anomaly, fill = weather_var)) +
  geom_hline(yintercept = 0) +
  geom_violin() +
  scale_fill_manual(values = c(precip_colour, temp_colour), guide = F) +
  labs(x = "Standardised longevity (binned in 0.2 increments)",
       y = "Mean weather anomaly") +
  facet_wrap(~weather_var, scales = "free") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank())

# litter
litplot <- ggplot(plotdat, aes(x = litter_bin, 
                                group = interaction(litter_bin, weather_var),
                                y = anomaly, fill = weather_var)) +
  geom_hline(yintercept = 0) +
  geom_violin() +
  scale_fill_manual(values = c(precip_colour, temp_colour), guide = F) +
  labs(x = "Standardised litter size (binned in 0.2 increments)",
       y = "Mean weather anomaly") +
  facet_wrap(~weather_var, scales = "free") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank())

# bodymass
bodplot <- ggplot(plotdat, aes(x = bodymass_bin, 
                               group = interaction(bodymass_bin, weather_var),
                               y = anomaly, fill = weather_var)) +
  geom_hline(yintercept = 0) +
  geom_violin() +
  scale_fill_manual(values = c(precip_colour, temp_colour), guide = F) +
  labs(x = "Standardised bodymass (binned in 0.2 increments)",
       y = "Mean weather anomaly") +
  facet_wrap(~weather_var, scales = "free") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank())

ggsave(longplot / litplot / bodplot,
       filename = "plots/weather_pop_growth/raw_anomaly_life_history_check.jpeg",
       width = 15, height = 18, units = "cm",dpi =500)


