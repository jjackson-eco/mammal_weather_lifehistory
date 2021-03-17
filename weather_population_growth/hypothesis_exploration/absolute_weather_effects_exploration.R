#########################################################
##                                                     ##
##     Global climate and population dynamics          ##
##                                                     ##
##  Exploring absolute weather effects + life-history  ##   
##                                                     ##
##                 Nov 9th 2020                        ##
##                                                     ##
#########################################################

## Investigating absolute weather effects rather than raw weather effects with 
## respect to life-history traits - we do not expect

rm(list = ls())
options(width = 100)

library(tidyverse)
library(patchwork)
library(viridis)

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

mam_temp <- mnanom_5km_GAM %>% 
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
  drop_na(litter, longevity, bodymass)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3. Exploration plots - Longevity ####

temp_colour <- viridis(20, option = "C")[13]
precip_colour <- viridis(20, option = "D")[10]

temp_plot <- mam_temp %>% 
  mutate(abs_temp = abs(coef_temp), abs_precip = abs(coef_precip),
         long_round = round(longevity, 1)) %>% 
  group_by(long_round) %>% 
  summarise(mn_temp = mean(abs_temp), se = sd(abs_temp)/sqrt(n()), n_obs = n()) %>% 
  ggplot(aes(x = long_round, y = mn_temp, size = n_obs)) +
  geom_errorbar(aes(ymax = mn_temp + se, 
                    ymin = mn_temp - se, 
                    size = NULL), 
                width = 0.01, show.legend = F,
                colour = temp_colour) +
  geom_point(colour = temp_colour) +
  scale_size_continuous(range = c(1,8), breaks = seq(1,80, by = 20),
                        name = "Number of\nrecords") +
  labs(x = "Standardised longevity", y = "Mean |temperature effect|") +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())


precip_plot <- mam_temp %>% 
  mutate(abs_temp = abs(coef_temp), abs_precip = abs(coef_precip),
         long_round = round(longevity, 1)) %>% 
  group_by(long_round) %>% 
  summarise(mn_precip = mean(abs_precip), se = sd(abs_precip)/sqrt(n()), n_obs = n()) %>% 
  ggplot(aes(x = long_round, y = mn_precip, size = n_obs)) +
  geom_errorbar(aes(ymax = mn_precip + se, 
                    ymin = mn_precip - se, size = NULL), 
                width = 0.01, show.legend = F,
                colour = precip_colour) +
  geom_point(colour = precip_colour) +
  scale_size_continuous(range = c(1,8), breaks = seq(1,80, by = 20),
                        name = "Number of\nrecords") +
  labs(x = "Standardised longevity", y = "Mean |precipitation effect|") +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())


ggsave(temp_plot + precip_plot,
       filename = "plots/phylogenetic_regression/absolute_weather_effect_lh.jpeg",
       width = 26, height = 10, units = "cm", dpi = 500)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 4. Exploration plots - all demographic rates ####

mam_temp %>% 
  dplyr::rename(Adult_bodymass = bodymass, Litter_size = litter,
                Maximum_longevity = longevity) %>% 
  mutate(Temperature = abs(coef_temp), Precipitation = abs(coef_precip)) %>% 
  pivot_longer(c(Temperature, Precipitation), 
               names_to = "weather", values_to = "coef") %>% 
  pivot_longer(c(Adult_bodymass, Litter_size,
                 Maximum_longevity),
               names_to = "lifehistory", 
               values_to = "lifehistory_value") %>% 
  mutate(lifehistory = gsub(pattern = "_",
                            replacement = " ", 
                            lifehistory),
         lifehistory_value = round(lifehistory_value, 1)) %>%
  group_by(weather, lifehistory, lifehistory_value) %>%
  summarise(mn_coef = mean(coef, na.rm = TRUE), 
            se = sd(coef, na.rm = TRUE)/sqrt(n()),
            n_records = n()) %>% 
  ggplot(aes(x = lifehistory_value, y = mn_coef, 
             colour = weather, size = n_records)) +
  geom_errorbar(aes(ymax = mn_coef + se, 
                    ymin = mn_coef - se,
                    size = NULL),
                width = 0.05, show.legend = F) +
  geom_point(show.legend = F) +
  facet_wrap(~ weather + lifehistory) +
  scale_colour_manual(values = c(precip_colour, temp_colour)) +
  labs(x = "Standardised life history value", y = "Mean |weather effect|") +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank()) +
  ggsave(filename = "plots/weather_pop_growth/absolute_weather_life_history_ALL.jpeg",
         width = 20, height = 18, units = "cm", dpi = 500)



