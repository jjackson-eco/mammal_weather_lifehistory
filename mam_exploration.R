#####################################################
##                                                 ##
##     Global climate and population dynamics      ##
##                                                 ##
##      Mammal Living planet data exploration      ##
##                                                 ##
##               April 7th 2020                    ##
##                                                 ##
#####################################################
rm(list = ls())
options(width = 100)

library(tidyverse)

# extra plotting 
library(grid)
library(gridExtra)
library(gganimate)
library(gifski)

# For the maps
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)
library(raster)
library(sf)
library(rasterVis)

## Load in the raw mammal data from the LPD, which we call mam
load("../rawdata/mam.RData", verbose = T)

# ID record summary
mam_IDsum <- mam %>% 
  group_by(ID) %>% 
  summarise(Binomial = Binomial[1], Order = Order[1],
            System = System[1], biome = biome[1], 
            Latitude = Latitude[1], Longitude = Longitude[1], 
            Specific_location = Specific_location[1],
            abundance_measure = abundance_measure[1],
            study_length = n(), middle_year = median(year),
            delta_sa = last(scaled_abundance) - first(scaled_abundance))

mam_spsum <- mam %>% 
  group_by(Binomial) %>% 
  summarise(study_length = n(), number_records = n_distinct(ID))

##__________________________________________________________________________________________________
#### 1. Mammal data locations ####

world_sf <- ne_coastline(scale = "medium", returnclass = "sf")

mam_map <- ggplot(data = world_sf) +
  geom_sf(size = 0.4) +
  geom_point(data = mam_IDsum, 
             aes(x = Longitude, y = Latitude, 
                 colour = Order),
             alpha = 0.5, size = 3) +
  guides(colour = guide_legend(title = NULL,
                               override.aes = list(size = 4, alpha = 1))) +
  theme_bw(base_size = 25) +
  theme(panel.grid.major = element_line(size = 0.25),
        legend.position = "bottom") +
  labs(x = "Longitude", y = "Latitude")

##__________________________________________________________________________________________________
#### 2. Database records features ####

# General database features
mam_datasum <- data.frame(Observations = nrow(mam),
                          Records = n_distinct(mam$ID),
                          Species = n_distinct(mam$Binomial),
                          Countries = n_distinct(mam_meta$Country))

General_sum <- tableGrob(mam_datasum,rows = NULL, theme = ttheme_minimal(base_size = 16))

# Order
spp_sum <- mam %>% 
  group_by(Order) %>% 
  summarise(No.species = n_distinct(Binomial))

# cl_spp <- ggplot(spp_sum, aes(x = Order, y = No.species, fill = Order)) +
#   geom_col(show.legend = FALSE) +
#   labs(x= NULL, y = "Number of Species") +
#   coord_flip() +
#   theme_bw(base_size = 16)

cl_obs <- ggplot(mam, aes(x = Order, fill = Order)) + 
  geom_bar(show.legend = F,) +
  labs(x= NULL, y = "Number of Observations") +
  coord_flip() +
  theme_bw(base_size = 16)

# Specific location
sl_obs <- ggplot(mam, aes(x = factor(Specific_location))) + 
  geom_bar(fill = "lightblue") +
  labs(x = "Specific location", y = "Number of Observations") +
  theme_bw(base_size = 16)

# Study length
studl <- ggplot(mam_IDsum, aes(x = study_length)) + 
  geom_histogram(bins = 30,fill = "lightblue") +
  geom_vline(linetype = "dashed", xintercept = median(mam_IDsum$study_length)) +
  labs(x = "Length of Study (years)", y = "Number of Records") +
  theme_bw(base_size = 16)

# Study length species
studl_spp <- ggplot(mam_spsum, aes(x = study_length)) + 
  geom_histogram(bins = 30,fill = "lightblue") +
  geom_vline(linetype = "dashed", xintercept = median(mam_spsum$study_length)) +
  labs(x = "Total years of Study", y = "Number of Species") +
  theme_bw(base_size = 16)

# Number of records species
rec_sp <- ggplot(mam_spsum, aes(x = number_records)) + 
  geom_histogram(bins = 30,fill = "lightblue") +
  geom_vline(linetype = "dashed", xintercept = median(mam_spsum$number_records)) +
  labs(x = "Number of LPD records", y = "Number of Species") +
  theme_bw(base_size = 16)

# Observations per year
obs_year <- ggplot(mam, aes(x = year)) + 
  geom_bar(show.legend = FALSE,fill = "lightblue")+
  labs(x = "Year", y = "Number of Observations") +
  theme_bw(base_size = 16)

# abundance measure
ab_measure <- ggplot(mam, aes(x = abundance_measure)) +
  geom_bar(fill = "lightblue") +
  labs(x = "Abundance measure", 
       y = "Number of observations") +
  theme_bw(base_size = 16)

# Changes in Scaled abundance?
mam_index <- dplyr::filter(mam, year >= 1970, year <= 2014)

m_ind <- ggplot(mam_index, aes(x = year, y = scaled_abundance, 
                               group = Order, colour = Order, 
                               fill = Order)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_smooth(method = "lm", alpha = 0.05) +
  labs(x = "Year", y = "Scaled abundance") +
  theme_bw(base_size = 16) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())

## PLotting 

# 1. The numbers
ggsave(grid.arrange(General_sum), 
       filename = "plots/mam_raw/mam_lpd_numbers.jpeg",
       width = 9, height = 2, units = "in",
       dpi = 400)

# 2. The distribution through time and length of study
ggsave(grid.arrange(studl, obs_year,
                    ncol = 1),
       filename = "plots/mam_raw/mam_years.jpeg", 
       width = 9, height = 8, units = "in", dpi = 400)

# 3. The map
ggsave(mam_map, 
       filename = "plots/mam_raw/mam_locations.jpeg", 
       width = 22, height = 14, units = "in", dpi = 400)

# 4. Species summaries
lay <- rbind(c(1,3),
             c(2,3))

ggsave(grid.arrange(studl_spp, rec_sp, cl_obs,
                    layout_matrix = lay),
       filename = "plots/mam_raw/mam_sp.jpeg", 
       width = 11, height = 10, units = "in", dpi = 400)

# 5. scaled abundance
ggsave(grid.arrange(ab_measure,m_ind, ncol = 1),
       filename = "plots/mam_raw/mam_scaled_abundance.jpeg",
       width = 11, height = 12, units = "in", dpi = 400)
