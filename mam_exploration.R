#####################################################
##                                                 ##
##     Global climate and population dynamics      ##
##                                                 ##
##          Mammal LPI data exploration            ##
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

## Load in the raw mammal data from the LPI, which we call mam
load("../rawdata/mam.RData", verbose = T)

# ID record summary
mam_IDsum <- mam %>% 
  group_by(ID) %>% 
  summarise(Binomial = Binomial[1], Order = Order[1],
            System = System[1], FW_biome = FW_biome[1], 
            Latitude = Latitude[1], Longitude = Longitude[1], 
            Specific_location = Specific_location[1],
            study_length = n(), middle_year = median(year),
            delta_sa = last(scaled_abundance) - first(scaled_abundance))

mam_spsum <- mam %>% 
  group_by(Binomial) %>% 
  summarise(study_length = n(), number_records = n_distinct(ID))

##__________________________________________________________________________________________________
#### 1. Mammal LPI data locations ####

world_sf <- ne_coastline(scale = "medium", returnclass = "sf")

ggplot(data = world_sf) +
  geom_sf(size = 0.1) +
  geom_point(data = mam_IDsum, 
             aes(x = Longitude, y = Latitude, 
                 colour = Order),
             alpha = 0.6, size = 2) +
  theme_bw(base_size = 25) +
  theme(panel.grid.major = element_line(size = 0.25)) +
  labs(x = "Longitude", y = "Latitude") +
  ggsave(filename = "plots/mam_raw/mam.jpeg", 
         width = 20, height = 12, units = "in", dpi = 400)

##__________________________________________________________________________________________________
#### 2. Database records features ####

# General database features
mam_datasum <- data.frame(Observations = nrow(mam),
                          Records = n_distinct(mam$ID),
                          Species = n_distinct(mam$Binomial),
                          Countries = n_distinct(mam_meta$Country))

General_sum <- tableGrob(mam_datasum,rows = NULL, theme = ttheme_minimal(base_size = 15))

# Order
spp_sum <- mam %>% 
  group_by(Order) %>% 
  summarise(No.species = n_distinct(Binomial))

cl_spp <- ggplot(spp_sum, aes(x = Order, y = No.species, fill = Order)) +
  geom_col(show.legend = FALSE) +
  labs(x= NULL, y = "Number of Species") +
  coord_flip() +
  theme_bw(base_size = 16)

cl_obs <- ggplot(mam, aes(x = Order, fill = Order)) + 
  geom_bar(show.legend = F,) +
  labs(x= NULL, y = "Number of Observations") +
  scale_x_discrete(labels = NULL) +
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
  labs(x = "Total years of Study", y = "Number of Species") +
  theme_bw(base_size = 16)

# Observations per year
obs_year <- ggplot(mam, aes(x = year)) + 
  geom_bar(show.legend = FALSE,fill = "lightblue")+
  labs(x = "Year", y = "Number of Observations") +
  theme_bw(base_size = 16)

# Changes in Scaled abundance?
mam_index <- dplyr::filter(mam, year >= 1970, year <= 2014)

m_ind <- ggplot(mam_index, aes(x = year, y = scaled_abundance, group = Order)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.2, colour = "lightblue") +
  geom_smooth(method = "lm", se = F) +
  labs(x = "Year", y = "Scaled abundance") +
  facet_wrap(~Order, ncol = 4) +
  theme_bw(base_size = 16) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())

## PLotting 
ggsave(grid.arrange(General_sum, sl_obs, studl, 
                    studl_spp, rec_sp, obs_year,
                    ncol = 3),
       filename = "plots/mam_raw/mam_summary.jpeg", 
       width = 15, height = 14, units = "in", dpi = 400)

ggsave(grid.arrange(cl_spp, cl_obs, ncol = 2, widths = c(6,5)),
       filename = "plots/mam_raw/mam_order.jpeg", 
       width = 9, height = 10, units = "in", dpi = 400)

ggsave(m_ind, filename = "plots/mam_raw/mam_scaled_abundance.jpeg",
       width = 9, height = 10, units = "in", dpi = 400)
