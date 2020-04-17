#####################################################
##                                                 ##
##     Global climate and population dynamics      ##
##                                                 ##
##          CHELSA weather exploration             ##
##                                                 ##
##               April 1st 2020                    ##
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

temp_aug1993 <- raster(x = "../rawdata/CHELSA_tmean_1993_08_V1.2.1.tif") 
precip_aug1993 <- raster(x = "../rawdata/precip_aug1993.tif") 

# coastline
world_sp <- ne_coastline(scale = "medium")
world_sf <- ne_coastline(scale = "medium", returnclass = "sf")

##__________________________________________________________________________________________________
#### 1. Raw raster data ####

# using rasterVis 

# 1a. Temperature
gplot(temp_aug1993) + 
  geom_tile(aes(fill=value)) +
  geom_polygon(data = world_sp,
               aes(x=long, y=lat, group = group),
               size = 0.2, fill = NA, colour = "black") +
  scale_fill_viridis_c(option = "inferno",
                       breaks = (seq(-60,60, by = 10) + 273.15)*10,
                       labels = seq(-60,60, by = 10)) +
  scale_x_continuous(breaks = seq(-180,180, by = 90), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(-90,90, by = 45), expand = c(0,0)) +
  guides(fill = guide_colorbar(title = expression(paste("Temperature", ~degree~C)),
                                 barheight = 21, barwidth = 3)) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 20) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank()) +
  ggsave(filename = "plots/chelsa_raw/temp_aug93.jpeg",
         width = 20, height = 12, units = "in",dpi = 400)

# 1b. Precipitation
gplot(precip_aug1993) + 
  geom_tile(aes(fill=value)) +
  geom_polygon(data = world_sp,
               aes(x=long, y=lat, group = group),
               size = 0.2, fill = NA, colour = "black") +
  scale_fill_viridis_c(option = "viridis",na.value = "white") +
  scale_x_continuous(breaks = seq(-180,180, by = 90), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(-90,90, by = 45), expand = c(0,0)) +
  guides(fill = guide_colorbar(title = "Total\nPrecipitation (mm)",
                               barheight = 21, barwidth = 3)) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 20) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank()) +
  ggsave(filename = "plots/chelsa_raw/precip_aug93.jpeg",
         width = 20, height = 12, units = "in",dpi = 400)

##__________________________________________________________________________________________________
#### 2. Raw CHELSA data for mammal study sites ####

# adding year/month as a factor
month_dat <- data.frame(month = 1:12, 
                        month_name = factor(month.name, levels = month.name)) %>% 
  expand_grid(year = 1979:2013, .) %>% 
  mutate(year_dec = year + (month/12),
         year_mon = paste0(year, " ", month_name))

month_dat$year_mon <- factor(month_dat$year_mon, levels = month_dat$year_mon)

mam_chelsa <- mam_chelsa %>% 
  left_join(x = ., y = month_dat,
            by = c("year", "month"))

#________________________________________________________
## 1a. Temperature
ggplot(data = world_sf) +
  geom_sf(size = 0.1) +
  geom_point(data = mam_chelsa, 
             aes(x = Longitude, y = Latitude, 
                 colour = temp),
             alpha = 0.6, size = 2) +
  scale_colour_viridis_c(option = "inferno") +
  guides(colour = guide_colorbar(title = expression(paste("Temperature", ~degree~C)),
                                 barheight = 25, barwidth = 5)) +
  facet_wrap(~month_name, ncol = 3) +
  theme_bw(base_size = 25) +
  theme(panel.grid.major = element_line(size = 0.25),
        strip.background = element_blank()) +
  labs(x = "Longitude", y = "Latitude") +
  ggsave(filename = "plots/chelsa_raw/temp_mam.jpeg", 
         width = 23, height = 19, units = "in", dpi = 400)

# Additions for animation
# labs(x = "Longitude", y = "Latitude", title = "{closest_state}") +
# transition_states(year_mon) +
# ease_aes('linear') 

# anim_save(filename = "plots/chelsa_raw/temp_mam.gif", 
#           animation = animate(temp_mam,fps = 6,
#                               width = 11, height = 9, 
#                               units = "in", res = 300))

#________________________________________________________
## 1b. Precipitation
ggplot(data = world_sf) +
  geom_sf(size = 0.1) +
  geom_point(data = mam_chelsa, 
             aes(x = Longitude, y = Latitude, 
                 colour = log(precip)),
             alpha = 0.6, size = 2) +
  scale_colour_viridis_c(option = "viridis") +
  guides(colour = guide_colorbar(title = "Precipitation\n(log scale)",
                                 barheight = 25, barwidth = 5)) +
  facet_wrap(~month_name, ncol = 3) +
  theme_bw(base_size = 25) +
  theme(panel.grid.major = element_line(size = 0.25),
        strip.background = element_blank()) +
  labs(x = "Longitude", y = "Latitude") +
  ggsave(filename = "plots/chelsa_raw/precip_mam.jpeg", 
         width = 23, height = 19, units = "in", dpi = 400)




