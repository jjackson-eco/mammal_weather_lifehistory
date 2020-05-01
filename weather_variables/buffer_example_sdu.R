#####################################################
##                                                 ##
##     Global climate and population dynamics      ##
##                                                 ##
##               Buffering example                 ##
##                                                 ##
##                 May 1st 2020                    ##
##                                                 ##
#####################################################
rm(list = ls())
options(width = 100)

library(tidyverse)
library(sp)
library(raster)
library(rgeos)
library(rgdal)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

###__________________________________________________________________________________________________
#### 1. Data set up ####

# countries
world_sf <- ne_countries(scale = "medium", returnclass = "sf") # defaults to WGS84 projection so all good

# study location
sdu <- tibble(stud = "SDU", Long = 10.425869, Lat = 55.368266)

# convert to a spatial dataframe and set the inital projection
coordinates(sdu) <- c("Long","Lat")
p_wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
proj4string(sdu) <- p_wgs
crs(sdu)

###__________________________________________________________________________________________________
#### 2. Creating buffers ####

# converting to a UTM coordinate system to calculate buffer in meters
p_utm <- CRS("+proj=utm +datum=WGS84")
utm <- spTransform(sdu, p_utm)

# creating the buffer polygons and converting back to wgs
sdu_5km <-  rgeos::gBuffer(utm, width = 5000, quadsegs = 50,
                       byid = TRUE, id = utm$stud)
sdu_5km <- sf::st_as_sf(spTransform(sdu_5km, p_wgs))

sdu_50km <-  rgeos::gBuffer(utm, width = 50000, quadsegs = 50,
                           byid = TRUE, id = utm$stud)
sdu_50km <- sf::st_as_sf(spTransform(sdu_50km, p_wgs))

###__________________________________________________________________________________________________
#### 3. Plot ####

ggplot(data = world_sf) +
  geom_sf(size = 0.1) +
  geom_sf(data = sdu_50km, 
          aes(fill = "50km buffer"),
          alpha = 0.2, size = 0.2) +
  geom_sf(data = sdu_5km,
          aes(fill = "5km buffer"), 
          alpha = 0.2, size = 0.2) +
  geom_point(data = as_tibble(sdu), 
             aes(x = Long, y = Lat), size = 0.2) +
  coord_sf(xlim = c(8, 14), ylim = c(54, 58)) +
  labs(x = NULL, y = NULL) +
  scale_fill_manual(name = NULL, values = c("lightgreen", "blue")) +
  theme_bw(base_size = 17) +
  theme(legend.position = c(0.83,0.85)) +
  ggsave(filename = "plots/developing_annual_weather_variables/buffer_example_sdu.jpeg",
         width = 7, height = 7, units = "in", dpi = 400)

  







