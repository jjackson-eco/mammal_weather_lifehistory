#####################################################
##                                                 ##
##     Global climate and population dynamics      ##
##                                                 ##
##     CHELSA weather data extraction - PILOT      ##
##                                                 ##
##               March 25th 2020                   ##
##                                                 ##
#####################################################
rm(list = ls())
options(width = 100)

# general
library(tidyverse)
library(grid)
library(gridExtra)

# spatial
library(sp)
library(rgeos)
library(raster)
library(exactextractr)
library(rgdal)
library(sf)

library(rnaturalearth)
library(rnaturalearthdata)

##__________________________________________________________________________________________________
#### 1. Set up data ####

#_________________________________________________
## 1a. Mammal living planet data
load("lpi_weather_pilot/mam.RData", verbose = T) # source different in the UCloud
load("../rawdata/mam.RData", verbose = T)
rm(mam_meta) # don't need the meta-data here

# Keeping only years from CHELSA and only 5 years of study beyond that.
# Expand for all month-year combinations to normalise the climate data after.
mam <- mam %>% 
  filter(year >= 1979 & year <= 2013) %>% # important to do this first
  group_by(ID) %>% 
  summarise(Binomial = Binomial[1], 
            Latitude = Latitude[1],
            Longitude = Longitude[1], 
            n = n()) %>% 
  ungroup() %>% 
  filter(n >= 5) %>% 
  dplyr::select(-n) %>% 
  expand_grid(., year = 1979:2013, month = 1:12) 

#__________________________________________________
## 1b. CHELSA raster file lists

# from the network drive
# netdr_t <- "/Volumes/PlantAnimalBiodemography/Climate/chelsa/temp/"
# netdr_p <- "/Volumes/PlantAnimalBiodemography/Climate/chelsa/prec/"

# from the UCloud - you select the files you want to use when creating the job
U_t <- "CHELSA/temperature"
U_p <- "CHELSA/precipitation"

files_temp <- paste(U_t, list.files(U_t), sep = "/")
files_precip <- paste(U_p, list.files(U_p), sep = "/")

files_df <- tibble(files_temp, files_precip) %>% 
  mutate(year = map_int(files_temp, ~ as.integer(strsplit(.x, "_")[[1]][3]))) %>% 
  mutate(month = map_int(files_temp, ~ as.integer(strsplit(.x, "_")[[1]][4])))

## 1c. Merge species-coordinate data with files_df data
mam_coord <- mam %>% 
  left_join(x = ., y = files_df, 
            by = c("year","month"))

##__________________________________________________________________________________________________
#### 2. Extract weather data at locations in the species data ####

## Going through the files_df to save computation
starttime <- Sys.time()
mam_chelsa <- bind_rows(lapply(X = 1:nrow(files_df), FUN = function(x){
  
  # 1. extract the right data
  c_species = dplyr::filter(mam_coord,
                            year == files_df[x,]$year,
                            month == files_df[x,]$month)
  
  chelsa_temp   = raster(x = files_df[x,]$files_temp)
  chelsa_precip = raster(x = files_df[x,]$files_precip)
  
  # incorporating NAs
  chelsa_precip[values(chelsa_precip) > 65000] = NA_real_
  
  # 2. Converting the species data to spatial data and aligning with CHELSA coordinate reference
  coordinates(c_species) = c("Longitude","Latitude")
  pbase = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  proj4string(c_species) = pbase
  c_species = spTransform(c_species, crs(chelsa_temp))
  
  #_________________________________________________________________________________________________
  # 3a. Extract the climate data from exact points- Takes a while
  cells_sp = cellFromXY(chelsa_temp, c_species)
  c_species$raster_cell = cells_sp
  c_species$temp = chelsa_temp[cells_sp] / 10 - 273.15 # Kelvin to celcius
  c_species$precip = chelsa_precip[cells_sp]
  
  #_________________________________________________________________________________________________
  # 3b. Extracting the climate data from a small buffer radius of 50m - if ~identical will use this
  ## Key advantage is using exactextractr, which uses C++ coding to do raster operations much faster.
  ## First convert to the Universal Transverse Mercator projection system to buffer in m accurately:
  ## Useful information on projection systems at https://rspatial.org/raster/spatial/6-crs.html
  p_utm = CRS("+proj=utm +datum=WGS84")
  utm = spTransform(c_species, p_utm)
  # do it by the study ID for future reference
  # 50 line segments for the approximation of a circle
  csp_buff_small = rgeos::gBuffer(utm, width = 50, quadsegs = 50,
                                  byid = TRUE, id = utm$ID)
  ## Convert back to WGS84 to extract from CHELSA and convert to sf type
  csp_buff_wgs_small = sf::st_as_sf(spTransform(csp_buff_small, 
                                                CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")))
  c_species$temp_50m = (exact_extract(x = chelsa_temp, y = csp_buff_wgs_small, 
                                      fun = "mean", progress = FALSE)/10) - 273.15
  c_species$precip_50m = exact_extract(x = chelsa_precip, y = csp_buff_wgs_small, 
                                       fun = "mean", progress = FALSE)
  
  #_________________________________________________________________________________________________
  # 3c. Large buffer radius of 5km - More biologically relevant?
  csp_buff = rgeos::gBuffer(utm, width = 5000, quadsegs = 50,
                            byid = TRUE, id = utm$ID)
  csp_buff_wgs = sf::st_as_sf(spTransform(csp_buff, 
                                          CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")))
  c_species$temp_5km = (exactextractr::exact_extract(x = chelsa_temp, y = csp_buff_wgs, 
                                                     fun = "mean", progress = FALSE)/10) - 273.15
  c_species$precip_5km = exactextractr::exact_extract(x = chelsa_precip, y = csp_buff_wgs, 
                                                      fun = "mean", progress = FALSE)
  
  #_________________________________________________________________________________________________
  # 3d. V Large buffer radius of 50km - More appropriate for species with large ranges?
  csp_buff_large = rgeos::gBuffer(utm, width = 50000, quadsegs = 50,
                                  byid = TRUE, id = utm$ID)
  csp_buff_wgs_large = sf::st_as_sf(spTransform(csp_buff_large, 
                                                CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")))
  c_species$temp_50km = (exactextractr::exact_extract(x = chelsa_temp, y = csp_buff_wgs_large, 
                                                      fun = "mean", progress = FALSE)/10) - 273.15
  c_species$precip_50km = exactextractr::exact_extract(x = chelsa_precip, y = csp_buff_wgs_large, 
                                                       fun = "mean", progress = FALSE)
  
  # 4. Return
  cat("\r", "Your job is ",round(x/nrow(files_df) *100, 2), "% complete          ")
  return(as_tibble(c_species))
  
}))
Sys.time() - starttime

##__________________________________________________________________________________________________
#### 3. Save ####

# This operation is now done in sections on the UCloud, so files are saved elsewhere.

# saveRDS(mam_chelsa, file = "lpi_weather_pilot/mam_chelsa.RDS") # On the UCloud

# Load data downloaded from the UCloud - in chunks to save computation time on the UCloud
mam_chelsa1979 <- readRDS("data/mam_chelsa_cloud/mam_chelsa1979.RDS")
mam_chelsa1986 <- readRDS("data/mam_chelsa_cloud/mam_chelsa1986.RDS") 
mam_chelsa1993 <- readRDS("data/mam_chelsa_cloud/mam_chelsa1993.RDS") 
mam_chelsa2000 <- readRDS("data/mam_chelsa_cloud/mam_chelsa2000.RDS") 
mam_chelsa2007 <- readRDS("data/mam_chelsa_cloud/mam_chelsa2007.RDS")

mam_chelsa <- bind_rows(mam_chelsa1979,
                        mam_chelsa1986,
                        mam_chelsa1993,
                        mam_chelsa2000,
                        mam_chelsa2007)

saveRDS(mam_chelsa, file = "data/mam_chelsa.RDS")

rm(mam_chelsa1979,
   mam_chelsa1986,
   mam_chelsa1993,
   mam_chelsa2000,
   mam_chelsa2007)

##__________________________________________________________________________________________________
#### 4. Plots to compare the CHELSA weather variables ####

# Data
temp_chelsa <- mam_chelsa %>% 
  dplyr::select(ID, Binomial, year, month, 
                temp_exact_cell = temp, 
                temp_buffer_50m = temp_50m, 
                temp_buffer_5km = temp_5km, 
                temp_buffer_50km = temp_50km) %>% 
  pivot_longer(data = ., 
               cols = starts_with("temp_buffer"),
               names_to = "scale", names_prefix = "temp_",
               values_to = "temperature") %>% 
  mutate(scale = factor(scale, levels = c("buffer_50m", "buffer_5km", "buffer_50km")))

precip_chelsa <- mam_chelsa %>% 
  dplyr::select(ID, Binomial, year, month, 
                precip_exact_cell = precip, 
                precip_buffer_50m = precip_50m, 
                precip_buffer_5km = precip_5km, 
                precip_buffer_50km = precip_50km) %>% 
  pivot_longer(data = ., 
               cols = starts_with("precip_buffer"),
               names_to = "scale", names_prefix = "precip_",
               values_to = "precipitation") %>% 
  mutate(scale = factor(scale, levels = c("buffer_50m", "buffer_5km", "buffer_50km")))

# Plots
temp_buffer_plot <- ggplot(temp_chelsa, aes(x =  temperature, y = temp_exact_cell)) +
  geom_point(alpha = 0.01, colour = "firebrick") +
  geom_abline(slope = 1, intercept = 0, size = 0.2, colour = "black") +
  labs(x = expression(paste("Average temperature", ~degree~C)),
       y = expression(paste("Exact raster cell temperature", ~degree~C))) +
  facet_wrap(~ scale, ncol = 3) +
  theme_bw(base_size = 17) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())

precip_buffer_plot <- ggplot(precip_chelsa, aes(x =  precipitation, y = precip_exact_cell)) +
  geom_point(alpha = 0.01, colour = "lightblue") +
  geom_abline(slope = 1, intercept = 0, size = 0.2, colour = "black") +
  labs(x = "Average precipitation (mm)",
       y = "Exact raster cell precipitation (mm)") +
  facet_wrap(~ scale, ncol = 3) +
  theme_bw(base_size = 17) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())

ggsave(grid.arrange(temp_buffer_plot,precip_buffer_plot, nrow = 2),
      filename = "plots/mam_chelsa/weather_buffer.jpeg",
      units = "in", width = 9, height = 8, dpi = 300)


