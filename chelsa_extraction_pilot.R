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

##__________________________________________________________________________________________________
#### 1. Set up data ####

#_________________________________________________
## 1a. Mammal living planet data
load("mam.RData", verbose = T) # source different in the UCloud
load("../rawdata/mam.RData", verbose = T)
rm(mam_meta) # don't need the meta-data here

# 5 or more years of data and CHELSA limits
mam <- mam %>% 
  group_by(ID) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  filter(n >= 5 &
         year >= 1979 & year <= 2013) %>%
  dplyr::select(ID, Binomial, Latitude, Longitude, year)

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
  full_join(x = ., y = files_df, by = "year") %>% 
  mutate(obj_temp = paste0("temp_",year,"_",month),
         obj_precip = paste0("precip_",year,"_",month))

##__________________________________________________________________________________________________
#### 2. Extract weather data at locations in the species data ####

## Going through the files_df to save computation
x <- Sys.time()
mam_chelsa <- bind_rows(lapply(X = 1:nrow(files_df), FUN = function(x){

  # 1. extract the right data
  c_species = dplyr::filter(mam_coord,
                            year == files_df[x,]$year,
                            month == files_df[x,]$month)

  chelsa_temp   = raster(x = files_df[x,]$files_temp)
  chelsa_precip = raster(x = files_df[x,]$files_precip)

  # 2. Converting the species data to spatial data and aligning with CHELSA coordinate reference
  if(nrow(c_species) > 0){
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
  c_species$precip[c_species$precip > 65000] = NA_real_}
  
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
  c_species$precip_50m[c_species$precip_50m > 65000] = NA_real_
  
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
  c_species$precip_5km[c_species$precip_50m > 65000] = NA_real_
  
  # 4. Return
  cat("\r", "Your job is ",round(x/nrow(files_df) *100, 2), "% complete          ")
  return(as_tibble(c_species))

}))
Sys.time() - x

##__________________________________________________________________________________________________
#### 3. Save ####
save(mam_chelsa, file = "lpi_weather_pilot/mam_chelsa.RData")

