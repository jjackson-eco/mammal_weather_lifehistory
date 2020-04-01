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

library(tidyverse)
library(grid)
library(gridExtra)
library(raster)

##__________________________________________________________________________________________________
#### 1. Set up data ####

## 1a. Mammal living planet index data
load("../mam.RData", verbose = T)
rm(mam_meta) # don't need the meta-data here

# 5 or more years of data and CHELSA limits
mam <- mam %>% 
  group_by(ID) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  filter(n >= 5 &
         year >= 1979 & year <= 2013) %>%
  dplyr::select(ID, Binomial, Latitude, Longitude, year)

## 1b. CHELSA raster file lists (via network drive) as a tibble
netdr_t <- "/Volumes/PlantAnimalBiodemography/Climate/chelsa/temp/"
netdr_p <- "/Volumes/PlantAnimalBiodemography/Climate/chelsa/prec/"

files_temp <- paste(netdr_t, list.files(netdr_t), sep = "/")
files_precip <- paste(netdr_p, list.files(netdr_p), sep = "/")

files_df <- tibble(files_temp, files_precip) %>% 
  mutate(year = map_int(files_temp, ~ as.integer(strsplit(.x, "_")[[1]][3]))) %>% 
  mutate(month = map_int(files_temp, ~ as.integer(strsplit(.x, "_")[[1]][4])))

## 1c. Merge species-coord data with files_df data
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
  
  cells_sp = cellFromXY(chelsa_temp, c_species)

  # 3. Extract the climate data - Takes a while
  c_species$raster_cell = cells_sp
  c_species$temp = chelsa_temp[cells_sp] / 10 - 273.15 # Kelvin to celcius
  c_species$precip = chelsa_precip[cells_sp]
  c_species$precip[c_species$precip > 65000] = NA_real_}
  
  # 4. Return
  cat("\r", "Your job is ",round(x/nrow(files_df) *100, 2), "% complete          ")
  return(as_tibble(c_species))

}))
Sys.time() - x

##__________________________________________________________________________________________________
#### 3. Save ####
save(mam_chelsa, file = "data/mam_chelsa.RData")
