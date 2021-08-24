#############################################################
##                                                         ##
##          Global climate and population dynamics         ##
##                                                         ##
##    Testing for spatial autocorrelation in GAM models    ##
##                                                         ##
##                    Aug 23rd 2021                        ##
##                                                         ##
#############################################################

## Using Morans I to test whether there is spatial autocorrelation in the 
## GAM population responses across the mammals

# https://rpubs.com/quarcs-lab/spatial-autocorrelation for an example

rm(list = ls())
options(width = 100)

library(tidyverse)
library(sf)
library(spdep)
library(viridis)
library(rnaturalearth)
library(patchwork)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 1. Data in a spatial format ####

load("data/mammal_analysis_data_GAM.RData", verbose = TRUE)

## Specify the coordinates
mam_sp <- mam_coef %>% 
  dplyr::select(id, coef_temp, coef_precip, Longitude, Latitude) %>% 
  drop_na(coef_precip) # Keeping only non-missing values from precipitation effects

coordinates(mam_sp) <- ~ Longitude + Latitude

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 2. Convert spatial points to a neighbours list for analyses ####

# return k nearest neighbours for each coordinate point
knea <- knearneigh(coordinates(mam_sp), longlat = TRUE)

# convert to a neighbours list
neighbours <- knn2nb(knea)

# plot 
world_sp <- ne_countries(scale = "medium", returnclass = "sp")
plot(world_sp, border = 'lightgrey')
plot(neighbours, coordinates(mam_sp), add=TRUE, col='blue') # should they not all be connected? try reading around to check

# convert neighbours list to a weights matrix for analysis
listw <- nb2listw(neighbours)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3. Morans I tests ####

morans_temp <- moran.mc(mam_sp$coef_temp, listw, nsim = 1000)
morans_precip <- moran.mc(mam_sp$coef_precip, listw, nsim = 1000)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 4. Exploring spatial autocorrelation in temperature more ####

## Plot of the spatially lagged temperature coefficients with the actual values
jpeg(filename = "plots/meta_regression/spatial_autocorrelation_temp.jpeg",
     width = 15, height = 20, units = "cm",res = 500)
par(mfrow = c(2,1))
plot(world_sp, border = 'lightgrey', main = "Nearest neighbours plot")
plot(neighbours, coordinates(mam_sp), add=TRUE, col='blue')
moran.plot(mam_sp$coef_temp, listw = nb2listw(neighbours, style = "W"),
           xlab = "Temperature effect", ylab = "Spatially lagged temperature effect") # not a convincing correlation
dev.off()

#______________________________________________________
## Local morans I for temp
listW_l <- nb2listw(neighbours, style = "W")
moranslocal_temp <- localmoran(mam_sp$coef_temp, listW_l)

## In base R - laaaaaame
# colour palette and transform the numeric variable in bins
moran_mapdat <- cbind(mam_sp, moranslocal_temp)
cols <- viridis(10, option = "B")

rank <- as.factor(as.numeric( cut(moran_mapdat$Ii, 10)))

plot(world_sp, border = 'lightgrey', lwd = 0.6)
points(moran_mapdat, bg = alpha(cols[rank], 0.5), pch = 21, col = NULL)


## ggplot using sf
moran_mapdat_sf <- mam_coef %>% 
  dplyr::select(id, coef_temp, coef_precip, Longitude, Latitude) %>% 
  drop_na(coef_precip) %>% 
  bind_cols(as.data.frame(moranslocal_temp)) %>% 
  mutate(sig = if_else(`Pr(z > 0)` < 0.05, "< 0.05", ">= 0.05")) %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), 
           crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

world_sf <- ne_countries(scale = "medium", returnclass = "sf")

ggplot(data = world_sf) +
  geom_sf(size = 0, fill = "lightgrey") + 
  geom_sf(data = moran_mapdat_sf, aes(colour = sig), alpha = 0.3, size = 2) +
  geom_sf(data = filter(moran_mapdat_sf, sig == "< 0.05"), 
          aes(colour = sig), alpha = 1, size = 2) +
  scale_colour_viridis_d(name = "Local morans I significance", begin = 0.3, end = 0.8) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position = "top") +
  ggsave("plots/meta_regression/spatial_autocorrelation_localtemp.jpeg",
         width = 20, height = 16, units = "cm",dpi = 500)

## Total number of significant local morans I values
table(moran_mapdat_sf$sig)

24/454

## Looking at the local morans I values and points from Asia in more detail

world_Ii <- ggplot(data = world_sf) +
  geom_sf(size = 0, fill = "lightgrey") + 
  geom_sf(data = moran_mapdat_sf, aes(colour = Ii), alpha = 0.2, size = 1) +
  geom_sf(data = filter(moran_mapdat_sf,sig == "< 0.05"), 
          aes(colour = Ii), alpha = 0.9, size = 1.5) +
  scale_colour_viridis_c(name = "Local morans I value", option = "C", begin = 0.2, end = 0.9) +
  guides(colour = guide_colorbar(title.position = "top")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position = "top") 

USA_Ii <- ggplot(data = world_sf) +
  geom_sf(size = 0, fill = "lightgrey") + 
  geom_sf(data = moran_mapdat_sf, aes(colour = Ii, shape = sig), alpha = 0.8, size = 1.5) +
  coord_sf(xlim = c(-120, -60), ylim = c(20, 60)) +
  scale_colour_viridis_c(name = "Local morans I value", option = "C", begin = 0.2, end = 0.9) +
  guides(colour = guide_colorbar(title.position = "top"),
         shape = guide_legend(title.position = "top",
                              title = "Significance",
                              override.aes = list(size = 5))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position = "top") 

ggsave(world_Ii + USA_Ii, 
       filename = "plots/meta_regression/spatial_autocorrelation_localvalues_temp.jpeg",
       width = 24, height = 18, units = "cm", dpi = 600)




