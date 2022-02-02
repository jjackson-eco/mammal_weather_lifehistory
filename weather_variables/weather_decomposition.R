#####################################################
##                                                 ##
##     Global climate and population dynamics      ##
##                                                 ##
##     CHELSA weather data - STL decomposition     ##
##                                                 ##
##               April 17th 2020                   ##
##                                                 ##
#####################################################
rm(list = ls())
options(width = 100)

library(tidyverse)
library(grid)
library(gridExtra)
library(xts)
library(MASS)

library(sf)

library(rnaturalearth)

# load the data
mam_chelsa <- readRDS("data/mam_chelsa.RDS")

##__________________________________________________________________________________________________
#### 1. Converting mam chelsa to have dates ####

## Test - temp for ID 28 - Cheetah
temp_test <- mam_chelsa %>% 
  filter(ID == 28) %>% 
  dplyr::select(ID, Binomial, 
                Longitude, Latitude, 
                year, month, temp) %>% 
  mutate(date = as.Date(paste0(year, "-", month, "-", "01")), # first day of the month - justified? matters?
         temp = as.numeric(scale(temp))) 

## 1a. Splitting temperature and precipitation for ease of merging later on

## All temperature data
temp_dat <- mam_chelsa %>% 
  dplyr::select(ID, Binomial, Longitude, 
                Latitude, year, month, 
                temp, temp_50m, 
                temp_5km, temp_50km) %>% 
  mutate(date = as.Date(paste0(year, "-", month, "-", "01"))) 

## All precipitation data
precip_dat <- mam_chelsa %>% 
  dplyr::select(ID, Binomial, Longitude, 
                Latitude, year, month, 
                precip, precip_50m, 
                precip_5km, precip_50km) %>% 
  mutate(date = as.Date(paste0(year, "-", month, "-", "01"))) 

# 1b Looking at NAs in precipitation - all from these NAs in the exact raster square
precip_dat %>% 
  filter(is.na(precip) == T) %>% 
  summarise(n_NA_precip = n(),studies = n_distinct(ID),
            precip50m_NAs = length(which(is.na(precip_50m) == T)),
            precip5km_NAs = length(which(is.na(precip_5km) == T)),
            precip50km_NAs = length(which(is.na(precip_50km) == T)))

world_sf <- ne_coastline(scale = "medium", returnclass = "sf")

pNA <- precip_dat %>% 
  filter(is.na(precip) == T) %>%
  group_by(ID) %>% 
  summarise_all(.funs = first) 

# All just by the coast
ggplot(data = world_sf) +
  geom_sf(size = 0.4) +
  geom_point(data = pNA, colour = "red", size = 2,
             aes(x = Longitude, y = Latitude)) 

# ID 2357 and 4920 - Seems to be that whole time-series are lost i.e. rasters line up
filter(precip_dat, ID == 2357)$precip_5km
filter(precip_dat, ID == 4920)$precip_5km

##__________________________________________________________________________________________________
#### 2. Additive seasonal Decomposition by Loess (STL) - Cheetah Temp Test ####

# 2a. Have a look at the temperature test
cheetah_test <- ggplot(temp_test, aes(x = date, y = temp)) +
  geom_line(colour = "firebrick", size = 0.2) +
  labs(x = NULL, y = "Scaled temperature") +
  theme_bw(base_size = 9)

ggsave(cheetah_test,
       filename = "plots/developing_annual_weather_variables/temp_scale_cheetah_example.jpeg",
       width = 5, height = 3, units = "in", dpi = 400)

# 2b. Convert to a timeseries
test_ts <- ts(zoo(temp_test$temp, order.by= temp_test$date),
               frequency=12, start=c(1979,1))

# 2c. STL - decomposition - 
# Magnitude of components are relatively robust to the changes in s and t windows.
test_stl <- stl(test_ts, s.window=7, t.window = 1000)

# 2d. Have a look
jpeg("plots/developing_annual_weather_variables/weather_decomposition_example_cheetah.jpeg",
     width = 5, height = 7, units = "in", res = 400)
plot(test_stl)
dev.off()

##__________________________________________________________________________________________________
#### 3. STL decomposition for each study ID ####

# 3a. using group_modify to do a function per study ID - very neat
temp_stl <- temp_dat %>% 
  pivot_longer(-c(ID,Binomial,Longitude, Latitude, year, month, date), 
               names_to = "temp_scale", values_to = "temp") %>% 
  group_by(ID, temp_scale) %>% 
  mutate(temp = as.numeric(scale(temp))) %>% 
  ungroup() %>% 
  group_by(ID, temp_scale) %>% 
  group_modify(~ {
    crr_ts = ts(zoo(.$temp, order.by= .$date),
                frequency=12, start=c(1979,1))
    crr_stl = stl(crr_ts, s.window=7, t.window = 1000)
    mutate(., 
           temp_season = crr_stl$time.series[,1],
           temp_trend = crr_stl$time.series[,2],
           temp_anomaly = crr_stl$time.series[,3],
           temp_season_range = diff(range(crr_stl$time.series[,1])),
           temp_anomaly_range = diff(range(crr_stl$time.series[,3])))}) %>% 
  ungroup()

precip_stl <- precip_dat %>% 
  pivot_longer(-c(ID,Binomial,Longitude, Latitude, year, month, date), 
               names_to = "precip_scale", values_to = "precip") %>% 
  group_by(ID, precip_scale) %>% 
  mutate(precip = as.numeric(scale(precip))) %>% 
  ungroup() %>% 
  filter(is.na(precip) == FALSE) %>% 
  group_by(ID, precip_scale) %>% 
  group_modify(~ {
    crr_ts = ts(zoo(.$precip, order.by= .$date),
                frequency=12, start=c(1979,1))
    crr_stl = stl(crr_ts, s.window=7, t.window = 1000)
    mutate(., 
           precip_season = crr_stl$time.series[,1],
           precip_trend = crr_stl$time.series[,2],
           precip_anomaly = crr_stl$time.series[,3],
           precip_season_range = diff(range(crr_stl$time.series[,1])),
           precip_anomaly_range = diff(range(crr_stl$time.series[,3])))}) %>% 
  ungroup()

## Join the data together
p_join <- dplyr::select(precip_stl,
                        ID, year, month, precip_scale, 
                        precip, precip_season, precip_trend, 
                        precip_anomaly, precip_season_range, 
                        precip_anomaly_range) %>% 
  mutate(scale = gsub(pattern = "precip", replacement = "scale",
                      x = precip_scale))

# This is a big file - can this be done more efficiently?
chelsa_stl <- temp_stl %>% 
  mutate(scale = gsub(pattern = "temp", replacement = "scale",
                      x = temp_scale)) %>% 
  left_join(x = ., y = p_join, by = c("ID", "year", "month", "scale"))

##__________________________________________________________________________________________________
#### 4. Plotting out decompositions for a random sample of studies ####

# Just with the exact raster values
set.seed(10)
sample_IDs <- sample(x = unique(temp_stl$ID), size = 10)

chelsa_stl_fp <- chelsa_stl %>% 
  filter(scale == "scale") %>% 
  mutate(Study = paste0(ID, ": ",gsub(pattern = "_", replacement = " ",
                         x = Binomial)))

# iterate through, pull out data and plot
for(i in sample_IDs){
  # the right data and some wrangling
  cdat_temp = filter(chelsa_stl_fp, ID == i) %>% 
    pivot_longer(cols = c(temp, temp_season, temp_trend, temp_anomaly),
                 names_to = "component", values_to = "value") %>% 
    mutate(component = factor(component, levels = c("temp", "temp_season", 
                                                    "temp_trend", "temp_anomaly")))
  cdat_precip = filter(chelsa_stl_fp, ID == i) %>% 
    pivot_longer(cols = c(precip, precip_season, precip_trend, precip_anomaly),
                 names_to = "component", values_to = "value") %>% 
    mutate(component = factor(component, levels = c("precip", "precip_season", 
                                                    "precip_trend", "precip_anomaly")))
  
  study = cdat_temp$Study[1]
  
  # the plot 
  cplot_temp = ggplot(cdat_temp, aes(x = date, y = value)) +
    geom_line(colour = "firebrick") +
    facet_wrap(~component,
               ncol = 1) +
    labs(x = NULL, y = "Value", title = study) +
    theme_bw(base_size = 17) +
    theme(strip.background = element_blank())
  
  cplot_precip = ggplot(cdat_precip, aes(x = date, y = value)) +
    geom_line(colour = "deepskyblue") +
    facet_wrap(~component,
               ncol = 1) +
    labs(x = NULL, y = "Value", title = study) +
    theme_bw(base_size = 17) +
    theme(strip.background = element_blank())
  
  # assign
  assign(x = paste0("temp_", i), value = cplot_temp)
  assign(x = paste0("precip_", i), value = cplot_precip)
}

ggsave(grid.arrange(temp_11512, temp_11754, temp_18256, temp_18286, temp_3443,
                    temp_5088, temp_5595, temp_5850, temp_6394, temp_8147, 
                    ncol = 5),
       filename = "plots/developing_annual_weather_variables/temp_stl_sample.jpeg", 
       width = 25, height = 22, units = "in", dpi = 400)

ggsave(grid.arrange(precip_11512, precip_11754, precip_18256, precip_18286, precip_3443,
                    precip_5088, precip_5595, precip_5850, precip_6394, precip_8147, 
                    ncol = 5),
       filename = "plots/developing_annual_weather_variables/precip_stl_sample.jpeg", 
       width = 25, height = 22, units = "in", dpi = 400)

# tidy the environment
rm(temp_11512, temp_11754, temp_18256, temp_18286, temp_3443,
   temp_5088, temp_5595, temp_5850, temp_6394, temp_8147,
   precip_11512, precip_11754, precip_18256, precip_18286, precip_3443,
   precip_5088, precip_5595, precip_5850, precip_6394, precip_8147)

##__________________________________________________________________________________________________
#### 4. Exploring the range of anomaly vs. season ####

# Want to have a look to see if the anomaly components that have been extracted 
# are substantial relative to the season component. 

temp_stl_range <- temp_stl %>% 
  group_by(ID,temp_scale) %>% 
  summarise(temp_season_range = temp_season_range[1],
            temp_anomaly_range = temp_anomaly_range[1])

precip_stl_range <- precip_stl %>% 
  group_by(ID,precip_scale) %>% 
  summarise(precip_season_range = precip_season_range[1],
            precip_anomaly_range = precip_anomaly_range[1])


ggplot(temp_stl_range, aes(x = temp_anomaly_range, y = temp_season_range)) + 
  geom_abline(slope = 1, intercept = 0) +
  geom_point(alpha = 0.2, size = 3, colour = "firebrick") +
  labs(x = "Range of temperature anomaly component",
       y = "Range of temperature seasonal component") +
  facet_wrap(~ temp_scale) +
  theme_bw(base_size = 17) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  ggsave(filename = "plots/developing_annual_weather_variables/temp_anomaly_season_scales.jpeg",
         width = 10, height = 10, units = "in", dpi = 400)

ggplot(precip_stl_range, aes(x = precip_anomaly_range, y = precip_season_range)) + 
  geom_abline(slope = 1, intercept = 0) +
  geom_point(alpha = 0.2, size = 3, colour = "deepskyblue") +
  labs(x = "Range of precipitation anomaly component",
       y = "Range of precipitation seasonal component") +
  facet_wrap(~ precip_scale) +
  theme_bw(base_size = 17) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  ggsave(filename = "plots/developing_annual_weather_variables/precip_anomaly_season_scales.jpeg",
         width = 10, height = 10, units = "in", dpi = 400)
  
# # check
# ID28_temp <- filter(temp_stl,  ID == 28 & temp_scale == "temp")

##__________________________________________________________________________________________________
#### 5. Save ####

saveRDS(chelsa_stl, "data/chelsa_stl.RDS") # add in precipitation when ready









