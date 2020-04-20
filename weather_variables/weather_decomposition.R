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
library(ggfortify)

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

## All temperature data
temp_dat <- mam_chelsa %>% 
  dplyr::select(ID, Binomial, Longitude, 
                Latitude, year, month, 
                temp, temp_50m, 
                temp_5km, temp_50km) %>% 
  mutate(date = as.Date(paste0(year, "-", month, "-", "01"))) # need to transform weather variables

##__________________________________________________________________________________________________
#### 2. Additive seasonal Decomposition by Loess (STL) - Cheetah Test ####

# 2a. Have a look at the temperature test
ggplot(temp_test, aes(x = date, y = temp)) +
  geom_line(colour = "firebrick") +
  labs(x = NULL, y = "Scaled temperature") +
  theme_bw(base_size = 16)

# 2b. Convert to a timeseries
test_ts <- ts(zoo(temp_test$temp, order.by= temp_test$date),
               frequency=12, start=c(1979,1))

# 2c. STL - decomposition - 
# Magnitude of components are relatively robust to the changes in s and t windows.
test_stl <- stl(test_ts, s.window=7, t.window = 1000)

# 2d. Have a look
plot(test_stl)

##__________________________________________________________________________________________________
#### 3. STL decomposition for each study ID ####

# 3a. using group_modify to do a function per study ID - very neat <---- WHEN ADDING PRECIP, TRY TO DO IT LONG FORMAT with one scale column
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
           season = crr_stl$time.series[,1],
           trend = crr_stl$time.series[,2],
           anomaly = crr_stl$time.series[,3],
           season_range = diff(range(crr_stl$time.series[,1])),
           anomaly_range = diff(range(crr_stl$time.series[,3])))}) %>% 
  ungroup()

# 3b. Plotting out decompositions for a random sample of studies
# Just for temp for now
set.seed(10)
sample_IDs <- sample(x = unique(temp_stl$ID), size = 10)

temp_stl_temp <- temp_stl %>% 
  filter(temp_scale == "temp") %>% 
  mutate(Study = paste0(ID, ": ",gsub(pattern = "_", replacement = " ",
                         x = Binomial)))

# iterate through, pull out data and plot
for(i in sample_IDs){
  # the right data and some wrangling
  cdat = filter(temp_stl_temp, ID == i) %>% 
    pivot_longer(cols = c(temp,season,trend,anomaly),
                 names_to = "component", values_to = "value") %>% 
    mutate(component = factor(component, levels = c("temp", "season", "trend", "anomaly")))
  
  study = cdat$Study[1]
  
  # the plot 
  cplot = ggplot(cdat, aes(x = date, y = value)) +
    geom_line(colour = "firebrick") +
    facet_wrap(~component,
               ncol = 1) +
    labs(x = NULL, y = "Value", title = study) +
    theme_bw(base_size = 17) +
    theme(strip.background = element_blank())
  
  # assign
  assign(x = paste0("stl_", i), value = cplot)
}

ggsave(grid.arrange(stl_11512, stl_11754, stl_18256, stl_18286, stl_3443,
                    stl_5088, stl_5595, stl_5850, stl_6394, stl_8147, 
                    ncol = 5),
       filename = "plots/mam_chelsa/temp_stl_sample.jpeg", 
       width = 25, height = 22, units = "in", dpi = 400)

# tidy the environment
rm(stl_11512, stl_11754, stl_18256, stl_18286, stl_3443,
   stl_5088, stl_5595, stl_5850, stl_6394, stl_8147)

##__________________________________________________________________________________________________
#### 4. Exploring the range of anomaly vs. season ####

# Want to have a look to see if the anomaly components that have been extracted 
# are substantial relative to the season component. 

temp_stl_range <- temp_stl %>% 
  group_by(ID,temp_scale) %>% 
  summarise(season_range = season_range[1],
            anomaly_range = anomaly_range[1])

ggplot(temp_stl_range, aes(x = anomaly_range, y = season_range)) + 
  geom_abline(slope = 1, intercept = 0) +
  geom_point(alpha = 0.2, size = 3) +
  labs(x = "Range of temperature anomaly component",
       y = "Range of temperature seasonal component") +
  facet_wrap(~ temp_scale) +
  theme_bw(base_size = 17) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  ggsave(filename = "plots/mam_chelsa/temp_anomaly_season_scales.jpeg",
         width = 10, height = 10, units = "in", dpi = 400)
  
# # check
# ID28_temp <- filter(temp_stl,  ID == 28 & temp_scale == "temp")

##__________________________________________________________________________________________________
#### 5. Save ####

saveRDS(temp_stl, "data/temp_stl.RDS") # add in precipitation when ready









