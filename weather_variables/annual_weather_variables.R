#####################################################
##                                                 ##
##     Global climate and population dynamics      ##
##                                                 ##
##        Annual CHELSA weather measures           ##
##                                                 ##
##               April 20th 2020                   ##
##                                                 ##
#####################################################
rm(list = ls())
options(width = 100)

library(tidyverse)
library(grid)
library(gridExtra)

temp_stl <- readRDS("data/temp_stl.RDS") # again testing with temperature for now

##__________________________________________________________________________________________________
#### 1. Annual weather measures ####

mam_chelsa_annual <- temp_stl %>% 
  group_by(ID, temp_scale) %>% 
  mutate(ts_anomaly_sd = sd(anomaly),
         ts_anomaly_mean = mean(anomaly)) %>% # adding summary stats for (6) for the full time series
  ungroup() %>% 
  group_by(ID, year, temp_scale) %>% 
  summarise(Binomial = Binomial[1], 
            Longitude = Longitude[1], 
            Latitude = Latitude[1],
            
            # 1) Mean climate variable
            mean_temp = mean(temp),
            
            # 2) Mean anomaly
            mean_temp_anomaly = mean(anomaly),
            
            # 3) Mean absolute anomaly
            mean_abtemp_anomaly = sd(abs(anomaly)),
            
            # 4) Maximum/minimum anomaly
            max_temp_anomaly = max(anomaly),
            min_temp_anomaly = min(anomaly),
            
            # 5) Anomaly variance
            anomaly_variance = var(anomaly),

            # 6) no. odd months
            no_odd_months_temp = length(which(abs(anomaly) > 
                                                ts_anomaly_mean + (2*ts_anomaly_sd))))

# test to demonstrate (6)
ts_ID28 <- filter(temp_stl, ID == 28 & temp_scale == "temp")
ID28 <- filter(temp_stl, ID == 28 & temp_scale == "temp", year == 1979)

ts_an_mean <- mean(ts_ID28$anomaly)
ts_an_sd <- sd(ts_ID28$anomaly)

ggplot(ts_ID28, aes(x = date)) +
  geom_segment(aes(xend = date, y = 0, yend = anomaly)) +
  geom_hline(yintercept = ts_an_mean + 2*ts_an_sd, colour = "red", 
             linetype = "dashed") +
  geom_hline(yintercept = ts_an_mean - 2*ts_an_sd, colour = "red", 
           linetype = "dashed") +
  labs(x = NULL, y = "Temperature anomaly", title = "Record 28: Cheetah") +
  theme_bw(base_size = 17) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggsave(filename = "plots/mam_chelsa/ID28_anomaly_sd.jpeg",
         width = 7, height = 4, units = "in", dpi = 400)
  

##__________________________________________________________________________________________________
#### 2. Save ####

saveRDS(mam_chelsa_annual, file = "data/mam_chelsa_annual.RDS") # <---- need to add precip















