####################################################
##                                                ##
##     Global climate and population dynamics     ##
##                                                ##
##      Annual weather and population growth      ##
##                                                ##
##               June 11th 2020                   ##
##                                                ##
####################################################
rm(list = ls())
options(width = 100)

library(tidyverse)
library(ggridges)
library(viridis)
library(grid)
library(gridExtra)

##__________________________________________________________________________________________________
#### 1. Load data ####

# mammal data
load("../rawdata/mammal.RData")
glimpse(mammal)

# annual weather anomaly - focus on just the mean anomaly in this script at a 5km range
mam_chelsa_annual <- readRDS("data/mam_chelsa_annual.RDS") %>% 
  dplyr::select(-c(4:6))
glimpse(mam_chelsa_annual)

##__________________________________________________________________________________________________
#### 2. Joining data ####

mammal_weather <- mammal %>% 
  left_join(., y = mam_chelsa_annual, by = c("ID", "year"))

##__________________________________________________________________________________________________
#### 3. Linear models for each variable and scale for each record ####

# 3a. set up iteration data
# Ignoring number of odd days vars for now - they follow a zero inflated pattern
iter_dat <- expand_grid(ID_block = unique(mammal_weather$ID_block),
                               scale = unique(mammal_weather$scale),
                               weather_var = colnames(mammal_weather)[24:39])

# 3b. weather coefficients for each variable
pgr_weather_res <- bind_rows(lapply(X = 1:nrow(iter_dat), function(x){
  
  crow = iter_dat[x,]
  
  # current data
  cdat = mammal_weather %>% 
    filter(ID_block == crow$ID_block, scale == crow$scale) %>% 
    dplyr::select(ID_block, year, ln_abundance,
                  weather_val = crow$weather_var,
                  pop_growth_rate)
  
  # record info
  rec_info = mammal_weather %>% 
    filter(ID_block == crow$ID_block, scale == crow$scale) %>% 
    dplyr::select(2:17) %>% 
    slice(1)
  
  # model
  if(length(which(is.na(cdat$weather_val) == T)) > 0){modcoef = rep(NA,4)}
  else{mod_weather = lm(pop_growth_rate ~ weather_val + ln_abundance + year, data = cdat)
       modcoef = coefficients(mod_weather)}
  
  # returning data
  cat('\r',"Your Job is",round((x/nrow(iter_dat))*100, 0),"% Complete       ")
  return(tibble(crow, coef_weather = modcoef[2], 
                coef_abun = modcoef[3], coef_trend = modcoef[4],
                rec_info))
}))
  
# 3c. Adding in weather variable labels
pgr_weather_res <- pgr_weather_res %>% 
  mutate(weather_var_lab = stringr::str_to_sentence(gsub("_", " ", weather_var))) %>% 
  mutate(weather_var_lab = gsub("emp", "emperature", weather_var_lab),
         weather_var_lab = gsub("recip", "recipitation", weather_var_lab))

##__________________________________________________________________________________________________
#### 4. Density rigge plots for the weather variables ####

# removing very large coefficients
pgr_plotdat_sm <- pgr_weather_res %>% 
  filter(coef_weather >= -0.05 & coef_weather < 0.05)

pgr_plotdat_lg <- pgr_weather_res %>% 
  filter(coef_weather >= -5 & coef_weather < 5)

pgr_weath_sm <- ggplot(pgr_plotdat_sm, aes(x = coef_weather, y = weather_var_lab, fill = stat(x))) +
  geom_vline(xintercept = 0) +
  geom_density_ridges_gradient(scale = 1.1) +
  scale_fill_viridis_c(option = "D", guide = F) +
  labs(x = "Weather variable coefficient", y = NULL) +
  theme_ridges(center_axis_labels = TRUE, font_size = 20, grid = F)

pgr_weath_lg <- ggplot(pgr_plotdat_lg, aes(x = coef_weather, y = weather_var_lab, fill = stat(x))) +
  geom_vline(xintercept = 0) +
  geom_density_ridges_gradient(scale = 1.1, ) +
  scale_fill_viridis_c(option = "D", guide = F) +
  labs(x = "Weather variable coefficient", y = NULL) +
  theme_ridges(center_axis_labels = TRUE, font_size = 20, grid = F) +
  theme(axis.text.y = element_blank())

ggsave(grid.arrange(pgr_weath_sm, pgr_weath_lg, ncol = 2, widths = c(6,4)),
       filename = "plots/weather_pop_growth/coef_weather_vars.jpeg",
       width = 15, height = 11, units = "in", dpi = 400)

##__________________________________________________________________________________________________
#### 5. Save data ####

saveRDS(pgr_weather_res, file = "data/pgr_weather/pgr_weather_res.RDS")





