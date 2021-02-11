####################################################
##                                                ##
##     Global climate and population dynamics     ##
##                                                ##
##           GAMM timeseries method               ##
##                                                ##
##                Dec 8th 2020                    ##
##                                                ##
####################################################

# Record-wise regressions linking weather to population growth rates, 
# accounting for autocorrelation with GAMM and formal AR(1) time-series analysis

rm(list = ls())
options(width = 100)

library(tidyverse)
library(viridis)
library(patchwork)
library(mgcv)      # Simon Wood to the rescue again. All Hail

temp_colour <- viridis(20, option = "C")[13]
precip_colour <- viridis(20, option = "D")[10]

##__________________________________________________________________________________________________
#### 1. Load data ####

# mammal data
load("../rawdata/mammal.RData")
glimpse(mammal)

# annual weather anomaly - focus on just the mean anomaly in this script at a 5km range
mam_chelsa_annual <- readRDS("data/mam_chelsa_annual.RDS") %>% 
  filter(scale == "scale_5km") %>% 
  dplyr::select(ID,year, weather_scale = scale, mean_temp_anomaly, mean_precip_anomaly)
glimpse(mam_chelsa_annual)

# Species names to merge
load("../rawdata/GBIF_species_names_mamUPDATE.RData", verbose = TRUE)

##__________________________________________________________________________________________________
#### 2. Joining data ####

# linking to weather data and species names 
mammal_weather <- mammal %>% 
  left_join(., y = mam_chelsa_annual, by = c("ID", "year")) %>% 
  left_join(., y = dplyr::select(lpd_gbif, Binomial, gbif_species = gbif.species.tree),
            by = "Binomial") %>% 
  mutate(year_s = as.numeric(scale(year)))

glimpse(mammal_weather)

##__________________________________________________________________________________________________
#### 3. White tailed-deer example ####

wtd <- filter(mammal_weather, ID_block == "22039_1")

ggplot(wtd, aes(x = year, y = pop_growth_rate)) + 
  geom_hline(yintercept = 1, linetype = "solid", size = 0.2) +
  geom_segment(aes(xend = year, yend = pop_growth_rate + mean_temp_anomaly),
               colour = temp_colour, arrow = arrow(length = unit(0.1, "cm"))) +
  # geom_segment(aes(xend = year, yend = pop_growth_rate + mean_precip_anomaly),
  #              colour = precip_colour, arrow = arrow(length = unit(0.1, "cm"))) +
  geom_point(size = 3) + geom_line() +
  geom_smooth(se = F, colour = "black", size = 0.5) +
  labs(x = "Year", y = "Per-capita population growth rate") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank())


ggplot(wtd, aes(x = mean_temp_anomaly, y = pop_growth_rate)) + 
  geom_hline(yintercept = 1, linetype = "solid", size = 0.2) +
  geom_point(size = 3) +
  labs(x = "Temperature anomaly", y = "Per-capita population growth rate") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank())

#____________________________________
## Temperature

wtd_gamarma <- gamm(pop_growth_rate ~ s(year, bs = "tp") + mean_temp_anomaly,
                  data = wtd, family = gaussian,
                  correlation = corARMA(form = ~ 1 | year, p = 1),
                  method = "REML")

wtd_linear <- lm(pop_growth_rate ~ mean_temp_anomaly + ln_abundance + year, data = wtd)
wtd_naive <- lm(pop_growth_rate ~ mean_temp_anomaly, data = wtd)

## The GAMM looks to be doing a similar job, whilst formally accounting for temporal autocorrelation

##__________________________________________________________________________________________________
#### 4. Models for each record ####

pgr_weather_gam <- mammal_weather %>% 
  group_by(ID_block) %>% 
  group_modify(~{
    
    # Temperature
    mod_temp = gamm(pop_growth_rate ~ s(year, bs = "tp", k = 5) + mean_temp_anomaly,
                    data = ., family = gaussian,
                    correlation = corARMA(form = ~ year, p = 1),
                    method = "REML")
    coef_tempmod = coef(mod_temp$gam)
    
    # Precipitation + dealing with NA values
    if(length(which(is.na(.$mean_precip_anomaly) == T)) == 0){
      mod_precip = gamm(pop_growth_rate ~ s(year, bs = "tp", k = 5) + mean_precip_anomaly,
                        data = ., family = gaussian,
                        correlation = corARMA(form = ~ year, p = 1),
                        method = "REML")
      coef_precipmod = coef(mod_precip$gam)}
    else{coef_precipmod = rep(NA,15)}     # Arbitrary long NA vector
    
    tibble(.[1,],
           coef_temp = unname(coef_tempmod[2]),
           coef_precip = unname(coef_precipmod[2]),
           n_obs = nrow(.))
  }) 

glimpse(pgr_weather_gam)

##__________________________________________________________________________________________________
#### 5. Comparison to linear effects ####

load("data/pgr_weather/mnanom_5km.RData")

lin_gam <- pgr_weather_gam %>% 
  dplyr::select(ID_block, Binomial, n_obs, coef_temp, coef_precip) %>% 
  left_join(x =., y = dplyr::select(mnanom_5km, ID_block, 
                                    lin_temp = coef_temp, lin_precip = coef_precip),
            by = "ID_block")
  

temp_compare <- ggplot(lin_gam, aes(x = lin_temp, y = coef_temp, size = n_obs)) + 
  geom_point(alpha = 0.6, colour = temp_colour) +
  geom_abline(slope = 1, intercept = 0) +
  scale_size_continuous(range = c(1,8), guide = F) +
  labs(x = "Linear temperature effect", y = "GAM temperature effect") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank())

precip_compare <- ggplot(lin_gam, aes(x = lin_precip, y = coef_precip, size = n_obs)) + 
  geom_point(alpha = 0.6, colour = precip_colour) +
  geom_abline(slope = 1, intercept = 0) +
  scale_size_continuous(range = c(1,8), guide = F) +
  labs(x = "Linear precipitation effect", y = "GAM precipitation effect") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank())

ggsave(temp_compare + precip_compare,
       filename = "plots/weather_pop_growth/linear_gam_comparison.jpeg",
       width = 25, height = 18, units = "cm", dpi = 400)

##__________________________________________________________________________________________________
#### 6. Saving data from GAM ####

mnanom_5km_GAM <- pgr_weather_gam
save(mnanom_5km_GAM, file = "data/pgr_weather/mnanom_5km_GAM.RData")



