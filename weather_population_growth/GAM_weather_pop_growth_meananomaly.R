####################################################
##                                                ##
##     Global climate and population dynamics     ##
##                                                ##
##           GAMM timeseries method               ##
##                                                ##
##                Aug 9th 2020                    ##
##                                                ##
####################################################

# Record-wise regressions linking weather to population growth rates, 
# accounting for autocorrelation with GAMM and formal AR(1) time-series analysis
# Models just for 5km buffer radius and the mean weather anomaly.

rm(list = ls())
options(width = 100)

library(tidyverse)
library(ggridges)
library(viridis)
library(patchwork)
library(gridExtra)
library(mgcv)      # Simon Wood to the rescue again. All Hail

temp_colour <- "#990a80"
precip_colour <- "#287f79"

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
  geom_hline(yintercept = 0, linetype = "solid", size = 0.2) +
  geom_segment(aes(xend = year, yend = pop_growth_rate + mean_temp_anomaly),
               colour = temp_colour, arrow = arrow(length = unit(0.1, "cm"))) +
  # geom_segment(aes(xend = year, yend = pop_growth_rate + mean_precip_anomaly),
  #              colour = precip_colour, arrow = arrow(length = unit(0.1, "cm"))) +
  geom_point(size = 3) + geom_line() +
  geom_smooth(se = F, colour = "black", size = 0.5) +
  labs(x = "Year", y = "Population growth rate") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank())


ggplot(wtd, aes(x = mean_temp_anomaly, y = pop_growth_rate)) + 
  geom_hline(yintercept = 0, linetype = "solid", size = 0.2) +
  geom_point(size = 3) +
  labs(x = "Temperature anomaly", y = "Population growth rate") +
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

# Correlation for plot labels
temp_cor <- cor.test(lin_gam$lin_temp, lin_gam$coef_temp)
temp_cor_val <-  c(round(temp_cor$estimate, 2),
                   ifelse(temp_cor$p.value < 0.001, "p < 0.001",
                          paste0("p = ", round(temp_cor$p.value, 2))))

precip_cor <- cor.test(lin_gam$lin_precip, lin_gam$coef_precip)
precip_cor_val <- c(round(precip_cor$estimate, 2),
                    ifelse(precip_cor$p.value < 0.001, "p < 0.001",
                           paste0("p = ", round(precip_cor$p.value, 2))))

temp_compare <- ggplot(lin_gam, aes(x = lin_temp, y = coef_temp, size = n_obs)) + 
  geom_point(alpha = 0.6, colour = temp_colour) +
  geom_abline(slope = 1, intercept = 0) +
  annotate('text', x = -5, y = 7,
           label = paste0("r = ", temp_cor_val[1], ", ", temp_cor_val[2])) +
  scale_size_continuous(range = c(1,8), guide = F) +
  labs(x = "Linear temperature effect", y = "GAM temperature effect") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank())

precip_compare <- ggplot(lin_gam, aes(x = lin_precip, y = coef_precip, size = n_obs)) + 
  geom_point(alpha = 0.6, colour = precip_colour) +
  geom_abline(slope = 1, intercept = 0) +
  annotate('text', x = -0.5, y = 5,
           label = paste0("r = ", precip_cor_val[1], ", ", precip_cor_val[2])) +
  scale_size_continuous(range = c(1,8), guide = F) +
  labs(x = "Linear precipitation effect", y = "GAM precipitation effect") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank())

ggsave(temp_compare + precip_compare,
       filename = "plots/weather_pop_growth/linear_gam_comparison.jpeg",
       width = 20, height = 10, units = "cm", dpi = 400)

##__________________________________________________________________________________________________
#### 6. Saving data from GAM ####

mnanom_5km_GAM <- pgr_weather_gam
save(mnanom_5km_GAM, file = "data/pgr_weather/mnanom_5km_GAM.RData")

##__________________________________________________________________________________________________
#### 7. Hypothesis exploration plots ####

## 7a. Coefficients overall
coef_ovrdat <- pgr_weather_gam %>% 
  dplyr::select(ID_block, starts_with("coef_")) %>% 
  pivot_longer(-ID_block) %>% 
  mutate(model = if_else(name == "coef_temp", "Temperature", "Precipitation"),
         label = if_else(name == "coef_temp", 
                         "Mean temperature anomaly",
                         "Mean precipitation anomaly")) 

ggplot(coef_ovrdat, aes(x = value, y = label, fill = model)) +
  geom_vline(xintercept = 0) +
  geom_density_ridges(alpha = 0.7, scale = 1, size = 0.3) +
  coord_cartesian(xlim = c(-5,5)) +
  scale_fill_manual(guide = F, values = c(precip_colour, temp_colour)) +
  labs(x = "Model coefficient", y = NULL) +
  theme_ridges(center_axis_labels = TRUE, font_size = 16) +
  ggsave(filename = "plots/weather_pop_growth/overall_coefficients_mnanom_5km_GAM.jpeg",
         height = 5, width = 6, units = "in", dpi = 400)

## 7b. Weather coefficients by Order
temp_sp <- ggplot(pgr_weather_gam, aes(x = coef_temp, y = Order, 
                                   fill = Order, height = stat(density))) +
  geom_vline(xintercept = 0) +
  geom_density_ridges(alpha = 0.7, scale = 2, stat = "density", size = 0.3) +
  coord_cartesian(xlim = c(-5,5)) +
  scale_fill_viridis_d(guide = F, option = "C") +
  labs(x = "Temperature anomaly coefficient", y = NULL) +
  theme_ridges(center_axis_labels = TRUE, font_size = 25) 

precip_sp <- ggplot(pgr_weather_gam, aes(x = coef_precip, y = Order, 
                                     fill = Order, height = stat(density))) +
  geom_vline(xintercept = 0) +
  geom_density_ridges(alpha = 0.7, scale = 1.5, stat = "density", size = 0.3) +
  coord_cartesian(xlim = c(-5,5)) +
  scale_fill_viridis_d(guide = F, option = "D") +
  labs(x = "Precipitation anomaly coefficient", y = NULL) +
  theme_ridges(center_axis_labels = TRUE, font_size = 25) +
  theme(axis.text.y = element_blank())

ggsave(temp_sp + precip_sp,
       filename = "plots/weather_pop_growth/coef_order_mnanom_5km_GAM.jpeg",
       width = 15, height = 13, units = "in", dpi = 400)

## 7c. Weather coefficients by biome
temp_biome <- ggplot(pgr_weather_gam, aes(x = coef_temp, y = biome, fill = biome,
                                      height = stat(density))) +
  geom_vline(xintercept = 0) +
  geom_density_ridges(alpha = 0.7, scale = 2, stat = "density", size = 0.3) +
  coord_cartesian(xlim = c(-5,5)) +
  scale_fill_viridis_d(guide = F, option = "C") +
  labs(x = "Temperature anomaly coefficient", y = NULL) +
  theme_ridges(center_axis_labels = TRUE, font_size = 20) 

precip_biome <- ggplot(pgr_weather_gam, aes(x = coef_precip, y = biome, fill = biome,
                                        height = stat(density))) +
  geom_vline(xintercept = 0) +
  geom_density_ridges(alpha = 0.7, scale = 1.5, stat = "density", size = 0.3) +
  coord_cartesian(xlim = c(-5,5)) +
  scale_fill_viridis_d(guide = F, option = "D") +
  labs(x = "Precipitation anomaly coefficient", y = NULL) +
  theme_ridges(center_axis_labels = TRUE, font_size = 20) +
  theme(axis.text.y = element_blank())

ggsave(temp_biome + precip_biome,
       filename = "plots/weather_pop_growth/coef_biome_mnanom_5km_GAM.jpeg",
       width = 45, height = 20, units = "cm", dpi = 400)

## 7d. Weather coefficients by latitude
pgr_lat <- pgr_weather_gam %>% 
  mutate(lat = abs(Latitude) - (abs(Latitude) %% 22.5))

temp_lat <- ggplot(pgr_lat, aes(x = coef_temp, y = factor(lat), 
                                fill = factor(lat), height = stat(density))) +
  geom_vline(xintercept = 0) +
  stat_density_ridges(quantile_lines = T, quantiles = 2, alpha = 0.7, scale = 1.1) +
  coord_cartesian(xlim = c(-5,5)) +
  scale_fill_viridis_d(guide = F, option = "C", begin = 0.1, end = 0.9) + 
  labs(x = "Temperature anomaly coefficient", y = "Absolute latitude") +
  theme_ridges(center_axis_labels = TRUE, font_size = 25) 

precip_lat <- ggplot(pgr_lat, aes(x = coef_precip, y = factor(lat), 
                                  fill = factor(lat), height = stat(density))) +
  geom_vline(xintercept = 0) +
  stat_density_ridges(quantile_lines = T, quantiles = 2, alpha = 0.7, scale = 1.1) +
  coord_cartesian(xlim = c(-5,5)) +
  scale_fill_viridis_d(guide = F, option = "D", begin = 0.2, end = 0.8) + 
  labs(x = "Precipitation anomaly coefficient", y = NULL) +
  theme_ridges(center_axis_labels = TRUE, font_size = 25) 

ggsave(grid.arrange(temp_lat, precip_lat, ncol = 2, widths = c(7,6)),
       filename = "plots/weather_pop_growth/coef_lat_mnanom_5km_GAM.jpeg",
       width = 15, height = 11, units = "in", dpi = 400)


