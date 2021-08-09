####################################################
##                                                ##
##     Global climate and population dynamics     ##
##                                                ##
##      Annual weather and population growth      ##
##                                                ##
##               March 16th 2021                  ##
##                                                ##
####################################################

# Record-wise regressions linking weather to population growth rates, 
# accounting for autocorrelation with GAMM and formal AR(1) time-series analysis.
# Models across all weather variables and spatial scales.

rm(list = ls())
options(width = 100)

library(tidyverse)
library(psych)
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
  dplyr::select(-c(4:6))
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
#### 3. GAMs for each variable and scale for each record ####

# 3a. set up iteration data
# Ignoring number of odd days vars for now - they follow a zero inflated pattern
iter_dat <- expand_grid(ID_block = unique(mammal_weather$ID_block),
                               scale = unique(mammal_weather$scale),
                               weather_var = colnames(mammal_weather)[25:40])

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
  else{mod_weather = gamm(pop_growth_rate ~ 
                            s(year, bs = "tp", k = 5) + weather_val,
                          data = cdat, 
                          family = gaussian,
                          correlation = corARMA(form = ~ year, p = 1),
                          method = "REML")
       modcoef = coef(mod_weather$gam)}
  
  # returning data
  cat('\r',"Your Job is",round((x/nrow(iter_dat))*100, 0),"% Complete       ")
  return(tibble(crow, coef_weather = modcoef[2], 
                rec_info))
}))
  
# 3c. Adding in weather variable labels
pgr_weather_res <- pgr_weather_res %>% 
  mutate(weather_var_lab = stringr::str_to_sentence(gsub("_", " ", weather_var))) %>% 
  mutate(weather_var_lab = gsub("emp", "emperature", weather_var_lab),
         weather_var_lab = gsub("recip", "recipitation", weather_var_lab))

##__________________________________________________________________________________________________
#### 4. Density ridge plots for the weather variables ####

# pgr_weather_res <- readRDS("data/pgr_weather/pgr_weather_res.RDS")

# removing very large coefficients
pgr_plotdat_sm <- pgr_weather_res %>% 
  filter(coef_weather >= -0.2 & coef_weather < 0.2)

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

ggsave(pgr_weath_sm,
       filename = "plots/weather_pop_growth/coef_weather_vars.jpeg",
       width = 8, height = 11, units = "in", dpi = 400)

##__________________________________________________________________________________________________
#### 5. Spatial scales ####

# 5a. Spatial scales consistent?
sp_res <- pgr_weather_res %>% 
  dplyr::select(ID_block, scale, coef_weather, weather_var) %>% 
  pivot_wider(names_from = scale, values_from = coef_weather) %>% 
  dplyr::select(starts_with("scale"))

jpeg(filename = "plots/weather_pop_growth/scale_weather_coef.jpeg",
     width = 7, height = 7, units = "in",res = 400)
pairs.panels(sp_res, smooth = FALSE, lm = TRUE,  ellipses = FALSE)
dev.off()

##__________________________________________________________________________________________________
#### 6. Save data ####

saveRDS(pgr_weather_res, file = "data/pgr_weather/pgr_weather_res.RDS")





