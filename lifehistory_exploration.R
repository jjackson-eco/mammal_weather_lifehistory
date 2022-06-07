#####################################################
##                                                 ##
##     Global climate and population dynamics      ##
##                                                 ##
##         Life-history data exploration           ##
##                                                 ##
##               October 26th 2020                 ##
##                                                 ##
#####################################################

## Life-history data: Demographic data from DSKI (ALHD, PanTHERIA and AnAge) 
## and body size from ALHD

## Exploration and Data summaries for analysis

rm(list = ls())
options(width = 100)

library(tidyverse)
library(patchwork)
library(psych)

##__________________________________________________________________________________________________
#### 1. Loading data ####

load("../rawdata/mam_dski.RData")
glimpse(mam_dski)

load("../rawdata/mam_alhd.RData")
glimpse(mam_alhd)

load("../rawdata/mam_comadre.RData")
glimpse(mam_comadre)

##__________________________________________________________________________________________________
#### 2. Summary of the DSKI data ####

# Species-level data from ALHD, PanTHERIA and AnAge

table(mam_dski$varname)

# variable names
mam_dski %>% 
  filter(demovar %in% c("Maximum recorded lifespan", "Litter/Clutch size")) %>% 
  group_by(varname) %>% summarise(n = n(), dem = demovar[1], db = sourcedb[1]) %>% 
  arrange(dem)

# unique species names to iterate over
spp <- unique(mam_dski$gbif.species.id)

# summary data per species from all 3 databases
longevity_litter <- bind_rows(lapply(spp, function(x){
  
  c_lit = filter(mam_dski, gbif.species.id == x & 
                   varname %in% c("litter_or_clutch_size_n", 
                                  "Litter.Clutch.size",
                                  "X151LitterSize"))
  c_lon = filter(mam_dski, gbif.species.id == x & 
                   varname %in% c("maximum_longevity_y",
                                "Maximum.longevity..yrs.",
                                "X171MaxLongevitym")) %>% 
    # change PanTHERIA max longevity from PanTHERIA to years from months
    mutate(value = if_else(varname == "X171MaxLongevitym", 
                           value/12, value))
  
  if(nrow(c_lit) == 0){lit = NA}
  else{lit = mean(c_lit$value)}
  
  if(nrow(c_lon) == 0){lon = NA}
  else{lon = max(c_lon$value)}
  
  return(tibble(c_lit[1,1:6], lit = lit, lon = lon))
  
}))

# Looking at large values
filter(longevity_litter, lit > 30) # Do Alaskan marmots have a litter of 1000 and Capybara a litter of 37? almost certainly not
filter(longevity_litter, lon > 100) # these seem possible - bowhead whale gets very old

longevity_litter <- filter(longevity_litter, lit < 30) %>% 
  mutate(IUCN = if_else(is.na(IUCNstatus), 
                        "Not Assessed", IUCNstatus),
         IUCN = factor(IUCN, levels = c("Not Assessed", "DD", "LC", "NT", 
                                        "VU", "EN", "CR", "EW", "EX")))

##__________________________________________________________________________________________________
#### 3. Plots of DSKI data ####

# basic plots of max longevity vs. litter size
longevity_litter %>% 
  filter(is.na(lit) == F & is.na(lon) == F) %>% 
  ggplot(aes(x = lon, y = lit)) +
  geom_point(size = 3.5, alpha = 0.25, colour = "#232CE0")  +
  coord_cartesian(xlim = c(0,100), ylim = c(1,10)) +
  labs(x = "Maximum longevity (years)", y = "Mean litter size") +
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank()) +
  ggsave(filename = "plots/lifehistory_raw/max_longevity_litter.jpeg",
         width = 18, height = 18, units = "cm", dpi = 500)

longevity_litter %>% 
  filter(is.na(lit) == F & is.na(lon) == F) %>% 
  group_by(order) %>% mutate(n = n()) %>% filter(n > 40) %>% 
  ggplot(aes(x = lon, y = lit)) +
  geom_point(size = 1.8, alpha = 0.6, colour = "#232CE0")  +
  labs(x = "Maximum longevity (years)", y = "Mean litter size") +
  facet_wrap(~ order, scales = "free", ncol = 5) +
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 10)) +
  ggsave(filename = "plots/lifehistory_raw/max_longevity_litter_order.jpeg",
         width = 30, height = 18, units = "cm", dpi = 500)

# longevity and litter size with IUCN
lonlit_IUCNsum <- longevity_litter %>% 
  group_by(IUCN) %>% 
  summarise(mn_lon = mean(lon, na.rm = T), md_lon = median(lon, na.rm = T),
            se_lon = sd(lon, na.rm = T)/sqrt(n()),
            mn_lit = mean(lit, na.rm = T), md_lit = median(lit, na.rm = T),
            se_lit = sd(lit, na.rm = T)/sqrt(n()))

# IUCN redlist and demography
IUCN_lon <- longevity_litter %>% 
  filter(is.na(lit) == F & is.na(lon) == F) %>% 
  ggplot(aes(x = IUCN, y = lon)) +
  geom_jitter(width = 0.1, size = 0.8, 
              alpha = 0.2, colour = "#232CE0") +
  geom_violin(alpha = 0) +
  geom_point(data = lonlit_IUCNsum, aes(y = md_lon), size = 4) +
  scale_x_discrete(labels = c("Not\nAssessed", "Data\nDeficient", 
                              "Least\nConcern", "Near\nThreatened", 
                              "Vulnerable", "Endangered", "Critically\nEndangered", 
                              "Extinct in\nthe Wild", "Extinct")) +
  labs(x = NULL, y = "Maximum longevity (years)") +
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 13))

IUCN_lit <- longevity_litter %>% 
  filter(is.na(lit) == F & is.na(lon) == F) %>% 
  ggplot(aes(x = IUCN, y = lit)) +
  geom_jitter(width = 0.1, size = 0.8, 
              alpha = 0.2, colour = "#232CE0") +
  geom_violin(alpha = 0) +
  geom_point(data = lonlit_IUCNsum, aes(y = md_lit), size = 4) +
  scale_x_discrete(labels = c("Not\nAssessed", "Data\nDeficient", 
                              "Least\nConcern", "Near\nThreatened", 
                              "Vulnerable", "Endangered", "Critically\nEndangered", 
                              "Extinct in\nthe Wild", "Extinct")) +
  scale_y_continuous(breaks = seq(0,20, by = 2)) +
  labs(x = "IUCN redlist status", y = "Litter size") +
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 13))

ggsave(IUCN_lon / IUCN_lit, 
      filename = "plots/lifehistory_raw/IUCN_lonlit.jpeg",
      width = 30, height = 30, units = "cm", dpi = 500)

##__________________________________________________________________________________________________
#### 4. Summary of ALHD body size data and plots ####

# Summary jitter plot
mam_alhd %>% 
  mutate(bodymass = as.numeric(scale(log(adult_body_mass_g)))) %>% 
  ggplot(aes(x = 1, y = bodymass)) +
  geom_violin(fill = "lightsteelblue") +
  geom_jitter(alpha = 0.3,  colour = "#232CE0", width = 0.13) +
  geom_hline(aes(yintercept = median(bodymass))) +
  scale_y_continuous(breaks = seq(-2,4, by = 1)) +
  labs(x = NULL, y = "Scaled adult body mass") +
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank()) +
  ggsave(filename = "plots/lifehistory_raw/adult_bodymass_distribution.jpeg",
         width = 15, height = 16, units = "cm", dpi = 500)

# summary data to aggregate
bodymass <- mam_alhd %>% 
  mutate(bodymass = as.numeric(scale(log(adult_body_mass_g)))) %>% 
  dplyr::select(species, gbif.species, gbif.species.id, bodymass)

##__________________________________________________________________________________________________
#### 5. Aggregate life-history data ####

lifehistory <- longevity_litter %>% 
  # scale columns for analysis
  mutate(litter = as.numeric(scale(log(lit))),
         longevity = as.numeric(scale(log(lon)))) %>% 
  dplyr::select(-c(lit, lon, IUCN)) %>% 
  filter(is.na(gbif.species.id) == F) %>% 
  full_join(x = ., y = dplyr::select(bodymass, gbif.species,
                                     gbif.species.id, bodymass),
            by = c("gbif.species", "gbif.species.id")) %>% 
  # sort out original species names for those not in long-lit (have to do separately cos of potential overlap) 
  left_join(x = ., y = dplyr::select(bodymass, species,gbif.species.id), 
            by = "gbif.species.id") %>% 
  mutate(species = if_else(is.na(species.x) == T, species.y, species.x)) %>% 
  dplyr::select(species, gbif.species, gbif.species.id, order, family, IUCNstatus, 
                litter, longevity, bodymass)

# how many complete records are there
nrow(filter(lifehistory, is.na(litter) == F & 
              is.na(longevity) == F & 
              is.na(bodymass) == F))

filter(lifehistory, is.na(bodymass) == F & 
         is.na(longevity) == T & 
         is.na(litter) == T)

##__________________________________________________________________________________________________
#### 6. Aggregated plots to show body mass patterns ####

bm_lon <- ggplot(lifehistory, aes(x = bodymass, y = longevity)) +
  geom_point(alpha = 0.3, colour = "#232CE0", size = 2) +
  labs(x = "Scaled adult bodymass", y = "Scaled maximum longevity") +
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank())

bm_lit <- ggplot(lifehistory, aes(x = bodymass, y = litter)) +
  geom_point(alpha = 0.3, colour = "#232CE0", size = 2) +
  labs(x = "Scaled adult bodymass", y = "Scaled mean litter size") +
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank())

ggsave(bm_lon + bm_lit, 
       filename = "plots/lifehistory_raw/bodymass_lonlit.jpeg",
       width = 30, height = 15, units = "cm", dpi = 500)

##__________________________________________________________________________________________________
#### 6. Save aggregated species-level data for analysis ####

save(lifehistory, file = "data/lifehistory.RData")

##__________________________________________________________________________________________________
#### 7. Combining with comadre data ####

# comadre at the species level
mam_comadre_spp <- mam_comadre %>% 
  group_by(gbif.species.id) %>% 
  summarise(gbif.species = gbif.species[1],
            generation_time = mean(generation_time, na.rm = T),
            life_expectancy = mean(life_expectancy, na.rm = T),
            adult_survival = mean(adult_survival, na.rm = T)) %>% 
  ungroup()

# combine
lifehistory_all <- lifehistory %>% 
  left_join(x = ., y = select(mam_comadre_spp, -2), by = "gbif.species.id")

# pairs plots
lifehistory_complete <- lifehistory_all %>% 
  na.omit() %>% 
  mutate_at(.vars = 10:12, .funs = function(x){log(x + 1)}) %>% 
  rename_all(.funs = str_to_sentence) %>% 
  rename(`Generation time` = Generation_time, 
         `Life expectancy` = Life_expectancy,
         `Adult survival` = Adult_survival)

jpeg("plots/lifehistory_raw/life_history_covariance.jpeg",
     width = 15, height = 15, units = "cm",res = 600)
pairs.panels(lifehistory_complete[,7:12], 
             ellipses = FALSE, lm = TRUE, method = "spearman")
dev.off()







