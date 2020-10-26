#####################################################
##                                                 ##
##     Global climate and population dynamics      ##
##                                                 ##
##    Demographic species knowledge index data     ##
##                                                 ##
##               October 26th 2020                 ##
##                                                 ##
#####################################################
rm(list = ls())
options(width = 100)

library(tidyverse)
library(patchwork)

##__________________________________________________________________________________________________
#### 1. Loading data ####

load("../rawdata/mam_dski.RData")
glimpse(mam_dski)

##__________________________________________________________________________________________________
#### 2. Summary of the ALHDB data ####

table(mam_dski$varname)

spp <- unique(mam_dski$species)

longevity_litter <- bind_rows(lapply(spp, function(x){
  
  c_lit = filter(mam_dski, species == x & varname == "litter_or_clutch_size_n")
  c_lon = filter(mam_dski, species == x & varname == "maximum_longevity_y")
  
  if(nrow(c_lit) == 0){lit = NA}
  else{lit = mean(c_lit$value)}
  
  if(nrow(c_lon) == 0){lon = NA}
  else{lon = max(c_lon$value)}
  
  return(tibble(c_lit[1,1:6], lit = lit, lon = lon))
  
}))


# Looking at large values
filter(longevity_litter, lit > 50) # Do Alaskan marmots have a littler of 1000? almost certainly not
filter(longevity_litter, lon > 100) # these seem possible - bowhead whale gets very old

longevity_litter <- filter(longevity_litter, lit < 50) %>% 
  mutate(IUCN = if_else(is.na(IUCNstatus), 
                        "Not Assessed", IUCNstatus),
         IUCN = factor(IUCN, levels = c("Not Assessed", "DD", "LC", "NT", 
                                        "VU", "EN", "CR", "EW", "EX")))

##__________________________________________________________________________________________________
#### 2. Plots of ALHDB data ####

# basic plots of max longevity vs. litter size
longevity_litter %>% 
  filter(is.na(lit) == F & is.na(lon) == F) %>% 
  ggplot(aes(x = lon, y = lit)) +
  geom_point(size = 3.5, alpha = 0.25, colour = "#232CE0")  +
  coord_cartesian(xlim = c(0,100), ylim = c(1,10)) +
  labs(x = "Maximum longevity (years)", y = "Mean litter size") +
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank()) +
  ggsave(filename = "plots/dski_raw/max_longevity_litter.jpeg",
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
  ggsave(filename = "plots/dski_raw/max_longevity_litter_order.jpeg",
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
      filename = "plots/dski_raw/IUCN_lonlit.jpeg",
      width = 30, height = 30, units = "cm", dpi = 500)




