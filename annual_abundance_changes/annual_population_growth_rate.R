#####################################################
##                                                 ##
##     Global climate and population dynamics      ##
##                                                 ##
##         Annual population growth rates          ##
##                                                 ##
##                 June 8th 2020                   ##
##                                                 ##
#####################################################

## Calculating and investigating annual abundance changes i.e. Raw population growth rates

rm(list = ls())
options(width = 100)

library(tidyverse)
library(viridis)
library(grid)
library(gridExtra)
library(psych)

##__________________________________________________________________________________________________
#### 1. Load data ####

load("../rawdata/mam_IDblocks.RData")
glimpse(mam_IDblocks)

load("../rawdata/mam.RData", verbose = T)

##__________________________________________________________________________________________________
#### 2. Calculate per-capita residual growth rate between time t and t+1 ####

# This will remove one year from each ID_block - Does this matter?
mammal <- mam_IDblocks %>% 
  group_by(ID_block) %>% 
  group_modify(~{
    t0 <- .$ln_abundance[-(length(.$ln_abundance))]
    t1 <- .$ln_abundance[-1]
    
    mutate(., pop_growth_rate = c(t1/t0,NA))
  }) %>% 
  ungroup() %>% 
  filter(is.na(pop_growth_rate) == F)

##__________________________________________________________________________________________________
#### 3. Exploratory plots ####

# 3a. Histogram
ggplot(mammal, aes(x = pop_growth_rate)) +
  geom_histogram(bins = 100, fill = "black") +
  labs(x = "Per-capita population growth rate", 
       y = "Frequency") +
  theme_bw(base_size = 12) +
  ggsave(filename = "plots/annual_abundance/pop_growth_rate_histogram.jpeg",
         width = 6, height = 4, units = "in", dpi = 400)

# 3b. Density dependence
ggplot(mammal, aes(x = ln_abundance, y = pop_growth_rate)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_point(alpha = 0.15, colour = viridis(10)[6], size = 2.5) +
  geom_smooth(method = "lm", se = F, colour = "black") +
  coord_cartesian(ylim = c(0,3)) +
  labs(x = "ln Abundance", y = "Per-capita population growth rate") +
  theme_bw(base_size = 18) + 
  # theme(panel.grid.major = element_blank(),
  #       panel.grid.minor = element_blank()) +
  ggsave(filename = "plots/annual_abundance/density_dependence_mam.jpeg",
         width = 7, height = 7, units = "in", dpi = 400)

# 3c. Abundance change distribution
mamquant <- tibble(quant = quantile(mammal$pop_growth_rate, probs = c(0.01,0.99)))

ggplot(mammal, aes(x = 1, y = pop_growth_rate)) +
  geom_hline(yintercept = 1) +
  geom_violin() +
  geom_jitter(alpha = 0.2, aes(colour = pop_growth_rate), show.legend = F) +
  geom_hline(data = mamquant, aes(yintercept = quant), linetype = "dashed") +
  scale_y_continuous(breaks = seq(0,8,by = 0.5)) +
  scale_color_viridis_c() +
  labs(x = NULL, y = "Annual abundance change") +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  ggsave("plots/annual_abundance/annual_abundance_change_distribution.jpeg",
         width = 12, height = 15, units = "cm",dpi = 500)

##__________________________________________________________________________________________________
#### 4. Unusual Abundance changes ####

# Pull out the ID blocks below-above the 1%-99% quantiles 
unusual_pop_growth <- mammal %>% 
  filter(pop_growth_rate > mamquant$quant[2] | 
           pop_growth_rate < mamquant$quant[1]) 

# and the raw abundance changes
unusual_dat <- mam_IDblocks %>% 
  filter(ID_block %in% unusual_pop_growth$ID_block == TRUE)

# Plotting for all populations
ggplot(unusual_dat, aes(x = year, y = ln_abundance, group = ID_block)) +
  geom_line() +
  geom_point(data = unusual_pop_growth, aes(colour = pop_growth_rate)) +
  scale_color_gradient(low = "#e66101", high = "#5e3c99", name = "Abundance\nchange",
                       guide = guide_colorbar(barwidth = 2, barheight = 10)) +
  labs(x = "Year", y = "ln Abundance") +
  facet_wrap(~ ID_block,nrow = 10) +
  theme_bw(base_size = 15) +
  theme(strip.background = element_blank(),
        panel.grid = element_blank()) +
  ggsave("plots/annual_abundance/unusual_abundance_changes.pdf",
         width = 300, height = 297, units = "mm")

# Any time-series with a raw 0 in it
raw0s <- mam_IDblocks %>% 
  group_by(ID_block) %>% 
  mutate(raw0 = any(raw_abundance == 0)) %>%
  ungroup() %>% 
  filter(raw0 == TRUE)

only0s <- filter(raw0s, raw_abundance == 0)

ggplot(raw0s, aes(x = year, y = raw_abundance)) +
  geom_line() +
  geom_point(size = 0.6) +
  geom_point(data = only0s,colour = "firebrick", size = 1) +
  facet_wrap(~ ID_block, scales = "free") +
  labs(x = "Year", y = "Abundance") +
  theme_bw(base_size = 12) +
  theme(strip.background = element_blank(),
        panel.grid = element_blank()) +
  ggsave("plots/annual_abundance/raw0_abundance.pdf",
         width = 450, height = 210, units = "mm")

raw0_ID <- unique(raw0s$ID)

meta_raw0 <- filter(mam_meta, ID %in% raw0_ID) %>% 
  dplyr::select(ID, Binomial, Common_name, Reference, Units, Method)
View(meta_raw0)

##__________________________________________________________________________________________________
#### 5. Save ####

save(mammal, file = "../rawdata/mammal.RData")

