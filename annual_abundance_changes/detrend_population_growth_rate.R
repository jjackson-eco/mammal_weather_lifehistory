######################################################
##                                                  ##
##     Global climate and population dynamics       ##
##                                                  ##
##   Linear Detrending + Per capita growth rates    ##
##            of abundance timeseries               ##
##                                                  ##
##                May 15th 2020                     ##
##                                                  ##
######################################################
rm(list = ls())
options(width = 100)

library(tidyverse)
library(viridis)
library(gridExtra)
library(grid)

##__________________________________________________________________________________________________
#### 1. Load and set up data ####

# Loading the raw mammal data with the blocks
load("../rawdata/mam_IDblocks.RData")
glimpse(mam_IDblocks)

# Checking number of blocks and records
mam_IDblocks %>% 
  group_by(ID) %>% 
  summarise(n_block = n_distinct(block), n_ID = 1) %>% 
  ungroup() %>% 
  mutate(num = n_ID*n_block) %>% 
  summarise(sum(num))
  
n_distinct(mam_IDblocks$ID_block) # Seem to match up ok

##__________________________________________________________________________________________________
#### 2. Linear detrend for each block of each study ####

# extracting the residuals after fitting a linear trend to each timeseries
mam_detrend <- mam_IDblocks %>% 
  group_by(ID_block) %>% 
  group_modify(~{
    mod = lm(ln_abundance ~ year, data = .) 
    
    resid_ab = mod$residuals + 10 # Centre around 10 for sensible population growth rate calculations
    
    mutate(., residual_abundance = resid_ab,
           coef = mod$coefficients[2])
  }) %>% 
  ungroup()

# Centering around 10 changes the values <- HOW DO WE SORT THIS OUT??

#test
# set.seed(10)
# dat <- tibble(year = 1991:2020, residab = rnorm(30)) %>% 
#   mutate(residab2 = residab + 100) %>% 
#   mutate(pgr  = c(.$residab[-1]/.$residab[-(length(.$residab))], NA),
#          pgr2 = c(.$residab2[-1]/.$residab2[-(length(.$residab2))], NA))

##__________________________________________________________________________________________________
#### 3. Rough look at coefficients of linear trend across taxa ####

mam_detrend %>% 
  group_by(Family) %>% 
  summarise(mn_coef = mean(coef),
            se_coef = sd(coef)/sqrt(n())) %>% 
  ungroup() %>% 
  arrange(desc(mn_coef)) %>% 
  mutate(Family = factor(Family, levels = .$Family)) %>% 
  ggplot(aes(x = Family, y = mn_coef)) +
  geom_hline(yintercept = 0, size = 0.5, linetype = "dashed") +
  geom_errorbar(aes(ymax = mn_coef + se_coef,
                    ymin = mn_coef - se_coef),
                width = 0.01) +
  geom_point() +
  labs(y = "Mean coefficient of linear abundance trend") +
  coord_flip() +
  theme_bw(base_size = 17) +
  theme(axis.text.y = element_text(size = 12)) +
  ggsave("plots/annual_abundance/linear_trend_and_detrending/linear_abundance_coefficients.jpeg",
         width = 9, height = 13, units = "in", dpi = 400)

##__________________________________________________________________________________________________
#### 4. Calculate per-capita residual growth rate between time t and t+1 ####

# This removes one year from each ID_block - Minimum of 4 years - Does this matter?
mammal <- mam_detrend %>% 
  group_by(ID_block) %>% 
  group_modify(~{
    resid_t0 = .$residual_abundance[-(length(.$residual_abundance))]
    resid_t1 = .$residual_abundance[-1]
    
    mutate(., pop_growth_rate = c(resid_t1/resid_t0,NA))
  }) %>% 
  ungroup() %>% 
  filter(is.na(pop_growth_rate) == F)

##__________________________________________________________________________________________________
#### 5. Looking at population growth rates ####

# 5a. Abundance and Population growth rates
# Look very normally distributed - linear relationship seems justified
# residual abundance
resid_ab <- ggplot(mammal, aes(x = residual_abundance)) +
  geom_histogram(bins = 30, colour = "black", size = 0.1,
                 fill = "lightblue") +
  labs(x = "Scaled residual abundance", y = "Frequency") +
  facet_wrap(~Order, scales = "free_y", ncol = 5) +
  theme_bw(base_size = 22)
  
# population growth rates
pop_growth <- ggplot(mammal, aes(x = pop_growth_rate)) +
  geom_histogram(bins = 30, colour = "black", size = 0.1,
                 fill = "lightblue") +
  labs(x = "Population growth rate", y = "Frequency") +
  facet_wrap(~Order, scales = "free_y", ncol = 5) +
  theme_bw(base_size = 22)

ggsave(grid.arrange(resid_ab, pop_growth, ncol = 1),
       filename = "plots/annual_abundance/linear_trend_and_detrending/abundance_histograms.jpeg",
       width = 17, height = 20, units = "in", dpi = 400)

# 5b. Evidence for density dependence
ggplot(mammal, aes(x = residual_abundance, y = pop_growth_rate)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_point(alpha = 0.15, colour = viridis(10)[6], size = 2.5) +
  geom_smooth(method = "lm", se = F, colour = "black") +
  labs(x = "Scaled residual abundance", y = "Population growth rate") +
  theme_bw(base_size = 18) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggsave(filename = "plots/annual_abundance/linear_trend_and_detrending/density_dependence_lpd.jpeg",
         width = 7, height = 7, units = "in", dpi = 400)

# 5c. Study length update
mammal %>% 
  group_by(ID) %>% 
  summarise(record_length = n()) %>% 
  ggplot(aes(x = record_length)) +
  geom_histogram(bins = 15, fill = "lightblue", colour = "black") +
  labs(x = "Length of abundance record (years)", y = "Frequency") +
  theme_bw(base_size = 15) +
  ggsave("plots/annual_abundance/record_length.jpeg",
         width = 6, height = 5, units = "in", dpi = 400)

# ##__________________________________________________________________________________________________
# #### 6. Saving data ####
# 
# mammal <- dplyr::select(mammal, -c(raw_abundance,scaled_abundance))
# 
# save(mammal, file = "data/mammal.RData")








