######################################################
##                                                  ##
##     Global climate and population dynamics       ##
##                                                  ##
## Linear Detrending of Annual Abundance Timeseries ##
##                                                  ##
##                May 12th 2020                     ##
##                                                  ##
######################################################
rm(list = ls())
options(width = 100)

library(tidyverse)
library(gridExtra)
library(grid)

##__________________________________________________________________________________________________
#### 1. Load and set up data ####

# Loading the raw mammal data with the blocks
load("../rawdata/mammal_data.RData")
glimpse(mammal)

# Checking number of blocks and records
mammal %>% 
  group_by(ID) %>% 
  summarise(n_block = n_distinct(block), n_ID = 1) %>% 
  ungroup() %>% 
  mutate(num = n_ID*n_block) %>% 
  summarise(sum(num))
  
n_distinct(mammal$ID_block) # Seem to match up ok

##__________________________________________________________________________________________________
#### 2. Linear detrend for each block of each study ####

# extracting the residuals after fitting a linear trend to each timeseries
mam_detrend <- mammal %>% 
  group_by(ID_block) %>% 
  group_modify(~{
    mod = lm(scaled_abundance ~ year, data = .) 
    
    mutate(., residual_abundance = mod$residuals,
           coef = mod$coefficients[2])
  }) %>% 
  ungroup()

# Look very normally distributed - linear relationship seems justified
ggplot(mam_detrend, aes(x = residual_abundance)) +
  geom_histogram(bins = 30, colour = "black", size = 0.1,
                 fill = "lightblue") +
  labs(x = "Residual scaled abundance", y = "Frequency") +
  facet_wrap(~Order, scales = "free_y", ncol = 5) +
  theme_bw(base_size = 22) +
  ggsave("plots/annual_abundance/Residaul_abundance_histogram.jpeg",
         width = 14, height = 12, units = "in", dpi = 400)

##__________________________________________________________________________________________________
#### 3. Rough look at coefficients across taxa ####

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
  ggsave("plots/annual_abundance/linear_abundance_coefficients.jpeg",
         width = 9, height = 13, units = "in", dpi = 400)

##__________________________________________________________________________________________________
#### 4. Saving data ####

mam_detrend <- dplyr::select(mam_detrend, -c(raw_abundance,scaled_abundance))

save(mam_detrend, file = "data/mam_detrend.RData")








