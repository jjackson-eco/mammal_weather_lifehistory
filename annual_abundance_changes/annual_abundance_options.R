#####################################################
##                                                 ##
##     Global climate and population dynamics      ##
##                                                 ##
##       Options for analyses for abundance        ##
##                                                 ##
##                 May 8th 2020                    ##
##                                                 ##
#####################################################
rm(list = ls())
options(width = 100)

library(tidyverse)
library(grid)
library(gridExtra)

# For the maps
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(cowplot)

##__________________________________________________________________________________________________
#### Load data and select example ####

load("../rawdata/mam_IDblocks.RData")

mam_IDblocks %>% 
  group_by(ID_block) %>% 
  mutate(n_obs = n()) %>% 
  ungroup() %>% 
  filter(n_obs == max(n_obs)) %>% 
  glimpse()

#____________________________________________
### White-tailed Deer ###
# 35 years of data from the white-tailed deer in New Brunswick Canada
wtd <- filter(mam_IDblocks, ID_block == "22039_1") %>% 
  mutate(ln_abundance = log(1 + raw_abundance))

world_sf <- ne_countries(scale = "medium", returnclass = "sf")

wtd_map <- ggplot(data = world_sf) +
  geom_sf(size = 0.1, fill = "#9DBF9E") +
  geom_point(aes(x = wtd$Longitude[1], y = wtd$Latitude[1]), 
             ize = 2, colour = "blue") +
  coord_sf(xlim = c(-145, -45), ylim = c(20, 55), expand = FALSE) +
  theme_void()

# abundance timeseries plot
wtd_ts <- ggplot(wtd, aes(x = year, y = raw_abundance)) +
  geom_point(size = 3) +
  labs(x = "Year (t)", y = "Number of Deer (N)") +
  theme_bw(base_size = 12)

## PLOTTING WITH AN INSET TO SHOW STUDY SITE <<----------------------
ggdraw(wtd_ts) +
  draw_plot({wtd_map},x = 0.54, y = 0.59,
            width = 0.46, height = 0.46) +
  ggsave(filename = "plots/annual_abundance/wtd_raw.jpeg",
         width = 7, height = 5, units = "in", dpi = 400)
  
##__________________________________________________________________________________________________
#### 1. Centered residual abundance + per capita population growth rate ####

# Linear model to detrend ln abundance
linear_mod <- lm(ln_abundance ~ year, data = wtd)
projY_vert <- fitted(linear_mod)

# Adding in residual abundance + re-centered residual abundance
wtd_residab <- wtd %>% 
  mutate(resid_ab = residuals(linear_mod),
         resid_ab_10 = residuals(linear_mod) + 10)

# residual abundance raw plot
resid_ab <- ggplot(wtd_residab, aes(x = year, y = ln_abundance)) +
  geom_segment(aes(xend = year, yend = projY_vert),
               linetype = "dashed", colour = "darkgreen") +
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = F, colour = "black") +
  labs(x = "Year (t)", y = "ln Number of Deer (X)") +
  theme_bw(base_size = 12)

# re-centered residual abundance + per-capita growth rate 
pgr_1989 <- tibble(x = 1989, xend = 1990, 
                   y = wtd_residab[wtd_residab$year == 1989,]$resid_ab_10,
                   yend = wtd_residab[wtd_residab$year == 1990,]$resid_ab_10) %>% 
  mutate(pgr = yend/y)
pgr_1989

resid_ab_recenter <- ggplot(wtd_residab, aes(x = year, y = resid_ab_10)) +
  geom_hline(yintercept = 10, linetype = "dashed") +
  geom_point(size = 3) +
  geom_segment(data = pgr_1989, aes(x = x, y = y, xend = xend, yend = yend)) +
  annotate("text", x = 1990, y = 9.7,
           label = expression(paste(r[1989], " = ", frac(C[1990],C[1989]), 
                                    " = ", frac(9.76, 10.2), " = ",0.96))) +
  labs(x = "Year (t)", y = "Centered residual abundance (C)") +
  theme_bw(base_size = 12)

ggsave(grid.arrange(resid_ab, resid_ab_recenter, ncol = 1),
       filename = "plots/annual_abundance/option1_center_resid_ab.jpeg",
       width = 5, height = 10)





