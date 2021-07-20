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
library(psych)

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
             size = 2, colour = "blue") +
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
  ggsave(filename = "plots/annual_abundance/annual_abundance_options/wtd_raw.jpeg",
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
pgr_1989_resid10 <- tibble(x = 1989, xend = 1990, 
                   y = wtd_residab[wtd_residab$year == 1989,]$resid_ab_10,
                   yend = wtd_residab[wtd_residab$year == 1990,]$resid_ab_10) %>% 
  mutate(pgr = yend/y)
pgr_1989_resid10

resid_ab_recenter <- ggplot(wtd_residab, aes(x = year, y = resid_ab_10)) +
  geom_hline(yintercept = 10, linetype = "dashed") +
  geom_point(size = 3) +
  geom_segment(data = pgr_1989_resid10, aes(x = x, y = y, xend = xend, yend = yend)) +
  annotate("text", x = 1990, y = 9.7,
           label = expression(paste(r[1989], " = ", frac(C[1990],C[1989]), 
                                    " = ", frac(9.76, 10.2), " = ",0.96))) +
  labs(x = "Year (t)", y = "Centered residual abundance (C)") +
  theme_bw(base_size = 12)

ggsave(grid.arrange(resid_ab, resid_ab_recenter, ncol = 1),
       filename = "plots/annual_abundance/annual_abundance_options/option1_center_resid_ab.jpeg",
       width = 7, height = 8, units = "in", dpi = 400)

##__________________________________________________________________________________________________
#### 2. Residual abundance and population change ####

# Using population change instead of per-capita growth rate
pgr_1989_resid <- tibble(x = 1989, xend = 1990, 
                           y = wtd_residab[wtd_residab$year == 1989,]$resid_ab,
                           yend = wtd_residab[wtd_residab$year == 1990,]$resid_ab) %>% 
  mutate(pgr = yend - y)
pgr_1989_resid

resid_ab_popchange <- ggplot(wtd_residab, aes(x = year, y = resid_ab)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 3) +
  geom_segment(data = pgr_1989_resid, aes(x = x, y = y, xend = xend, yend = yend)) +
  annotate("text", x = 1990, y = -0.31,
           label = expression(paste(delta[1989], " = ", R[1990] - R[1989], 
                                    " = ", -0.24 -0.19, " = ", -0.43))) +
  labs(x = "Year (t)", y = "Residual abundance (R)") +
  theme_bw(base_size = 12)

ggsave(grid.arrange(resid_ab, resid_ab_popchange, ncol = 1),
       filename = "plots/annual_abundance/annual_abundance_options/option2_resid_ab_popchange.jpeg",
       width = 7, height = 8, units = "in", dpi = 400)

##__________________________________________________________________________________________________
#### 3. ln Abundance and per-capita population growth rates ####

pgr_1989 <- tibble(x = 1989, xend = 1990, 
                   y = wtd_residab[wtd_residab$year == 1989,]$ln_abundance,
                   yend = wtd_residab[wtd_residab$year == 1990,]$ln_abundance) %>% 
  mutate(pgr = yend/y)
pgr_1989

# ln abundance plot
ggplot(wtd_residab, aes(x = year, y = ln_abundance)) +
  geom_segment(data = pgr_1989, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_point(size = 3) + 
  annotate("text", x = 1990, y = 11.25,
           label = expression(paste(r[1989], " = ",frac(X[1990],X[1989]), 
                                    " = ", frac(11.6, 12.1), " = ", 0.96))) +
  labs(x = "Year (t)", y = "ln Number of Deer (X)") +
  theme_bw(base_size = 12) +
  ggsave(filename = "plots/annual_abundance/annual_abundance_options/option3_ln_ab.jpeg",
         width = 7, height = 5, units = "in", dpi = 400)

##__________________________________________________________________________________________________
#### 4. R population growth rates - ln of raw abundance changes ####

pgr_1989 <- tibble(x = 1989, xend = 1990, 
                   y = wtd_residab[wtd_residab$year == 1989,]$raw_abundance,
                   yend = wtd_residab[wtd_residab$year == 1990,]$raw_abundance) %>% 
  mutate(pgr = log(yend/y))
pgr_1989

# ln abundance plot
ggplot(wtd_residab, aes(x = year, y = raw_abundance)) +
  geom_segment(data = pgr_1989, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_point(size = 3) + 
  annotate("text", x = 1997, y = 150000,
           label = expression(paste(R[1989], " = ",ln, frac(N[1990],N[1989]), 
                                    " = ", ln, frac(110626, 177356), " = ", -0.47))) +
  labs(x = "Year (t)", y = "Number of Deer (N)") +
  theme_bw(base_size = 12) +
  ggsave(filename = "plots/annual_abundance/annual_abundance_options/option4_R.jpeg",
         width = 7, height = 5, units = "in", dpi = 400)


##__________________________________________________________________________________________________
#### 5. Multivariate time-series analysis ####

# residual abundance timeseries style
ggplot(wtd_residab, aes(x = year, y = resid_ab)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(size = 1) +
  labs(x = "Year (t)", y = "Stationary abundance (y)") +
  theme_bw(base_size = 12) +
  ggsave("plots/annual_abundance/annual_abundance_options/option5_timeseries.jpeg",
         width = 7, height = 5, units = "in", dpi = 400)

##__________________________________________________________________________________________________
#### 6. Calculating 1-4 for all records ####

mam_pgr <- mam_IDblocks %>% 
  mutate(ln_abundance = log(raw_abundance + 1),
         ln_abundance = if_else(ln_abundance == 0, 0.1, ln_abundance)) %>% 
  group_by(ID_block) %>% 
  group_modify(~{
    
    # model for residuals
    mod = lm(ln_abundance ~ year, data = .)
    
    # method 1 - centered residuals
    resid_10 = residuals(mod) + 10
    t1_res10 = resid_10[-1]
    t0_res10 = resid_10[-(nrow(.))]
    
    pgr_res10 = t1_res10/t0_res10
    
    # method 2 - pop difference
    resid_ab = residuals(mod)
    t1_res = resid_ab[-1]
    t0_res = resid_ab[-(nrow(.))]
    
    pgr_res = t1_res - t0_res
    
    # method 3 - ln abundance
    t1 = .$ln_abundance[-1]
    t0 = .$ln_abundance[-(nrow(.))]
    
    pgr = t1/t0
    
    # method 4 - R
    t1R = .$raw_abundance[-1]
    t0R = .$raw_abundance[-(nrow(.))]
    
    pgr_R = log(t1/t0)
    
    # years with pop growth rates
    year = .$year[-(nrow(.))] 
    
    mutate(., option1 = c(pgr_res10,NA),
           option2 = c(pgr_res, NA),
           option3 = c(pgr, NA),
           option4 = c(pgr_R, NA))
  }) %>% 
  ungroup() %>% 
  filter(is.na(option3) == F & option3 < 2) 

glimpse(mam_pgr)

## Plotting
jpeg(filename = "plots/annual_abundance/annual_abundance_options/pop_growth_rate_correlations.jpeg",
     width = 10, height = 10, units = "in",res = 400)
pairs.panels(dplyr::select(mam_pgr, contains("option")),
             smooth = FALSE, lm = TRUE,  ellipses = FALSE)
dev.off()

# looking at pgrs from option 3
hist(mam_pgr$option3)

ggplot(mam_pgr, aes(x = option3, y = option4)) +
  geom_point()








