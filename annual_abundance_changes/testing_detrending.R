#####################################################
##                                                 ##
##     Global climate and population dynamics      ##
##                                                 ##
##    Testing linear detrending with simulation    ##
##                                                 ##
##                 May 8th 2020                    ##
##                                                 ##
#####################################################
rm(list = ls())
options(width = 100)

library(tidyverse)
library(grid)
library(gridExtra)

##__________________________________________________________________________________________________
#### 1. Set up data and model ####

# Testing a trend of 1.2
set.seed(100)
dat <- tibble(x = -10:10, y = 1.2*x - 4) %>% 
  mutate(y = y + rnorm(20, sd = 10))

mod <- lm(y ~ x, data = dat)

mod$coefficients[2]

##__________________________________________________________________________________________________
#### 2. Finding orthogonal (perpedicular) line and projections ####

# Orthogonal slope and y intercepts
orthSlope <- (-1) / mod$coefficients[2]
orthIntercepts <- dat$y - (orthSlope*dat$x)

# Projected x and y values based on the orthoganal
projX <- (mod$coefficients[1] - orthIntercepts) / # diff between line y intercept and orthogonal y intercepts
  (orthSlope - mod$coefficients[2]) # difference between the two slopes
projY <- orthIntercepts + (orthSlope * projX)

# Add in vertical fitted points from the mod too i.e. where each point is on the line
dat_proj <- mutate(dat, projX = projX, projY = projY,
                   projY_vert = fitted(mod)) 

##__________________________________________________________________________________________________
#### 3. Plots of the different residuals ####

# first plot to look at the orthogonal intercepts and orthogonal slope
ggplot(dat, aes(x,y)) + 
  geom_smooth(method = "lm", se = F,
              colour = "black") +
  geom_vline(xintercept = 0) +
  geom_abline(intercept = 0,
              slope = orthSlope) +
  geom_segment(aes(x = 0, xend = x,
                   y = orthIntercepts, yend = y),
               linetype = "dashed", colour = "red") +
  geom_point() +
  coord_cartesian(xlim = c(-10,10), ylim = c(-10,10)) 

# basic plot
basic <- ggplot(dat_proj, aes(x, y)) +
  geom_smooth(method = "lm", se = F, colour = "black") +
  geom_point(size = 2) +
  coord_cartesian(xlim = c(-15,15), ylim = c(-15,15)) +
  theme_bw(base_size = 12)

# standard residuals
vert_resid <- ggplot(dat_proj, aes(x, y)) +
  geom_segment(aes(xend = x, yend = projY_vert),
               linetype = "dashed", colour = "darkgreen") +
  geom_smooth(method = "lm", se = F, colour = "black") +
  geom_point(aes(x = x, y = projY_vert), 
             shape = 1, colour = "darkgreen") +
  geom_point(size = 2) +
  coord_cartesian(xlim = c(-15,15), ylim = c(-15,15)) +
  theme_bw(base_size = 12)

# orthogonal residuals
orth_resid <- ggplot(dat_proj, aes(x, y)) +
  geom_segment(aes(xend = projX, yend = projY),
               linetype = "dashed", colour = "blue") +
  geom_smooth(method = "lm", se = F, colour = "black") +
  geom_point(aes(x = projX, y = projY), 
             shape = 1, colour = "blue") +
  geom_point(size = 2) +
  coord_cartesian(xlim = c(-15,15), ylim = c(-15,15)) +
  theme_bw(base_size = 12)

ggsave(grid.arrange(basic, vert_resid, orth_resid, ncol = 3),
       filename = "plots/annual_abundance/detrending_test.jpeg",
       width = 12, height = 4, units = "in", dpi = 400)


##__________________________________________________________________________________________________
#### 4. Comparing orthogonal and vertical residuals ####

dat_resid <- dat_proj %>% 
  mutate(orth_resid = sqrt((projX - x)^2 + (projY - y)^2), # equation for the cartesian (euclidian) distance between two points
         vert_resid = resid(mod))

ggplot(dat_resid, aes(x = vert_resid, y = orth_resid)) +
  geom_point()

# These cartesian distances can't be negative
# Also probably don't need this orthogonal approach, as there isn't error in x i.e. The independent variable
# of years doesn't have error around it, 
# Therefore, it's probably sensible to use vertical residuals. They are also linked tightly, particularly in
# examples with few data points like this.




