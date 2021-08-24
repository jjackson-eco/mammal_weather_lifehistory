##############################################################
##                                                          ##
##          Global climate and population dynamics          ##
##                                                          ##
## Spatially autocorrelated meta-regression for Temperature ##
##                                                          ##
##                     Aug 23rd 2021                        ##
##                                                          ##
##############################################################

## Using a nearest neighbours spatial weighting to fit a spatially autocorrelated
## meta-regression model on the GAM weather coefficients.

rm(list = ls())
options(width = 100)

# general
library(tidyverse)
library(viridis)
library(cowplot)
library(flextable)

# spatial
library(sf)
library(spdep)

# bayesian
library(brms)
library(tidybayes)
library(ggdist)
library(ggridges)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 1. Load data ####

load("data/mammal_analysis_data_GAM.RData", verbose = TRUE)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 2. Nearest neighbours spatial weights matrix for analysis ####

# specify the coordinates of the data
coordinates(mam_coef) <- ~ Longitude + Latitude

# return k nearest neighbours for each coordinate point
knea <- knearneigh(coordinates(mam_coef), longlat = TRUE)

# convert to a neighbours list
neighbours <- knn2nb(knea)

# Spatial weighting matrix
Wmat <- nb2mat(neighbours, style = "W")

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3. Temperature model ####

## Base model
set.seed(666)
temp_base <- brm(
  coef_temp ~ 1 + sample_size + biome + (1| species),  
  data = mam_coef, family = gaussian(),
  prior = c( # lagsar gets a flat prior bound between 0 and 1
    prior(normal(0, 0.3), class =  Intercept),
    prior(normal(0, 0.3), class = b, coef = "sample_size"),
    prior(exponential(8), class = sd, group = "species")),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 4000, warmup = 2000
)

## spatial model
set.seed(666)
temp_sp <- brm(
  coef_temp ~ 1 + sample_size + biome + sar(Wmat, type = "lag") + (1| species),  
  data = mam_coef, family = gaussian(),
  data2 = list(Wmat = Wmat),
  prior = c( # lagsar gets a flat prior bound between 0 and 1
    prior(normal(0, 0.3), class =  Intercept),
    prior(normal(0, 0.3), class = b, coef = "sample_size"),
    prior(exponential(8), class = sd, group = "species")),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 4000, warmup = 2000
)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 4. Model comparisons ####

## Model comparisons
temp_base <- add_criterion(temp_base, criterion = c("loo","waic"))
temp_sp <- add_criterion(temp_sp, criterion = c("loo","waic"))

mod_comp_temp <- as.data.frame(loo_compare(temp_base, temp_sp, criterion = "loo"))

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 5. Posterior plots for spatial ####

temp_colour <- "#990a80"

sp_temp_plot <- temp_sp %>%
  gather_draws(`lagsar|sd.*|b_Intercept|b_sample_size|sigma`, regex = TRUE) %>% #tidybayes
  ungroup() %>% 
  mutate(spatial = if_else(.variable == "lagsar", "yes", "no")) %>% 
  ggplot(aes(y = .variable, x = .value, fill = spatial)) + 
  stat_halfeye(show.legend = FALSE) +
  geom_vline(xintercept = 0, size = 0.8) +
  scale_fill_manual(values = c("grey", temp_colour)) +
  scale_x_continuous(breaks = seq(-0.75,0.75,by = 0.25)) +
  scale_y_discrete(labels = c(expression(paste("Global intercept ", bar(alpha))),
                              expression(paste("Sample size ", beta[N])),
                              "Spatial autocorrelation (SAR)",
                              expression(paste("Species level variance ", sigma[SPECIES])),
                              "Population-level variance")) +
  labs(x = "Posterior estimate", y = NULL) +
  theme_ridges(center_axis_labels = TRUE, grid = T, line_size = 0.3) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))


ggsave(plot = sp_temp_plot, 
       filename =  "plots/meta_regression/spatial_autocorrelation_posterior.jpeg",
       device = "jpeg", width = 15, height = 14, units = "cm", dpi = 500)

## model comparisons plot
mod_comp_temp %>% 
  mutate(model = c("Base model", "Spatial autocorrelation (SAR)")) %>% 
  dplyr::select(model, elpd_loo, se_elpd_loo, elpd_diff, se_diff, looic) %>% 
  flextable(cwidth = 1.2) %>% 
  set_header_labels(model = "Model",
                    elpd_loo = "LOO elpd",
                    se_elpd_loo = "LOO elpd error",
                    elpd_diff = "elpd difference",
                    se_diff = "elpd error difference",
                    looic = "LOO information criterion") %>% 
  colformat_num(digits = 2) %>% 
  save_as_image("plots/meta_regression/spatial_autocorrelation_model_comparison.png")
  









