####################################################
##                                                ##
##     Global climate and population dynamics     ##
##                                                ##
##         brms Phylogenetic regression           ##
##                                                ##
##                Oct 16th 2020                   ##
##                                                ##
####################################################

## Species level variation, phylogenetic signal and spatial effects,
## for mean temp/precip anomalies at 5km resolution

rm(list = ls())
options(width = 100)

library(tidyverse)
library(ape)
library(brms)
library(nlme)
library(patchwork)
library(ggridges)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 1. Load data ####

## Weather effects
load("data/pgr_weather/mnanom_5km.RData", verbose = TRUE)

## Phylogenetic
load("../Phlogenetics/mam_lpd_pruned.RData", verbose = TRUE)

## Mammal abundance data
load("../rawdata/mam_IDblocks.RData")

## Weather anomaly data 
mam_chelsa_annual <- readRDS("data/mam_chelsa_annual.RDS") %>% 
  filter(scale == "scale_5km") %>% 
  dplyr::select(ID,year, weather_scale = scale, mean_temp_anomaly, mean_precip_anomaly)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 2. Mammal coefficient data ####

mam_temp <- mnanom_5km %>% 
  ungroup() %>% 
  dplyr::select(ID, ID_block, Order, biome, Latitude, coef_temp, coef_precip) %>% 
  left_join(x = ., 
            y = dplyr::select(LPD_tree_update, ID, spp = checked_speciesname),
            by = "ID") %>% 
  mutate(spp = gsub(" ", "_", spp),
         species = spp,
         ## z transform the coefficients
         coef_temp = as.numeric(scale(coef_temp)),
         coef_precip = as.numeric(scale(coef_precip)),
         #z transform absolute latitude for models
         lat = as.numeric(scale(abs(Latitude)))) 

## restrict for precipitation - dropping NAs
mam_precip <- drop_na(mam_temp, coef_precip)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3. Phylogenetic covariance matrix ####

# Trim tree to right data
mamMCC_temp <- keep.tip(mamMCC_pruned, mam_temp$spp)
mamMCC_precip <- keep.tip(mamMCC_pruned, mam_precip$spp)

# Covariance matrix - Brownian motion model
A_temp <- ape::vcv.phylo(mamMCC_temp)
A_precip <- ape::vcv.phylo(mamMCC_precip)

# Phylogenetic distance matrix
Dmat <- cophenetic(mamMCC_temp)

# implied covariance - quite unrealistic?
plot(Dmat, A_temp, type = "b", xlab = "Phylogenetic distance", ylab = "Brownian covariance")

# Exponential covariance decay from OU model (simple_phylogenetic_regression.R) -
eta_sq <- 3.90
rho_sq <- 5.32

K <- Dmat
for(i in 1:nrow(K)){
  for(j in 1:ncol(K)){
    K[i,j] <- eta_sq*exp(-rho_sq*Dmat[i,j])
  }
}

library(plot.matrix)

jpeg(filename = "plots/phylogenetic_regression/distance_matrices_temp.jpeg", 
     width = 20, height = 8, units = "in", res = 500)
par(mfrow = c(1,3))
plot(K[rownames(A_temp),rownames(A_temp)], main = "OU distance", border = NA)
plot(Dmat[rownames(A_temp),rownames(A_temp)], main = "Raw distance", border = NA)
plot(A_temp[rownames(A_temp),rownames(A_temp)], main = "Brownian motion", border = NA)
dev.off()

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 4. Exponential prior tests ####

plot(NULL, xlab = "Value", ylab = "Density", ylim = c(0,10), xlim = c(0,3))

for(i in 1:15){
  vec = density(rexp(1e4, rate = i))
  lines(vec, col = i)
  text(x = 0.05, y = max(unlist(vec[2])), labels = paste0("Exponential(",i, ")"), col = i)}

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 5. Single chain tests ####

# Simple phylogenetic model
set.seed(666)
temp_test <- brm(
  coef_temp ~ 1 + (1|gr(spp, cov = A_temp)) + (1| species), 
  data = mam_temp, family = gaussian(), 
  data2 = list(A_temp = A_temp),
  prior = c(
    prior(normal(0, 0.5), class =  Intercept),
    prior(exponential(1), class = sd, group = "spp"),
    prior(exponential(1), class = sd, group = "species"),
    prior(exponential(1), class = sigma)),
  control = list(adapt_delta = 0.95),
  chains = 1, iter = 4000, warmup = 500
)

plot(temp_test)

precip_test <- brm(
  coef_precip ~ 1 + (1|gr(spp, cov = A_precip)) + (1| species), 
  data = mam_precip, family = gaussian(), 
  data2 = list(A_precip = A_precip),
  prior = c(
    prior(normal(0, 0.5), class =  Intercept),
    prior(exponential(3), class = sd, group = "spp"),
    prior(exponential(3), class = sd, group = "species"),
    prior(exponential(3), class = sigma)),
  control = list(adapt_delta = 0.99),
  chains = 1, iter = 4000, warmup = 500
)

plot(precip_test)

# # Custom covariance matrix
# set.seed(666)
# temp_testOU <- brm(
#   coef_temp ~ 1 + (1|gr(spp, cov = K)) + (1| species), 
#   data = mam_temp, family = gaussian(), 
#   data2 = list(K = K),
#   prior = c(
#     prior(normal(0, 0.25), class =  Intercept),
#     prior(exponential(0.1), class = sd, group = "spp"),
#     prior(exponential(0.1), class = sd, group = "species"),
#     prior(exponential(0.1), class = sigma)),
#   control = list(adapt_delta = 0.95),
#   chains = 1, iter = 4000, warmup = 500
# )
# 
# plot(temp_testOU)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 6. Full models ####

##___________________________________________________________________________
## 6a. Temperature 

set.seed(666)

# Base model
temp_base <- brm(
  coef_temp ~ 1 + (1|gr(spp, cov = A_temp)) + (1| species), 
  data = mam_temp, family = gaussian(), 
  data2 = list(A_temp = A_temp),
  prior = c(
    prior(normal(0, 0.5), class =  Intercept),
    prior(exponential(1), class = sd, group = "spp"),
    prior(exponential(1), class = sd, group = "species"),
    prior(exponential(1), class = sigma)),
  control = list(adapt_delta = 0.97),
  chains = 4, cores = 4, iter = 2000, warmup = 500
)

# biome
temp_biome <- brm(
  coef_temp ~ 1 + biome + (1|gr(spp, cov = A_temp)) + (1| species), 
  data = mam_temp, family = gaussian(), 
  data2 = list(A_temp = A_temp),
  prior = c(
    prior(normal(0, 0.5), class =  Intercept),
    prior(normal(0, 0.5), class =  b),
    prior(exponential(1), class = sd, group = "spp"),
    prior(exponential(1), class = sd, group = "species"),
    prior(exponential(1), class = sigma)),
  control = list(adapt_delta = 0.99),
  chains = 4, cores = 4, iter = 2000, warmup = 500
)

# latitude
temp_lat <- brm(
  coef_temp ~ 1 + lat + (1|gr(spp, cov = A_temp)) + (1| species), 
  data = mam_temp, family = gaussian(), 
  data2 = list(A_temp = A_temp),
  prior = c(
    prior(normal(0, 0.5), class =  Intercept),
    prior(normal(0, 0.5), class =  b),
    prior(exponential(1), class = sd, group = "spp"),
    prior(exponential(1), class = sd, group = "species"),
    prior(exponential(1), class = sigma)),
  control = list(adapt_delta = 0.97),
  chains = 4, cores = 4, iter = 2000, warmup = 500
)

## Model comparisons
temp_base <- add_criterion(temp_base, criterion = c("loo", "waic"))
temp_biome <- add_criterion(temp_biome, criterion = c("loo", "waic"))
temp_lat <- add_criterion(temp_lat, criterion = c("loo", "waic"))

## Biome model looks promising
loo_compare(temp_base, temp_biome, temp_lat, criterion = "loo")

# Pareto k outliers
temp_pareto <- as.data.frame(loo(temp_base, pointwise = TRUE)$pointwise)

mam_temp %>% 
  mutate(pareto = temp_pareto$influence_pareto_k) %>% 
  ggplot(aes(x = coef_temp, y = pareto)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_text(aes(label = ID, colour = pareto), show.legend = F) +
  labs(x = "Linear temperature coefficient", y = "Pareto k estimate") +
  theme_bw(base_size = 14) + theme(panel.grid = element_blank()) +
  ggsave(filename = "plots/phylogenetic_regression/repeated_obs_spatial/pareto_k_temp.jpeg",
         width = 6, height = 6, units = "in", dpi = 400)

# biome 2 - no bat Myotis_emarginatus
mam2 <- filter(mam_temp, spp != "Myotis_emarginatus")
mamMCC2 <- keep.tip(mamMCC_temp, mam2$spp)
A_temp2 <- ape::vcv.phylo(mamMCC2)

temp_biome <- brm(
  coef_temp ~ 1 + biome + (1|gr(spp, cov = A_temp2)) + (1| species), 
  data = mam2, family = gaussian(), 
  data2 = list(A_temp2 = A_temp2),
  prior = c(
    prior(normal(0, 0.5), class =  Intercept),
    prior(normal(0, 0.5), class =  b),
    prior(exponential(1), class = sd, group = "spp"),
    prior(exponential(1), class = sd, group = "species"),
    prior(exponential(1), class = sigma)),
  control = list(adapt_delta = 0.99),
  chains = 4, cores = 4, iter = 2000, warmup = 500
)

biome_post <- posterior_samples(temp_biome, pars = "biome")
colnames(biome_post) <- c("Desert", "Flooded.grassland", "Mediterranean", "Montane", 
                          "Temperate.broadleaf", "Temperate.coniferous","Temperate.grassland", 
                          "Tropical.coniferous", "Tropical.broadleaf", "Tropical.grassland", 
                          "Tropical.moist", "Tundra")


## Plot the biome effects out
tibble(biome_post, sim = 1:nrow(biome_post)) %>% 
  pivot_longer(-sim, names_to = "biome", values_to = "post") %>%
  mutate(biome = gsub("[.]", " ", biome)) %>% 
  ggplot(aes(x = post, y = biome, fill = stat(x))) +
  geom_vline(xintercept = 0) +
  geom_density_ridges_gradient(scale = 1.1) +
  scale_fill_viridis_c(option = "C", guide = F) +
  labs(x = "Posterior temperature effect", y = NULL) +
  theme_ridges(center_axis_labels = TRUE, font_size = 20, grid = F) +
  ggsave(filename = "plots/phylogenetic_regression/repeated_obs_spatial/biome_effects.jpeg",
         width = 7, height = 10, units = "in", dpi = 500)

