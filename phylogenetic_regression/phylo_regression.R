####################################################
##                                                ##
##     Global climate and population dynamics     ##
##                                                ##
##         brms Phylogenetic regression           ##
##                                                ##
##                Oct 16th 2020                   ##
##                                                ##
####################################################

## Species level variation and phylogenetic signal for mean temp/precip anomalies at 5km resolution

rm(list = ls())
options(width = 100)

library(tidyverse)
library(ape)
library(brms)
library(nlme)
library(patchwork)

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
    prior(exponential(15), class = sd, group = "spp"),
    prior(exponential(3), class = sd, group = "species"),
    prior(exponential(3), class = sigma)),
  control = list(adapt_delta = 0.95),
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
#### 6. Full model tests ####

# Simple phylogenetic model
set.seed(666)
temp_mod <- brm(
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


