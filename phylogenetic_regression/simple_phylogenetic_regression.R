####################################################
##                                                ##
##     Global climate and population dynamics     ##
##                                                ##
##    Simple brms Phylogenetic Meta-regression    ##
##                                                ##
##                Oct 14th 2020                   ##
##                                                ##
####################################################

## Simplest single observation for each species for mean temp/precip anomalies at 5km resolution
## Testing brms

rm(list = ls())
options(width = 100)

library(tidyverse)
library(ape)
library(brms)
library(rethinking)
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
  dplyr::select(ID, Order, coef_temp, coef_precip) %>% 
  left_join(x = ., 
            y = dplyr::select(LPD_tree_update, ID, checked_speciesname),
            by = "ID") %>% 
  mutate(checked_speciesname = gsub(" ", "_", checked_speciesname)) %>% 
  # Average coefficient for each species
  group_by(checked_speciesname) %>% 
  summarise(Order = Order[1],
            coef_temp = mean(coef_temp), 
            coef_precip = mean(coef_precip),
            .groups = "drop") %>% 
  dplyr::rename(spp = checked_speciesname)

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

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 4. Prior predictive simulation ####

# sample from a Gaussian for the intercept i.e. the coefficient of the temperature effect
set.seed(100)
ppd_coefs <- tibble(weak_int = rnorm(1000, mean = 0, sd = 10),
                        reg_int = rnorm(1000, mean = 0, sd = 1),
                    sim = 1:1000)


p1 <- ggplot(data.frame(x = seq(-1,1,length.out = 10), 
                        y = seq(1,max(mam_IDblocks$ln_abundance), length.out = 10)), 
                        aes(x = x,y = y)) +
  geom_blank() +
  geom_abline(data = ppd_coefs,
              aes(slope = weak_int, intercept = 0, group = sim),
              alpha = 0.2) +
  ggtitle("Weak prior ~ Normal(0,10)") +
  labs(x = "Temperature anomaly", y = "ln Abundance") +
  theme_bw(base_size = 13) + theme(panel.grid = element_blank())


p2 <- ggplot(data.frame(x = seq(-1,1,length.out = 10), 
                        y = seq(1,max(mam_IDblocks$ln_abundance), length.out = 10)), 
             aes(x = x,y = y)) +
  geom_blank() +
  geom_abline(data = ppd_coefs,
              aes(slope = reg_int, intercept = 0, group = sim),
              alpha = 0.2) +
  ggtitle("Regularising prior ~ Normal(0,1)") +
  labs(x = "Temperature anomaly", y = "ln Abundance") +
  theme_bw(base_size = 13) + theme(panel.grid = element_blank())

p1 | p2

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 5. MCMC single chain test ####

temp_test <- brm(
  coef_temp ~ 1 + (1|gr(spp, cov = A_temp)), 
  data = mam_temp, family = gaussian(), 
  data2 = list(A_temp = A_temp),
  prior = c(
    prior(normal(0, 1), class =  Intercept),
    prior(exponential(0.5), class = sd, group = "spp"),
    prior(exponential(0.5), class = sigma)),
  control = list(adapt_delta = 0.95),
  chains = 1, iter = 2000, warmup = 500
)

precip_test <- brm(
  coef_precip ~ 1 + (1|gr(spp, cov = A_precip)), 
  data = mam_precip, family = gaussian(), 
  data2 = list(A_precip = A_precip),
  prior = c(
    prior(normal(0, 1), class =  Intercept),
    prior(exponential(0.5), class = sd, group = "spp"),
    prior(exponential(0.5), class = sigma)),
  control = list(adapt_delta = 0.95),
  chains = 1, iter = 2000, warmup = 500
)

# MCMC plots
plot(temp_test)
plot(precip_test)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 6. Sampling from the posterior ####

temp_simple <- brm(
  coef_temp ~ 1 + (1|gr(spp, cov = A_temp)), 
  data = mam_temp, family = gaussian(), 
  data2 = list(A_temp = A_temp),
  prior = c(
    prior(normal(0, 0.5), class =  Intercept),
    prior(exponential(0.5), class = sd, group = "spp"),
    prior(exponential(0.5), class = sigma)),
  control = list(adapt_delta = 0.95),
  chains = 4, cores = 4, iter = 2000, warmup = 500
)

precip_simple <- brm(
  coef_precip ~ 1 + (1|gr(spp, cov = A_precip)), 
  data = mam_precip, family = gaussian(), 
  data2 = list(A_precip = A_precip),
  # prior = c(
  #   prior(normal(0, 0.5), class =  Intercept),
  #   prior(exponential(0.5), class = sd, group = "spp"),
  #   prior(exponential(0.5), class = sigma)),
  control = list(adapt_delta = 0.95),
  chains = 4, cores = 4, iter = 2000, warmup = 500
)

summary(temp_simple)
summary(precip_simple)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 7. Phylogenetic signal ####

# Hypothesis test - 
hyp <- "sd_spp__Intercept^2 / (sd_spp__Intercept^2 + sigma^2) = 0"

## Doesn't look good - effectively 0
(hyp_temp <- hypothesis(temp_simple, hyp, class = NULL))
(hyp_precip <- hypothesis(precip_simple, hyp, class = NULL))

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 8. Gaussian process with rethinking - Temperature ####
  
# Correlation matrix as a distance matrix
Dmat <- cophenetic(mamMCC_temp)
spp_obs <- mam_temp$spp
Dmat <- Dmat[spp_obs, spp_obs]/max(Dmat)
  
# data as list
temp_list <- list(
    N_spp = nrow(mam_temp),
    C = as.vector(standardize(mam_temp$coef_temp)),
    O = as.numeric(as.factor(mam_temp$Order)),
    Dmat = Dmat
)
  
# Model with our exponential decay function
temp_OU_model <- ulam(alist(
  # The model
  C ~ multi_normal(mu, SIGMA),
  mu <- a + betaO*O,
  
  # Our adaptive prior - UL gaussian process L1
  matrix[N_spp,N_spp]:SIGMA <- cov_GPL1(Dmat, eta_sq, rho_sq, 0.01),
  
  # Our fixed priors
  a ~ normal(0,1),
  betaO ~ normal(0,2),
  eta_sq ~ half_normal(1,0.25),
  rho_sq ~ half_normal(3,0.25)),
  data = temp_list, 
  chains = 4, cores = 4)
  

# Look at the decay
post <- extract.samples(temp_OU_model)
plot(NULL, xlim = c(0, max(temp_list$Dmat)), ylim = c(0,5))

# posterior
for(i in 1:50){
  curve(post$eta_sq[i]*exp(-post$rho_sq[i]*x), add = TRUE, col = rangi2)
}

# prior
eta <- abs(rnorm(1e3, 1, 0.25))
rho <- abs(rnorm(1e3, 3, 0.25))
d_seq <- seq(0,1,length.out = 50)
K <- sapply(d_seq, function(x) eta*exp(-rho*x))

lines(d_seq, colMeans(K), lwd = 2)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 9. Simple PGLS for Temperature ####


tempsimple_gls <- gls(coef_temp ~ 1, 
                       data = mam_temp, 
                       correlation = corPagel(1,mamMCC_temp))
summary(tempsimple_gls)





