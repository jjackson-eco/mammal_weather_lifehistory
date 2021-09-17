####################################################
##                                                ##
##     Global climate and population dynamics     ##
##                                                ##
##     State-Space models for weather effects     ##
##                                                ##
##                Sep 15th 2021                   ##
##                                                ##
##       Christie Le Coeur & John Jackson         ##
##                                                ##
####################################################

## Expanding CLs work to implement state space models across all records in the
## database. See 'SSM_weather_CL.R' for the initial models.

rm(list = ls())
options(width = 100)

library(tidyverse)
library(jagsUI)
library(patchwork)

##__________________________________________________________________________________________________
#### 1. Load data ####

# mammal data
load("../rawdata/mammal.RData")
glimpse(mammal)

# annual weather anomaly - focus on just the mean anomaly in this script at a 5km range
mam_chelsa_annual <- readRDS("data/mam_chelsa_annual.RDS") %>% 
  filter(scale == "scale_5km") %>% 
  dplyr::select(ID,year, weather_scale = scale, mean_temp_anomaly, mean_precip_anomaly)
glimpse(mam_chelsa_annual)

# Species names to merge
load("../rawdata/GBIF_species_names_mamUPDATE.RData", verbose = TRUE)

##__________________________________________________________________________________________________
#### 2. Joining data ####

# linking to weather data and species names 
mammal_weather <- mammal %>% 
  left_join(., y = mam_chelsa_annual, by = c("ID", "year")) %>% 
  left_join(., y = dplyr::select(lpd_gbif, Binomial, gbif_species = gbif.species.tree),
            by = "Binomial") %>% 
  mutate(year_s = as.numeric(scale(year))) 

glimpse(mammal_weather)

##__________________________________________________________________________________________________
#### 3. State-Space models ####

## ID data
IDs <- unique(mammal_weather$ID_block)

## Results data
ssm_coefdat <- data.frame(ID_block = IDs, 
                          coef_temp_ssm = NA, 
                          coef_precip_ssm = NA)

ssm_preddat <- mammal_weather %>% 
  dplyr::select(ID_block, year, raw_abundance2) %>% 
  mutate(fitted_temp = NA, lwr_temp = NA, upr_temp = NA,
         fitted_precip = NA, lwr_precip = NA, upr_precip = NA) %>% 
  as.data.frame()

## State-space models
for (i in 440:length(IDs)){
  
  record <- mammal_weather[mammal_weather$ID_block==IDs[i],]
  print(paste0("######## " , i, ") ", record$ID_block[1], " ", record$Binomial[1], " ########"))
  
  ##____________________________________________________________________________
  #### 3a. Temperature ####
  
  # Specify model in BUGS
  cat(file="ssmtemp.jags", "
  model {
  
  # Priors and constraints
  logN[1] ~ dnorm(5.6, 0.01)           # Prior for initial population size
  b0 ~ dnorm(0, 0.001)                 # Prior for mean growth rate
  beta ~ dnorm(0, 0.001)               # Prior for slope parameter
  sigma.proc ~ dunif(0, 1)             # Prior for sd of state process
  sigma2.proc <- pow(sigma.proc, 2)
  tau.proc <- pow(sigma.proc, -2)
  sigma.obs ~ dunif(0, 1)              # Prior for sd of observation process
  sigma2.obs <- pow(sigma.obs, 2)
  tau.obs <- pow(sigma.obs, -2)
  
  # Likelihood
  # State process
  for (t in 1:(T-1)){
    r[t] <- b0 + beta * x[t] + epsilon[t]   # Linear model for the population growth rate
    epsilon[t] ~ dnorm(0, tau.proc)         # Random noise of the population growth rate
    logN[t+1] <- logN[t] + r[t]
  }
  # Observation process
  for (t in 1:T) {
    y[t] ~ dnorm(logN[t], tau.obs)
  }
  
  # Population sizes on real scale
  for (t in 1:T) {
    N[t] <- exp(logN[t])
  }
  }
  ")
  
  # Bundle data
  jags.data <- list(y=log(record$raw_abundance2), 
                    T=length(record$year), 
                    x=record$mean_temp_anomaly)
  
  # Initial values
  inits <- function(){list(sigma.proc=runif(1, 0, 1), 
                           sigma.obs=runif(1, 0, 1), 
                           b0=rnorm(1), 
                           logN=c(rnorm(1, 5.6, 0.1), 
                                  rep(NA,(length(record$year)-1))))}
  
  # Parameters monitored
  parameters <- c("r",  "sigma2.obs", "sigma2.proc", "N", "b0", "beta")
  
  # MCMC settings
  ni <- 200000
  nt <- 6
  nb <- 100000
  nc <- 3
  na <- 5000
  
  # Call JAGS from R (BRT 3 min)
  ssmtemptot <- jags(jags.data, inits, parameters, "ssmtemp.jags", 
                     n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, n.adapt=na)
  tempbetaSSM <- ssmtemptot$mean$beta
  
  ##____________________________________________________________________________
  #### 3b. Precipitation ####
  
  # Specify model in BUGS language
  cat(file="ssmprec.jags", "
  model {

  # Priors and constraints
  logN[1] ~ dnorm(5.6, 0.01)           # Prior for initial population size
  b0 ~ dnorm(0, 0.001)                 # Prior for mean growth rate
  #beta ~ dnorm(0, 0.001)               # Prior for slope parameter
  beta ~ dunif(-10, 10)               # Prior for slope parameter

  sigma.proc ~ dunif(0, 1)             # Prior for sd of state process
  sigma2.proc <- pow(sigma.proc, 2)
  tau.proc <- pow(sigma.proc, -2)
  sigma.obs ~ dunif(0, 1)              # Prior for sd of observation process
  sigma2.obs <- pow(sigma.obs, 2)
  tau.obs <- pow(sigma.obs, -2)

  # Likelihood
  # State process
  for (t in 1:(T-1)){
    r[t] <- b0 + beta * x[t] + epsilon[t]   # Linear model for the population growth rate
    epsilon[t] ~ dnorm(0, tau.proc)         # Random noise of the population growth rate
    logN[t+1] <- logN[t] + r[t]
  }
  # Observation process
  for (t in 1:T) {
    y[t] ~ dnorm(logN[t], tau.obs)
  }

  # Population sizes on real scale
  for (t in 1:T) {
    N[t] <- exp(logN[t])
  }
  }
  ")
  
  # Bundle data
  jags.data <- list(y=log(record$raw_abundance2),
                    T=length(record$year),
                    x=record$mean_precip_anomaly)
  
  # Initial values
  inits <- function(){list(sigma.proc=runif(1, 0, 1),
                           sigma.obs=runif(1, 0, 1),
                           b0 = rnorm(1),
                           logN=c(rnorm(1, 5.6, 0.1), 
                                  rep(NA,(length(record$year)-1))))}
  
  # Parameters monitored
  parameters <- c("r", "sigma2.obs", "sigma2.proc", "N", "b0", "beta")
  
  # Call JAGS from R (BRT 3 min)
  ssmprectot <- jags(jags.data, inits, 
                     parameters, "ssmprec.jags", 
                     n.chains=nc, n.thin=nt, 
                     n.iter=ni, n.burnin=nb, 
                     n.adapt=na)
  precbetaSSM <- ssmprectot$mean$beta
  
  ##____________________________________________________________________________
  #### 3c. Gather the data ####
  
  ## Coefficient data 
  ssm_coefdat[i,"coef_temp_ssm"] <- tempbetaSSM
  ssm_coefdat[i,"coef_precip_ssm"] <- precbetaSSM
  
  ## Year predictions
  ssm_preddat[which(ssm_preddat$ID_block == record$ID_block[1]), "fitted_temp"] <- 
    colMeans(ssmtemptot$sims.list$N, na.rm = T)
  ssm_preddat[which(ssm_preddat$ID_block == record$ID_block[1]), "lwr_temp"] <- 
    apply(ssmtemptot$sims.list$N, 2, FUN = function(x) quantile(x, 0.025))
  ssm_preddat[which(ssm_preddat$ID_block == record$ID_block[1]), "upr_temp"] <- 
    apply(ssmtemptot$sims.list$N, 2, FUN = function(x) quantile(x, 0.975))
  
  ssm_preddat[which(ssm_preddat$ID_block == record$ID_block[1]), "fitted_precip"] <- 
    colMeans(ssmprectot$sims.list$N, na.rm = T)
  ssm_preddat[which(ssm_preddat$ID_block == record$ID_block[1]), "lwr_precip"] <- 
    apply(ssmprectot$sims.list$N, 2, FUN = function(x) quantile(x, 0.025))
  ssm_preddat[which(ssm_preddat$ID_block == record$ID_block[1]), "upr_precip"] <- 
    apply(ssmprectot$sims.list$N, 2, FUN = function(x) quantile(x, 0.975))
  
  print(paste0("######## Your job is ", round(i/length(IDs), 3)*100, "% Complete ########"))
}

##__________________________________________________________________________________________________
#### 4. Saving the data ####

save(ssm_coefdat, ssm_preddat, file = "data/pgr_weather/mnanom_5km_ssm.RData")

##__________________________________________________________________________________________________
#### 5. In-sample prediction ####

load("data/mammal_analysis_data_GAM.RData")
load("data/pgr_weather/mnanom_5km_ssm.RData")

temp_colour <- "#d45371"
precip_colour <- "#30738e"

temp_insample <- ssm_preddat %>% 
  mutate(ln_raw = log(raw_abundance2), 
         ln_temp_fit = log(fitted_temp)) %>% 
  group_by(ID_block) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  ggplot(aes(x = ln_raw, y = ln_temp_fit)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(alpha = 0.07, colour = temp_colour, size = 0.8) +
  labs(x = "ln Abundance", y = "ln Fitted abundance", tag = "a)") +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank())

precip_insample <- ssm_preddat %>% 
  mutate(ln_raw = log(raw_abundance2), 
         ln_precip_fit = log(fitted_precip)) %>% 
  group_by(ID_block) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  ggplot(aes(x = ln_raw, y = ln_precip_fit)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(alpha = 0.07, colour = precip_colour, size = 0.8) +
  labs(x = "ln Abundance", y = "ln Fitted abundance", tag = "b)") +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank())

## Length of record and difference
temp_lengthp <- ssm_preddat %>% 
  mutate(ln_raw = log(raw_abundance2), 
         ln_temp_fit = log(fitted_temp)) %>% 
  group_by(ID_block) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  mutate(diff = ln_raw - ln_temp_fit) %>% 
  ggplot(aes(x = n, y = diff)) +
  geom_point(alpha = 0.1, colour = temp_colour) +
  labs(x = "Length of record", y = "Abundance prediction difference", tag = "c)") +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank())

precip_lengthp <- ssm_preddat %>% 
  mutate(ln_raw = log(raw_abundance2), 
         ln_precip_fit = log(fitted_precip)) %>% 
  group_by(ID_block) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  mutate(diff = ln_raw - ln_precip_fit) %>% 
  ggplot(aes(x = n, y = diff)) +
  geom_point(alpha = 0.1, colour = precip_colour) +
  labs(x = "Length of record", y = "Abundance prediction difference", tage = "d)") +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank())

ggsave("plots/weather_pop_growth/ssm_insample_prediction.jpeg",
       (temp_insample + precip_insample)/
       (temp_lengthp + precip_lengthp),
       width = 20, height = 20, units = "cm", dpi = 1200)

##__________________________________________________________________________________________________
#### 5. Method comparison ####

mam_coef_comp <- mam_coef %>% 
  left_join(x = ., y = ssm_coefdat, by = c("id_block" = "ID_block")) %>% 
  na.omit(coef_temp_ssm, coef_precip_ssm)

# Pearsons regression
temp_cor <- cor.test(mam_coef_comp$coef_temp_raw, mam_coef_comp$coef_temp_ssm)
precip_cor <- cor.test(mam_coef_comp$coef_precip_raw, mam_coef_comp$coef_precip_ssm)

temp_cor_plot <- ggplot(mam_coef_comp, aes(x = coef_temp_raw, coef_temp_ssm)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point(alpha = 0.7, size = 1.5, colour = temp_colour) +
  annotate('text', x = 0, y = -8,
           label = paste0("r = ", round(temp_cor$estimate, 2), ", p < 0.001")) +
  labs(x = "GAM temperature effect", y = "SSM temperature effect", tage = "a)") +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank())
  
precip_cor_plot <- ggplot(mam_coef_comp, aes(x = coef_precip_raw, coef_precip_ssm)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point(alpha = 0.7, size = 1.5, colour = precip_colour) +
  annotate('text', x = 3, y = -1.25,
           label = paste0("r = ", round(precip_cor$estimate, 2), ", p < 0.001")) +
  labs(x = "GAM precipitation effect", y = "SSM precipitation effect", tage = "b)") +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank())

ggsave("plots/weather_pop_growth/ssm_method_comparison.jpeg",
       temp_cor_plot + precip_cor_plot,
       width = 24, height = 10, units = "cm", dpi = 1200)



