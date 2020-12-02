####################################################
##                                                ##
##     Global climate and population dynamics     ##
##                                                ##
##                 Precipitation                  ##
##                                                ##
##        brms Phylogenetic meta-regression       ##
##        life-history and spatial effects        ##
##                                                ##
##                Nov 13th 2020                   ##
##                                                ##
####################################################

## Investigating the effects of biome and life-history on population responses
## Implementing the meta-regression framework for precipitation

rm(list = ls())
options(width = 100)

library(tidyverse)
library(tidybayes)
library(ape)
library(brms)
library(rethinking)
library(nlme)
library(patchwork)
library(ggridges)
library(ggdist)
library(viridis)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 1. Load data ####

## Weather effects
load("data/pgr_weather/mnanom_5km.RData", verbose = TRUE)
glimpse(mnanom_5km)

## Phylogenetic
load("../Phlogenetics/mam_lpd_pruned.RData", verbose = TRUE)
str(mamMCC_pruned)

## Life history data
load("data/lifehistory.RData", verbose = TRUE)

## Species names to link data
load("../rawdata/GBIF_species_names_mamUPDATE.RData", verbose = TRUE)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 2. Mammal coefficient data ####

mam_precip <- mnanom_5km %>% 
  left_join(x = ., y = dplyr::select(lpd_gbif, Binomial, gbif.species.id), 
            by = "Binomial") %>% 
  left_join(x = ., y = dplyr::select(lifehistory, c(3,6:9)),
            by = "gbif.species.id") %>% 
  mutate(species = gsub(" ", "_", gbif_species),
         phylo = species,
         # z transform the coefficients
         coef_temp = as.numeric(scale(coef_temp)),
         coef_precip = as.numeric(scale(coef_precip)),
         # z transform absolute latitude for models
         lat = as.numeric(scale(abs(Latitude))),
         # observation-level term for residual variance (not sure if needed)
         OLRE = 1:n(),
         # iucn change  
         IUCNstatus = if_else(is.na(IUCNstatus) == T, "NotAss", IUCNstatus),
         sample_size = as.numeric(scale(log(n_obs)))) %>% 
  dplyr::select(id = ID, id_block = ID_block, n_obs, sample_size, 
                order = Order, species, phylo, 
                biome, lat, iucn = IUCNstatus, litter,
                longevity, bodymass, coef_temp, 
                coef_precip) %>% 
  drop_na(litter, longevity, bodymass) %>% 
  ## Restrict for Precipitation data
  drop_na(coef_precip)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3. Phylogenetic covariance matrix ####

# Trim tree to right data
mamMCC_precip <- keep.tip(mamMCC_pruned, mam_precip$phylo)

# Covariance matrix - Brownian motion model
A_precip <- ape::vcv.phylo(mamMCC_precip)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 4. Precipitation models ####

## Base model
set.seed(666)
precip_base <- brm(
  coef_precip ~ 1 + sample_size + (1|gr(phylo, cov = A_precip)) + (1| species),  
  data = mam_precip, family = gaussian(), 
  data2 = list(A_precip = A_precip),
  prior = c(
    prior(normal(0, 0.1), class =  Intercept),
    prior(normal(0, 0.2), class = b, coef = "sample_size"),
    prior(exponential(6), class = sd, group = "phylo"),
    prior(exponential(5), class = sd, group = "species")),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 4, cores = 4, iter = 2000, warmup = 300
)

#_______________________________________________________________________________
### 4a. Univariate predictor models

## Longevity
set.seed(666)
precip_lon <- brm(
  coef_precip ~ 1 + longevity + sample_size + (1|gr(phylo, cov = A_precip)) + (1| species),  
  data = mam_precip, family = gaussian(), 
  data2 = list(A_precip = A_precip),
  prior = c(
    prior(normal(0, 0.1), class =  Intercept),
    prior(normal(0, 0.2), class = b, coef = "sample_size"),
    prior(normal(0, 0.2), class = b, coef = "longevity"),
    prior(exponential(6), class = sd, group = "phylo"),
    prior(exponential(5), class = sd, group = "species")),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 4, cores = 4, iter = 2000, warmup = 300
)

## Bodymass
set.seed(666)
precip_bod <- brm(
  coef_precip ~ 1 + bodymass + sample_size + (1|gr(phylo, cov = A_precip)) + (1| species),  
  data = mam_precip, family = gaussian(), 
  data2 = list(A_precip = A_precip),
  prior = c(
    prior(normal(0, 0.1), class =  Intercept),
    prior(normal(0, 0.2), class = b, coef = "sample_size"),
    prior(normal(0, 0.2), class = b, coef = "bodymass"),
    prior(exponential(6), class = sd, group = "phylo"),
    prior(exponential(5), class = sd, group = "species")),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 4, cores = 4, iter = 2000, warmup = 300
)

## Litter size
set.seed(666)
precip_lit <- brm(
  coef_precip ~ 1 + litter + sample_size + (1|gr(phylo, cov = A_precip)) + (1| species),  
  data = mam_precip, family = gaussian(), 
  data2 = list(A_precip = A_precip),
  prior = c(
    prior(normal(0, 0.1), class =  Intercept),
    prior(normal(0, 0.2), class = b, coef = "sample_size"),
    prior(normal(0, 0.2), class = b, coef = "litter"),
    prior(exponential(6), class = sd, group = "phylo"),
    prior(exponential(5), class = sd, group = "species")),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 4, cores = 4, iter = 2000, warmup = 300
)

## Biome
set.seed(666)
precip_biome <- brm(
  coef_precip ~ 1 + biome + sample_size + (1|gr(phylo, cov = A_precip)) + (1| species),  
  data = mam_precip, family = gaussian(), 
  data2 = list(A_precip = A_precip),
  prior = c(
    prior(normal(0, 0.1), class =  Intercept),
    prior(normal(0, 0.2), class = b, coef = "sample_size"),
    prior(normal(0, 0.2), class = b),
    prior(exponential(6), class = sd, group = "phylo"),
    prior(exponential(5), class = sd, group = "species")),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 4, cores = 4, iter = 2000, warmup = 300
)

#_______________________________________________________________________________
### 4b. Univariate model comparisons

precip_base <- add_criterion(precip_base, criterion = c("loo", "waic"))
precip_lon <- add_criterion(precip_lon, criterion = c("loo", "waic"))
precip_bod <- add_criterion(precip_bod, criterion = c("loo", "waic"))
precip_lit <- add_criterion(precip_lit, criterion = c("loo", "waic"))
precip_biome <- add_criterion(precip_biome, criterion = c("loo", "waic"))

loo_compare(precip_base, precip_lon, precip_bod, precip_lit, precip_biome, criterion = "loo")

#_______________________________________________________________________________
### 4a. Multivariate predictor models with longevity

## All life-history - not very appropriate because of covariance
set.seed(666)
precip_lh <- brm(
  coef_precip ~ 1 + longevity + bodymass + litter + sample_size + (1|gr(phylo, cov = A_precip)) + (1| species),  
  data = mam_precip, family = gaussian(), 
  data2 = list(A_precip = A_precip),
  prior = c(
    prior(normal(0, 0.1), class =  Intercept),
    prior(normal(0, 0.2), class = b),
    prior(exponential(6), class = sd, group = "phylo"),
    prior(exponential(5), class = sd, group = "species")),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 4, cores = 4, iter = 2000, warmup = 300
)

## Longevity effect controlling for bodymass
set.seed(666)
precip_lonbm <- brm(
  coef_precip ~ 1 + longevity + bodymass + sample_size + (1|gr(phylo, cov = A_precip)) + (1| species),  
  data = mam_precip, family = gaussian(), 
  data2 = list(A_precip = A_precip),
  prior = c(
    prior(normal(0, 0.1), class =  Intercept),
    prior(normal(0, 0.2), class = b),
    prior(exponential(6), class = sd, group = "phylo"),
    prior(exponential(5), class = sd, group = "species")),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 4, cores = 4, iter = 2000, warmup = 300
)

## Longevity effect excluding sample size
set.seed(666)
precip_lon_ns <- brm(
  coef_precip ~ 1 + longevity + (1|gr(phylo, cov = A_precip)) + (1| species),  
  data = mam_precip, family = gaussian(), 
  data2 = list(A_precip = A_precip),
  prior = c(
    prior(normal(0, 0.1), class =  Intercept),
    prior(normal(0, 0.2), class = b),
    prior(exponential(6), class = sd, group = "phylo"),
    prior(exponential(5), class = sd, group = "species")),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 4, cores = 4, iter = 2000, warmup = 300
)

#_______________________________________________________________________________
### 4b. Full model comparisons

precip_lh <- add_criterion(precip_lh, criterion = c("loo", "waic"))
precip_lonbm <- add_criterion(precip_lonbm, criterion = c("loo", "waic"))
precip_lon_ns <- add_criterion(precip_lon_ns, criterion = c("loo", "waic"))

loo_compare(precip_base, precip_lon, precip_bod, 
            precip_lit, precip_biome, precip_lh,
            precip_lonbm, precip_lon_ns, criterion = "loo")

# Looks like best support is for just the model 
# with longevity and controlling for sample size

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 5. Posterior plots ####

precip_colour <- viridis(20, option = "D")[10]

#_______________________________________________________________________________
### 5a. Main model effects

maineff <- precip_lon %>%
  gather_draws(`b_[Ils].*`, regex=TRUE) %>% ### <<<<---------------------------------------------------------- COOL NEW FUNCTION
  ungroup() %>% 
  mutate(group = str_replace_all(.variable, c("b_" = "",
                                              "sd_" = "",
                                              "__Intercept" = ""))) %>% 
  ggplot(aes(y = group, x = .value)) + 
  geom_vline(xintercept = 0, size = 0.8) +
  stat_halfeye(fill = precip_colour) +
  scale_y_discrete(labels = c("Population intercept", "Longevity", "Sample size")) +
  labs(x = "Posterior precipitation effect", y = NULL) +
  theme_ridges(center_axis_labels = TRUE, grid = T, line_size = 0.3) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 12))

#_______________________________________________________________________________
### 5b. Variance effects

vareff <- precip_lon %>%
  gather_draws(`sd_.*|sigma`, regex=TRUE) %>%  
  ungroup() %>% 
  mutate(group = str_replace_all(.variable, c("b_" = "",
                                              "sd_" = "",
                                              "__Intercept" = ""))) %>%
  ggplot(aes(y = group, x = .value)) + 
  geom_vline(xintercept = 0, size = 0.8) +
  stat_halfeye(fill = precip_colour) +
  scale_y_discrete(labels = c("Phylogenetic covariance", 
                              "Population-level variance", 
                              "Species-level variance")) +
  labs(x = "Posterior estimate", y = NULL) +
  theme_ridges(center_axis_labels = TRUE, grid = T, line_size = 0.3) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 12))

#_______________________________________________________________________________
### 5c. Save

ggsave(maineff / vareff,
       filename = "plots/phylogenetic_regression/meta_regression_full/precip_effects_all.jpeg",
       width = 17, height = 18, units = "cm", dpi = 600)


##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 6. Posterior predictions ####

preddat <- expand_grid(longevity = seq(-2.1,2, length.out = 50),
                       sample_size = mean(mam_precip$sample_size),
                       phylo = unique(mam_precip$phylo)) %>% 
  mutate(species = phylo)

post_pred <- brms::posterior_predict(precip_lon, newdata = preddat) 

# pull out posterior summaries for each longevity
postpred_dat <- bind_rows(lapply(unique(preddat$longevity), function(x){
  
  cpos = which(preddat$longevity == x)
  
  # posterior mean and median
  post_mn = mean(post_pred[,cpos])
  post_md = median(post_pred[,cpos])
  
  # prediction intervals - 90% from rethinking package
  cPI = PI(post_pred[,cpos])
  cQuant = quantile(post_pred[,cpos], c(0.025, 0.975))
  
  return(tibble(longevity = x, 
                postmn = post_mn, postmd = post_md,
                lwrPI = cPI[1], uprPI = cPI[2], 
                lwr = cQuant[1], upr = cQuant[2]))
}))

# posterior simulation runs
postsimrums <- preddat %>%
  add_fitted_draws(precip_lon, n = 100) %>% 
  ungroup() %>% 
  group_by(longevity, .draw) %>% 
  summarise(mnpost = mean(.value))

ggplot(mam_precip, aes(x = longevity, y = coef_precip)) +
  geom_hline(yintercept = 0) +
  geom_point(aes(size = sample_size), alpha = 0.3) +
  geom_line(data = postsimrums, aes(y = mnpost, group = .draw), 
            alpha = 0.3, colour = precip_colour) +
  geom_line(data = postpred_dat, aes(y = postmn)) +
  scale_size_continuous(range = c(1,8), guide = F) +
  labs(x = "Standardised longevity", y = "Precipitation effect") +
  coord_cartesian(ylim = c(-2,2)) +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank()) +
  ggsave(filename = "plots/phylogenetic_regression/meta_regression_full/longevity_effect_precip.jpeg",
         width = 15, height = 15, units = "cm", dpi = 600)





