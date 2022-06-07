####################################################
##                                                ##
##     Global climate and population dynamics     ##
##                                                ##
##      Absolute precipitation effect with GAM    ##
##                                                ##
##       Varying length of population records     ##
##                                                ##
##                 Feb 14th 2022                  ##
##                                                ##
####################################################

## Models regressing biome and life-history on absolute population responses
## from GAM ARMA models. Implementing the meta-regression framework for precipitation

rm(list = ls())
options(width = 100)

## Change the .libPaths and R_LIBS_USER to the right thing if you're on a uni computer
if(Sys.info()["nodename"] == "BIO-W-LT-000083" |
   Sys.info()["nodename"] == "BIO-W-DT-02108") {
  .libPaths("C:/Users/zool2541/R-4.1.1/library/")
  .libPaths("!\\\\zoo-suitcase/home$/zool2541/My Documents/R/win-library/4.1")}

## Packages
library(tidyverse)
library(tidybayes)
library(ape)
library(brms)
library(nlme)
library(patchwork)
library(ggridges)
library(ggdist)
library(viridis)
library(flextable)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 1. Load data ####

## Weather effects
load("data/pgr_weather/mnanom_5km_GAM_varied_record_length.RData", verbose = TRUE)

## Phylogenetic
load("../Phlogenetics/mam_lpd_pruned.RData", verbose = TRUE)
str(mamMCC_pruned)

## Life history data
load("data/lifehistory.RData", verbose = TRUE)

## Species names to link data
load("../rawdata/GBIF_species_names_mamUPDATE.RData", verbose = TRUE)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 2. Mammal coefficient data ####

# Short records
mam_precip_5yr <- mnanom_5km_GAM_5yr %>% 
  ungroup() %>%  ## <- The groups ruin the z-transformations
  left_join(x = ., y = dplyr::select(lpd_gbif, Binomial, gbif.species.id), 
            by = "Binomial") %>% 
  left_join(x = ., y = dplyr::select(lifehistory, c(3,6:9)),
            by = "gbif.species.id") %>% 
  mutate(species = gsub(" ", "_", gbif_species),
         phylo = species,
         # z transform the coefficients and take absolute value
         abs_temp = abs(as.numeric(scale(coef_temp))),
         abs_precip = abs(as.numeric(scale(coef_precip))),
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
                longevity, bodymass, abs_temp, 
                abs_precip) %>% 
  drop_na(litter, longevity, bodymass, abs_precip) %>% 
  filter(phylo != "Damaliscus_korrigum" &
           abs_precip <= 5)

# Long records
mam_precip_20yr <- mnanom_5km_GAM_20yr %>% 
  ungroup() %>%  ## <- The groups ruin the z-transformations
  left_join(x = ., y = dplyr::select(lpd_gbif, Binomial, gbif.species.id), 
            by = "Binomial") %>% 
  left_join(x = ., y = dplyr::select(lifehistory, c(3,6:9)),
            by = "gbif.species.id") %>% 
  mutate(species = gsub(" ", "_", gbif_species),
         phylo = species,
         # z transform the coefficients
         abs_temp = abs(as.numeric(scale(coef_temp))),
         abs_precip = abs(as.numeric(scale(coef_precip))),
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
                longevity, bodymass, abs_temp, 
                abs_precip) %>% 
  drop_na(litter, longevity, bodymass, abs_precip)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3. Phylogenetic covariance matrices ####

## Short Records
# Trim tree to right data
mamMCC_precip_5yr <- keep.tip(mamMCC_pruned, mam_precip_5yr$phylo) 

# Covariance matrix - Brownian motion model
A_precip_5yr <- ape::vcv.phylo(mamMCC_precip_5yr)

## Long Records
# Trim tree to right data
mamMCC_precip_20yr <- keep.tip(mamMCC_pruned, mam_precip_20yr$phylo) 

# Covariance matrix - Brownian motion model
A_precip_20yr <- ape::vcv.phylo(mamMCC_precip_20yr)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 4. >5 year records - Precipitation models ####

# Gamma tuning plot
plot(density(rgamma(1000, shape = 3, scale = 0.7)))

tibble(sim = rgamma(1000, shape = 3, scale = 0.7)) %>% 
  ggplot(aes(x = sim)) + geom_density(fill = "lightgrey") +
  geom_density(data = mam_precip_5yr, 
               aes(x = abs_precip), 
               fill = "blue", alpha = 0.6) +
  theme_classic()

## Base model
set.seed(666)
precip_base_5yr <- brm(
  abs_precip ~ 1 + sample_size + (1|gr(phylo, cov = A_precip_5yr)) + (1| species),  
  data = mam_precip_5yr, 
  family = Gamma(link = "log"), 
  data2 = list(A_precip_5yr = A_precip_5yr),
  prior = c(
    prior(normal(0, 0.2), class =  Intercept),
    prior(normal(0, 0.2), class = b, coef = "sample_size"),
    prior(exponential(14), class = sd, group = "phylo"),
    prior(exponential(6), class = sd, group = "species"),
    prior(gamma(3,0.8), class = shape)),
  control = list(adapt_delta = 0.99,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 4000, warmup = 2000
)

## Life-history univariate only
set.seed(666)
precip_lh_uni_5yr <- brm(
  abs_precip ~ 1 + longevity + bodymass + litter +
    sample_size + (1|gr(phylo, cov = A_precip_5yr)) + (1| species),  
  data = mam_precip_5yr, 
  family = Gamma(link = "log"), 
  data2 = list(A_precip_5yr = A_precip_5yr),
  prior = c(
    prior(normal(0, 0.2), class =  Intercept),
    prior(normal(0, 0.2), class = b),
    prior(normal(0, 0.2), class = b, coef = "sample_size"),
    prior(exponential(14), class = sd, group = "phylo"),
    prior(exponential(6), class = sd, group = "species"),
    prior(gamma(3,0.8), class = shape)),
  control = list(adapt_delta = 0.99,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 4000, warmup = 2000
)

#_______________________________________________________________________________
### 4b. Full model comparisons

precip_base_5yr <- add_criterion(precip_base_5yr, criterion = c("loo","waic"))
precip_lh_uni_5yr <- add_criterion(precip_lh_uni_5yr, criterion = c("loo","waic"))

mod_comp_5yr <- as.data.frame(loo_compare(precip_base_5yr,precip_lh_uni_5yr, 
                                          criterion = "loo"))

mod_comp_5yr
save(mod_comp_5yr, file = "results/local_gamma_models/precipitation_model_comparisons_5yr.RData")

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 5. >20 year records - Precipitation models ####

tibble(sim = rgamma(1000, shape = 2, scale = 0.6)) %>% 
  ggplot(aes(x = sim)) + geom_density(fill = "lightgrey") +
  geom_density(data = mam_precip_20yr, 
               aes(x = abs_precip), 
               fill = "blue", alpha = 0.6) +
  theme_classic()

## Base model
set.seed(666)
precip_base_20yr <- brm(
  abs_precip ~ 1 + sample_size + (1|gr(phylo, cov = A_precip_20yr)) + (1| species),  
  data = mam_precip_20yr, 
  family = Gamma(link = "log"), 
  data2 = list(A_precip_20yr = A_precip_20yr),
  prior = c(
    prior(normal(0, 0.5), class =  Intercept),
    prior(normal(0, 0.5), class = b, coef = "sample_size"),
    prior(exponential(11), class = sd, group = "phylo"),
    prior(exponential(2), class = sd, group = "species"),
    prior(gamma(2,0.6), class = shape)),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 4000, warmup = 2000
)

## Life-history univariate only
set.seed(666)
precip_lh_uni_20yr <- brm(
  abs_precip ~ 1 + longevity + bodymass + litter +
    sample_size + (1|gr(phylo, cov = A_precip_20yr)) + (1| species),  
  data = mam_precip_20yr, 
  family = Gamma(link = "log"), 
  data2 = list(A_precip_20yr = A_precip_20yr),
  prior = c(
    prior(normal(0, 0.3), class =  Intercept),
    prior(normal(0, 0.2), class = b),
    prior(normal(0, 0.2), class = b, coef = "sample_size"),
    prior(exponential(11), class = sd, group = "phylo"),
    prior(exponential(2), class = sd, group = "species"),
    prior(gamma(2,0.6), class = shape)),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 4000, warmup = 2000
)

## Longevity
set.seed(666)
precip_longevity_20yr <- brm(
  abs_precip ~ 1 + longevity +
    sample_size + (1|gr(phylo, cov = A_precip_20yr)) + (1| species),  
  data = mam_precip_20yr, 
  family = Gamma(link = "log"), 
  data2 = list(A_precip_20yr = A_precip_20yr),
  prior = c(
    prior(normal(0, 0.3), class =  Intercept),
    prior(normal(0, 0.2), class = b),
    prior(normal(0, 0.2), class = b, coef = "sample_size"),
    prior(exponential(11), class = sd, group = "phylo"),
    prior(exponential(2), class = sd, group = "species"),
    prior(gamma(2,0.6), class = shape)),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 4000, warmup = 2000
)

## Litter
set.seed(666)
precip_litter_20yr <- brm(
  abs_precip ~ 1 + litter +
    sample_size + (1|gr(phylo, cov = A_precip_20yr)) + (1| species),  
  data = mam_precip_20yr, 
  family = Gamma(link = "log"), 
  data2 = list(A_precip_20yr = A_precip_20yr),
  prior = c(
    prior(normal(0, 0.3), class =  Intercept),
    prior(normal(0, 0.2), class = b),
    prior(normal(0, 0.2), class = b, coef = "sample_size"),
    prior(exponential(11), class = sd, group = "phylo"),
    prior(exponential(2), class = sd, group = "species"),
    prior(gamma(2,0.6), class = shape)),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 4000, warmup = 2000
)


#_______________________________________________________________________________
### 5b. Full model comparisons

precip_base_20yr <- add_criterion(precip_base_20yr, criterion = c("loo","waic"))
precip_lh_uni_20yr <- add_criterion(precip_lh_uni_20yr, criterion = c("loo","waic"))
precip_longevity_20yr <- add_criterion(precip_longevity_20yr, criterion = c("loo","waic"))
precip_litter_20yr <- add_criterion(precip_litter_20yr, criterion = c("loo","waic"))

mod_comp_20yr <- as.data.frame(loo_compare(precip_base_20yr, precip_lh_uni_20yr,
                                           precip_longevity_20yr, precip_litter_20yr,
                                           criterion = "loo"))

mod_comp_20yr
save(mod_comp_20yr, file = "results/local_gamma_models/precipitation_model_comparisons_20yr.RData")


##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 6. >20 year records - No phylogeny ####

# Restricted sample size prevents meaningful phylogenetic analyses?

## Base model
set.seed(666)
precip_base_20yr_nophylo <- brm(
  abs_precip ~ 1 + sample_size + (1| species),  
  data = mam_precip_20yr, 
  family = Gamma(link = "log"), 
  prior = c(
    prior(normal(0, 0.8), class =  Intercept),
    prior(normal(0, 0.8), class = b, coef = "sample_size"),
    prior(exponential(2), class = sd, group = "species"),
    prior(gamma(2,0.5), class = shape)),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 4000, warmup = 2000
)

## Life-history univariate only
set.seed(666)
precip_lh_uni_20yr_nophylo <- brm(
  abs_precip ~ 1 + longevity + bodymass + litter +
    sample_size + (1| species),  
  data = mam_precip_20yr, 
  family = Gamma(link = "log"), 
  prior = c(
    prior(normal(0, 0.3), class =  Intercept),
    prior(normal(0, 0.2), class = b),
    prior(normal(0, 0.2), class = b, coef = "sample_size"),
    prior(exponential(2), class = sd, group = "species"),
    prior(gamma(2,0.6), class = shape)),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 4000, warmup = 2000
)

## Longevity
set.seed(666)
precip_longevity_20yr_nophylo <- brm(
  abs_precip ~ 1 + longevity + 
    sample_size + (1| species),  
  data = mam_precip_20yr, 
  family = Gamma(link = "log"), 
  prior = c(
    prior(normal(0, 0.3), class =  Intercept),
    prior(normal(0, 0.2), class = b),
    prior(normal(0, 0.2), class = b, coef = "sample_size"),
    prior(exponential(2), class = sd, group = "species"),
    prior(gamma(2,0.6), class = shape)),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 4000, warmup = 2000
)

## Litter
set.seed(666)
precip_litter_20yr_nophylo <- brm(
  abs_precip ~ 1 + litter + 
    sample_size + (1| species),  
  data = mam_precip_20yr, 
  family = Gamma(link = "log"), 
  prior = c(
    prior(normal(0, 0.3), class =  Intercept),
    prior(normal(0, 0.2), class = b),
    prior(normal(0, 0.2), class = b, coef = "sample_size"),
    prior(exponential(2), class = sd, group = "species"),
    prior(gamma(2,0.6), class = shape)),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 4000, warmup = 2000
)


#_______________________________________________________________________________
### 6b. Full model comparisons

precip_base_20yr_nophylo <- add_criterion(precip_base_20yr_nophylo, criterion = c("loo","waic"))
precip_lh_uni_20yr_nophylo <- add_criterion(precip_lh_uni_20yr_nophylo, criterion = c("loo","waic"))
precip_longevity_20yr_nophylo <- add_criterion(precip_longevity_20yr_nophylo, criterion = c("loo","waic"))
precip_litter_20yr_nophylo <- add_criterion(precip_litter_20yr_nophylo, criterion = c("loo","waic"))

mod_comp_20yr_nophylo <- as.data.frame(loo_compare(precip_base_20yr_nophylo, precip_lh_uni_20yr_nophylo,
                                                   precip_longevity_20yr_nophylo, precip_litter_20yr_nophylo,
                                           criterion = "loo")) 

mod_comp_20yr_nophylo

save(precip_longevity_20yr_nophylo, file = "../precip_longevity_20yr_nophylo.RData")

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 7. Model plots ####

# Colours
temp_colour <- "#d45371"
precip_colour <- "#30738e"

##_____________________________________________________________
## 7a. Model comparison##

# 5 year
load("results/local_gamma_models/precipitation_model_comparisons_5yr.RData")

# Full set of predictor variables tested
base_predictors <- tibble(predictors = c("base", "longevity + bodymass + litter"),
                         precip_name = c("precip_base_5yr","precip_lh_uni_5yr"))

mod_comp_5yr %>% 
  mutate(model = rownames(.)) %>% 
  left_join(x = ., y = base_predictors, 
            by = c("model" = "precip_name")) %>% 
  dplyr::select(model, predictors, elpd_loo, se_elpd_loo, elpd_diff, se_diff, looic) %>% 
  flextable(cwidth = 1.5) %>% 
  set_header_labels(model = "Model",
                    predictos = "Life-history predictors",
                    elpd_loo = "LOO elpd",
                    se_elpd_loo = "LOO elpd error",
                    elpd_diff = "elpd difference",
                    se_diff = "elpd error difference",
                    looic = "LOO information criterion") %>% 
  colformat_double(digits = 2) %>% 
  save_as_image("results/local_gamma_models/precip_model_selection_5yr.png")

## 20 year
load("results/local_gamma_models/precipitation_model_comparisons_20yr.RData")

# Full set of predictor variables tested
lh_predictors <- tibble(predictors = c("base", "longevity", "litter",
                                         "longevity + bodymass + litter"),
                          precip_name = c("precip_base_20yr_nophylo",
                                          "precip_longevity_20yr_nophylo",
                                          "precip_litter_20yr_nophylo",
                                          "precip_lh_uni_20yr_nophylo"))

mod_comp_20yr_nophylo %>% 
  mutate(model = rownames(.)) %>% 
  left_join(x = ., y = lh_predictors, 
            by = c("model" = "precip_name")) %>% 
  dplyr::select(model, predictors, elpd_loo, se_elpd_loo, elpd_diff, se_diff, looic) %>% 
  flextable(cwidth = 1.5) %>% 
  set_header_labels(model = "Model",
                    predictos = "Life-history predictors",
                    elpd_loo = "LOO elpd",
                    se_elpd_loo = "LOO elpd error",
                    elpd_diff = "elpd difference",
                    se_diff = "elpd error difference",
                    looic = "LOO information criterion") %>% 
  colformat_double(digits = 2) %>% 
  save_as_image("results/local_gamma_models/precip_model_selection_20yr.png")


##_____________________________________________________________
## 7a. Posterior prediction ##

load("../precip_longevity_20yr_nophylo.RData")

# Binned longevity data
longevity_binned_data <- mam_precip_20yr %>% 
  mutate(lon_bin = as.numeric(as.character(cut(longevity, 
                                               breaks = c(seq(-2.1, 2, by = 0.2), 2), 
                                               labels = seq(-2.1, 2, by = 0.2))))) %>% 
  group_by(lon_bin) %>% 
  summarise(mntemp = mean(abs_temp),
            setemp = sd(abs_temp)/sqrt(n()),
            mnprecip = mean(abs_precip, na.rm = T),
            seprecip = sd(abs_precip, na.rm = T)/sqrt(n())) 

# prediction data
preddat_lon <- expand_grid(longevity = seq(-2.1,2, length.out = 50),
                           sample_size = mean(mam_precip_20yr$sample_size),
                           phylo = unique(mam_precip_20yr$phylo)) %>% 
  mutate(species = phylo)

brms::stancode(precip_longevity_20yr_nophylo) #  mu[n] = shape * exp(-(mu[n]))

# posterior predictions
precip_pred_lon <-  brms::posterior_predict(precip_longevity_20yr_nophylo, 
                                            newdata = preddat_lon, 
                                            type = "response") 

# summary data for each longevity value
posterior_summary_longevity <- bind_rows(lapply(unique(preddat_lon$longevity), function(x){
  
  cpos = which(preddat_lon$longevity == x) # vector of positions corresponding to where the values are equal to x
  
  # posterior mean 
  post_mn_precip = mean(precip_pred_lon[,cpos])
  
  # prediction intervals - 80% from rethinking package
  cQuant_precip = quantile(precip_pred_lon[,cpos], c(0.05, 0.95))
  
  # return data
  return(tibble(longevity = x, 
                post_mn_precip = post_mn_precip,
                lwr_precip = cQuant_precip[1], upr_precip = cQuant_precip[2]))
}))


# The plot
lon_precip_bin <- ggplot(longevity_binned_data, aes(x = lon_bin, y = mnprecip)) +
  geom_errorbar(aes(ymax = mnprecip + seprecip, ymin = mnprecip - seprecip),
                width = 0.03, colour = precip_colour) +
  geom_point(size = 4, colour = precip_colour) +
  geom_smooth(data = posterior_summary_longevity,
              aes(x = longevity, y = post_mn_precip, 
                  ymax = upr_precip, ymin = lwr_precip),
              colour = "black", fill = precip_colour, alpha = 0.3,
              stat = "identity") +
  coord_cartesian(ylim = c(0,4)) +
  scale_size_continuous(range = c(2,8), guide = "none") +
  labs(x = "Standardised longevity", y = "|Precipitation effect on abundance|",
       tag = "b)") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank()) 

save(lon_precip_bin, file = "../precip_longevity_plot_20yr.RData")




