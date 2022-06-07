####################################################
##                                                ##
##     Global climate and population dynamics     ##
##                                                ##
##      Absolute temperature effect with GAM      ##
##                                                ##
##       Varying length of population records     ##
##                                                ##
##                 May 16th 2022                  ##
##                                                ##
####################################################

## Models regressing biome and life-history on absolute population responses
## from GAM ARMA models. Implementing the meta-regression framework for temperature

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
mam_temp_5yr <- mnanom_5km_GAM_5yr %>% 
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
  drop_na(litter, longevity, bodymass) %>% 
  filter(phylo != "Damaliscus_korrigum" &
           abs_temp <= 5)

# Long records
mam_temp_20yr <- mnanom_5km_GAM_20yr %>% 
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
  drop_na(litter, longevity, bodymass)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3. Phylogenetic covariance matrices ####

## Short Records
# Trim tree to right data
mamMCC_temp_5yr <- keep.tip(mamMCC_pruned, mam_temp_5yr$phylo) 

# Covariance matrix - Brownian motion model
A_temp_5yr <- ape::vcv.phylo(mamMCC_temp_5yr)

## Long Records
# Trim tree to right data
mamMCC_temp_20yr <- keep.tip(mamMCC_pruned, mam_temp_20yr$phylo) 

# Covariance matrix - Brownian motion model
A_temp_20yr <- ape::vcv.phylo(mamMCC_temp_20yr)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 4. >5 year records - Temperature models ####

# Gamma tuning plot
plot(density(rgamma(1000, shape = 3, scale = 0.7)))

tibble(sim = rgamma(1000, shape = 0.5, scale = 1.5)) %>% 
  ggplot(aes(x = sim)) + geom_density(fill = "lightgrey") +
  geom_density(data = mam_temp_5yr, 
               aes(x = abs_temp), 
               fill = "blue", alpha = 0.6) +
  theme_classic()

## Base model
set.seed(666)
temp_base_5yr <- brm(
  abs_temp ~ 1 + sample_size + (1|gr(phylo, cov = A_temp_5yr)) + (1| species),  
  data = mam_temp_5yr, 
  family = Gamma(link = "log"), 
  data2 = list(A_temp_5yr = A_temp_5yr),
  prior = c(
    prior(normal(0, 0.2), class =  Intercept),
    prior(normal(0, 0.2), class = b, coef = "sample_size"),
    prior(exponential(14), class = sd, group = "phylo"),
    prior(exponential(6), class = sd, group = "species"),
    prior(gamma(0.5,1.5), class = shape)),
  control = list(adapt_delta = 0.99,
                 max_treedepth = 15),
  chains = 3, cores = 3, iter = 4000, warmup = 2000
)

## Life-history univariate only
set.seed(666)
temp_lh_uni_5yr <- brm(
  abs_temp ~ 1 + longevity + bodymass + litter +
    sample_size + (1|gr(phylo, cov = A_temp_5yr)) + (1| species),  
  data = mam_temp_5yr, 
  family = Gamma(link = "log"), 
  data2 = list(A_temp_5yr = A_temp_5yr),
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

temp_base_5yr <- add_criterion(temp_base_5yr, criterion = c("loo","waic"))
temp_lh_uni_5yr <- add_criterion(temp_lh_uni_5yr, criterion = c("loo","waic"))

mod_comp_5yr <- as.data.frame(loo_compare(temp_base_5yr,temp_lh_uni_5yr, 
                                          criterion = "loo"))

mod_comp_5yr
save(mod_comp_5yr, file = "results/local_gamma_models/temperature_model_comparisons_5yr.RData")

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 5. >20 year records - Temperature models ####

tibble(sim = rgamma(1000, shape = 2, scale = 0.6)) %>% 
  ggplot(aes(x = sim)) + geom_density(fill = "lightgrey") +
  geom_density(data = mam_temp_20yr, 
               aes(x = abs_temp), 
               fill = "blue", alpha = 0.6) +
  theme_classic()

## Base model
set.seed(666)
temp_base_20yr <- brm(
  abs_temp ~ 1 + sample_size + (1|gr(phylo, cov = A_temp_20yr)) + (1| species),  
  data = mam_temp_20yr, 
  family = Gamma(link = "log"), 
  data2 = list(A_temp_20yr = A_temp_20yr),
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
temp_lh_uni_20yr <- brm(
  abs_temp ~ 1 + longevity + bodymass + litter +
    sample_size + (1|gr(phylo, cov = A_temp_20yr)) + (1| species),  
  data = mam_temp_20yr, 
  family = Gamma(link = "log"), 
  data2 = list(A_temp_20yr = A_temp_20yr),
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
temp_longevity_20yr <- brm(
  abs_temp ~ 1 + longevity +
    sample_size + (1|gr(phylo, cov = A_temp_20yr)) + (1| species),  
  data = mam_temp_20yr, 
  family = Gamma(link = "log"), 
  data2 = list(A_temp_20yr = A_temp_20yr),
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
temp_litter_20yr <- brm(
  abs_temp ~ 1 + litter +
    sample_size + (1|gr(phylo, cov = A_temp_20yr)) + (1| species),  
  data = mam_temp_20yr, 
  family = Gamma(link = "log"), 
  data2 = list(A_temp_20yr = A_temp_20yr),
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

temp_base_20yr <- add_criterion(temp_base_20yr, criterion = c("loo","waic"))
temp_lh_uni_20yr <- add_criterion(temp_lh_uni_20yr, criterion = c("loo","waic"))
temp_longevity_20yr <- add_criterion(temp_longevity_20yr, criterion = c("loo","waic"))
temp_litter_20yr <- add_criterion(temp_litter_20yr, criterion = c("loo","waic"))

mod_comp_20yr <- as.data.frame(loo_compare(temp_base_20yr, temp_lh_uni_20yr,
                                           temp_longevity_20yr, temp_litter_20yr,
                                           criterion = "loo"))

mod_comp_20yr
save(mod_comp_20yr, file = "results/local_gamma_models/temperature_model_comparisons_20yr.RData")


##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 6. >20 year records - No phylogeny ####

# Restricted sample size prevents meaningful phylogenetic analyses?

## Base model
set.seed(666)
temp_base_20yr_nophylo <- brm(
  abs_temp ~ 1 + sample_size + (1| species),  
  data = mam_temp_20yr, 
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
temp_lh_uni_20yr_nophylo <- brm(
  abs_temp ~ 1 + longevity + bodymass + litter +
    sample_size + (1| species),  
  data = mam_temp_20yr, 
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
temp_longevity_20yr_nophylo <- brm(
  abs_temp ~ 1 + longevity + 
    sample_size + (1| species),  
  data = mam_temp_20yr, 
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
temp_litter_20yr_nophylo <- brm(
  abs_temp ~ 1 + litter + 
    sample_size + (1| species),  
  data = mam_temp_20yr, 
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

temp_base_20yr_nophylo <- add_criterion(temp_base_20yr_nophylo, criterion = c("loo","waic"))
temp_lh_uni_20yr_nophylo <- add_criterion(temp_lh_uni_20yr_nophylo, criterion = c("loo","waic"))
temp_longevity_20yr_nophylo <- add_criterion(temp_longevity_20yr_nophylo, criterion = c("loo","waic"))
temp_litter_20yr_nophylo <- add_criterion(temp_litter_20yr_nophylo, criterion = c("loo","waic"))

mod_comp_20yr_nophylo <- as.data.frame(loo_compare(temp_base_20yr_nophylo, temp_lh_uni_20yr_nophylo,
                                                   temp_longevity_20yr_nophylo, temp_litter_20yr_nophylo,
                                                   criterion = "loo")) 

mod_comp_20yr_nophylo
save(mod_comp_20yr_nophylo, file = "results/local_gamma_models/temperature_model_comparisons_20yr.RData")

save(temp_litter_20yr_nophylo, file = "../temp_litter_20yr_nophylo.RData")

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 7. Model plots ####

# Colours
temp_colour <- "#d45371"
precip_colour <- "#30738e"

##_____________________________________________________________
## 7a. Model comparison##

# 5 year
load("results/local_gamma_models/temperature_model_comparisons_5yr.RData")

# Full set of predictor variables tested
base_predictors <- tibble(predictors = c("base", "longevity + bodymass + litter"),
                          temp_name = c("temp_base_5yr","temp_lh_uni_5yr"))

mod_comp_5yr %>% 
  mutate(model = rownames(.)) %>% 
  left_join(x = ., y = base_predictors, 
            by = c("model" = "temp_name")) %>% 
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
  save_as_image("results/local_gamma_models/temp_model_selection_5yr.png")

## 20 year
load("results/local_gamma_models/temperature_model_comparisons_20yr.RData")

# Full set of predictor variables tested
lh_predictors <- tibble(predictors = c("base", "longevity", "litter",
                                       "longevity + bodymass + litter"),
                        temp_name = c("temp_base_20yr_nophylo",
                                        "temp_longevity_20yr_nophylo",
                                        "temp_litter_20yr_nophylo",
                                        "temp_lh_uni_20yr_nophylo"))

mod_comp_20yr_nophylo %>% 
  mutate(model = rownames(.)) %>% 
  left_join(x = ., y = lh_predictors, 
            by = c("model" = "temp_name")) %>% 
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
  save_as_image("results/local_gamma_models/temp_model_selection_20yr.png")


##_____________________________________________________________
## 7a. Posterior prediction ##

load("../temp_litter_20yr_nophylo.RData")

# Binned litter data
litter_binned_data <- mam_temp_20yr %>% 
  mutate(lit_bin = as.numeric(as.character(cut(litter, 
                                               breaks = seq(-1.2, 2.8, by = 0.2), 
                                               labels = seq(-1.2, 2.6, by = 0.2))))) %>% 
  group_by(lit_bin) %>% 
  summarise(mntemp = mean(abs_temp),
            setemp = sd(abs_temp)/sqrt(n()),
            mnprecip = mean(abs_precip, na.rm = T),
            seprecip = sd(abs_precip, na.rm = T)/sqrt(n())) 

# prediction data
preddat_lit <- expand_grid(litter = seq(-1.2,3, length.out = 50),
                           sample_size = mean(mam_temp_20yr$sample_size),
                           phylo = unique(mam_temp_20yr$phylo)) %>% 
  mutate(species = phylo)

brms::stancode(temp_litter_20yr_nophylo) #  mu[n] = shape * exp(-(mu[n]))

# posterior predictions
temp_pred_lit <-  brms::posterior_predict(temp_litter_20yr_nophylo, 
                                            newdata = preddat_lit, 
                                            type = "response") 

# summary data for each longevity value
posterior_summary_litter <- bind_rows(lapply(unique(preddat_lit$litter), function(x){
  
  cpos = which(preddat_lit$litter == x) # vector of positions corresponding to where the values are equal to x
  
  # posterior mean 
  post_mn_temp = mean(temp_pred_lit[,cpos])
  
  # prediction intervals - 80% from rethinking package
  cQuant_temp = quantile(temp_pred_lit[,cpos], c(0.05, 0.95))
  
  # return data
  return(tibble(litter = x, 
                post_mn_temp = post_mn_temp,
                lwr_temp = cQuant_temp[1], upr_temp = cQuant_temp[2]))
}))


# The plot
lit_temp_bin <- ggplot(litter_binned_data, aes(x = lit_bin, y = mntemp)) +
  geom_errorbar(aes(ymax = mntemp + setemp, ymin = mntemp - setemp),
                width = 0.03, colour = temp_colour) +
  geom_point(size = 4, colour = temp_colour) +
  geom_smooth(data = posterior_summary_litter,
              aes(x = litter, y = post_mn_temp, 
                  ymax = upr_temp, ymin = lwr_temp),
              colour = "black", fill = temp_colour, alpha = 0.3,
              stat = "identity") +
  coord_cartesian(ylim = c(0,4)) +
  scale_size_continuous(range = c(2,8), guide = "none") +
  labs(x = "Standardised litter size", y = "|Temperature effect on abundance|",
       tag = "a)") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank()) 

load("../precip_longevity_plot_20yr.RData")
lon_precip_bin <- lon_precip_bin +
  labs(tag = "b)")


ggsave(lit_temp_bin + lon_precip_bin,
       filename = "plots/manuscript_figures/Supplementary figures/20yr_posteriorpred.jpeg",
       width = 24, height = 12, units = "cm", dpi = 1000)


