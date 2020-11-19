####################################################
##                                                ##
##     Global climate and population dynamics     ##
##                                                ##
##                   Temperature                  ##
##                                                ##
##        brms Phylogenetic meta-regression       ##
##        life-history and spatial effects        ##
##                                                ##
##                 Nov 9th 2020                   ##
##                                                ##
####################################################

## Investigating the effects of biome and life-history on population responses
## Implementing the meta-regression framework for temperature

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

mam_temp <- mnanom_5km %>% 
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
  drop_na(litter, longevity, bodymass)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3. Phylogenetic covariance matrix ####

# Trim tree to right data
mamMCC_temp <- keep.tip(mamMCC_pruned, mam_temp$phylo) 

# Covariance matrix - Brownian motion model
A_temp <- ape::vcv.phylo(mamMCC_temp)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 4. Temperature models ####

## Base model
set.seed(666)
temp_base <- brm(
  coef_temp ~ 1 + sample_size + (1|gr(phylo, cov = A_temp)) + (1| species),  
  data = mam_temp, family = gaussian(),
  data2 = list(A_temp = A_temp),
  prior = c(
    prior(normal(0, 0.2), class =  Intercept),
    prior(normal(0, 0.5), class = b, coef = "sample_size"),
    prior(exponential(5), class = sd, group = "phylo"),
    prior(exponential(5), class = sd, group = "species")),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 4, cores = 4, iter = 3000, warmup = 500
)

#_______________________________________________________________________________
### 4a. Univariate predictor models

## Longevity
set.seed(666)
temp_lon <- brm(
  coef_temp ~ 1 + longevity + sample_size + (1|gr(phylo, cov = A_temp)) + (1| species),  
  data = mam_temp, family = gaussian(),
  data2 = list(A_temp = A_temp),
  prior = c(
    prior(normal(0, 0.2), class =  Intercept),
    prior(normal(0, 0.2), class = b, coef = "longevity"),
    prior(normal(0, 0.5), class = b, coef = "sample_size"),
    prior(exponential(5), class = sd, group = "phylo"),
    prior(exponential(5), class = sd, group = "species")),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 4, cores = 4, iter = 3000, warmup = 500
)

## Bodymass
set.seed(666)
temp_bod <- brm(
  coef_temp ~ 1 + bodymass + sample_size + (1|gr(phylo, cov = A_temp)) + (1| species),  
  data = mam_temp, family = gaussian(),
  data2 = list(A_temp = A_temp),
  prior = c(
    prior(normal(0, 0.2), class =  Intercept),
    prior(normal(0, 0.2), class = b, coef = "bodymass"),
    prior(normal(0, 0.5), class = b, coef = "sample_size"),
    prior(exponential(5), class = sd, group = "phylo"),
    prior(exponential(5), class = sd, group = "species")),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 4, cores = 4, iter = 3000, warmup = 500
)

## Litter size
set.seed(666)
temp_lit <- brm(
  coef_temp ~ 1 + litter + sample_size + (1|gr(phylo, cov = A_temp)) + (1| species),  
  data = mam_temp, family = gaussian(),
  data2 = list(A_temp = A_temp),
  prior = c(
    prior(normal(0, 0.2), class =  Intercept),
    prior(normal(0, 0.15), class = b, coef = "litter"),
    prior(normal(0, 0.7), class = b, coef = "sample_size"),
    prior(exponential(5), class = sd, group = "phylo"),
    prior(exponential(5), class = sd, group = "species")),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 4, cores = 4, iter = 3000, warmup = 500
)

## Biome
set.seed(666)
temp_biome <- brm(
  coef_temp ~ 1 + biome + sample_size + (1|gr(phylo, cov = A_temp)) + (1| species),  
  data = mam_temp, family = gaussian(),
  data2 = list(A_temp = A_temp),
  prior = c(
    prior(normal(0, 0.2), class =  Intercept),
    prior(normal(0, 0.2), class = b),
    prior(normal(0, 0.5), class = b, coef = "sample_size"),
    prior(exponential(5), class = sd, group = "phylo"),
    prior(exponential(5), class = sd, group = "species")),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 4, cores = 4, iter = 3000, warmup = 500
)

## Biome no sample size
set.seed(666)
temp_biome_ns <- brm(
  coef_temp ~ 1 + biome + (1|gr(phylo, cov = A_temp)) + (1| species),  
  data = mam_temp, family = gaussian(),
  data2 = list(A_temp = A_temp),
  prior = c(
    prior(normal(0, 0.2), class =  Intercept),
    prior(normal(0, 0.2), class = b),
    prior(exponential(5), class = sd, group = "phylo"),
    prior(exponential(5), class = sd, group = "species")),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 4, cores = 4, iter = 3000, warmup = 500
)

#_______________________________________________________________________________
### 4b. Univariate model comparisons

## Model comparisons
temp_base <- add_criterion(temp_base, criterion = c("loo","waic"))
temp_lon <- add_criterion(temp_lon, criterion = c("loo","waic"))
temp_bod <- add_criterion(temp_bod, criterion = c("loo","waic"))
temp_lit <- add_criterion(temp_lit, criterion = c("loo","waic"))
temp_biome <- add_criterion(temp_biome, criterion = c("loo","waic"))
temp_biome_ns <- add_criterion(temp_biome_ns, criterion = c("loo","waic"))

loo_compare(temp_base, temp_lon, temp_bod, 
            temp_lit, temp_biome, criterion = "loo")

#_______________________________________________________________________________
### 4c. Full models of interest

## Longevity and biome
set.seed(666)
temp_lon_biome <- brm(
  coef_temp ~ 1 + longevity + biome + sample_size + (1|gr(phylo, cov = A_temp)) + (1| species),  
  data = mam_temp, family = gaussian(),
  data2 = list(A_temp = A_temp),
  prior = c(
    prior(normal(0, 0.2), class =  Intercept),
    prior(normal(0, 0.2), class = b, coef = "longevity"),
    prior(normal(0, 0.2), class = b), # For the rest of the beta terms i.e. the biome effect
    prior(normal(0, 0.5), class = b, coef = "sample_size"),
    prior(exponential(5), class = sd, group = "phylo"),
    prior(exponential(5), class = sd, group = "species")),
  control = list(adapt_delta = 0.99,
                 max_treedepth = 15),
  chains = 4, cores = 4, iter = 3000, warmup = 500
)

# ## Longevity and biome and no sample size
# set.seed(666)
# temp_lon_biome_ns <- brm(
#   coef_temp ~ 1 + longevity + biome + (1|gr(phylo, cov = A_temp)) + (1| species),  
#   data = mam_temp, family = gaussian(),
#   data2 = list(A_temp = A_temp),
#   prior = c(
#     prior(normal(0, 0.2), class =  Intercept),
#     prior(normal(0, 0.2), class = b, coef = "longevity"),
#     prior(normal(0, 0.1), class = b), # For the rest of the beta terms i.e. the biome effect
#     prior(exponential(5), class = sd, group = "phylo"),
#     prior(exponential(5), class = sd, group = "species")),
#   control = list(adapt_delta = 0.99,
#                  max_treedepth = 15),
#   chains = 4, cores = 4, iter = 3000, warmup = 500
# )

## Litter size and biome
set.seed(666)
temp_lit_biome <- brm(
  coef_temp ~ 1 + litter + biome + sample_size + (1|gr(phylo, cov = A_temp)) + (1| species),  
  data = mam_temp, family = gaussian(),
  data2 = list(A_temp = A_temp),
  prior = c(
    prior(normal(0, 0.2), class =  Intercept),
    prior(normal(0, 0.2), class = b, coef = "litter"),
    prior(normal(0, 0.2), class = b), # For the rest of the beta terms i.e. the biome effect
    prior(normal(0, 0.5), class = b, coef = "sample_size"),
    prior(exponential(5), class = sd, group = "phylo"),
    prior(exponential(5), class = sd, group = "species")),
  control = list(adapt_delta = 0.99,
                 max_treedepth = 15),
  chains = 4, cores = 4, iter = 3000, warmup = 500
)

## Longevity and biome interaction- Why is the longevity effect negative
set.seed(666)
temp_lon_biome_interaction <- brm(
  coef_temp ~ 1 + longevity*biome + sample_size + (1|gr(phylo, cov = A_temp)) + (1| species),  
  data = mam_temp, family = gaussian(),
  data2 = list(A_temp = A_temp),
  prior = c(
    prior(normal(0, 0.2), class =  Intercept),
    prior(normal(0, 0.1), class = b, coef = "longevity"),
    prior(normal(0, 0.1), class = b), # For the rest of the beta terms i.e. the biome effect and interaction - MORE CONSERVATIVE THAN BEFORE
    prior(normal(0, 0.5), class = b, coef = "sample_size"),
    prior(exponential(5), class = sd, group = "phylo"),
    prior(exponential(5), class = sd, group = "species")),
  control = list(adapt_delta = 0.99,
                 max_treedepth = 15),
  chains = 4, cores = 4, iter = 3000, warmup = 500
)


#_______________________________________________________________________________
### 4d. Full model comparisons

temp_lon_biome <- add_criterion(temp_lon_biome, criterion = c("loo","waic"))
temp_lon_biome_ns <- add_criterion(temp_lon_biome_ns, criterion = c("loo","waic"))
temp_lit_biome <- add_criterion(temp_lit_biome, criterion = c("loo","waic"))
temp_lon_biome_interaction <- add_criterion(temp_lon_biome_interaction, criterion = c("loo","waic"))

loo_compare(temp_base, temp_lon, 
            temp_bod, temp_lit, 
            temp_biome, temp_biome_ns,
            temp_lon_biome, #temp_lon_biome_ns,
            temp_lit_biome, temp_lon_biome_interaction,
            criterion = "loo")

## Really looks like the longevity effect is strongest + partly the biome - Plot it

#_______________________________________________________________________________
### 4e. Save pertinent posterior samples

set.seed(666)
post_tlb <- posterior_samples(temp_lon_biome)[1:10,1:18]
save(post_tlb, file = "data/temp_lon_biome_posteriorsamples.RData")

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 5. Posterior plots ####

## Getting right data separate
mod_post <- post_tlb[,-(grep("biome", colnames(post_tlb)))]
biome_post <- post_tlb[,grep("biome", colnames(post_tlb))]

temp_colour <- viridis(20, option = "C")[13]

#_______________________________________________________________________________
### 5a. Main model effects

maineff <- temp_lon_biome %>%
  gather_draws(`b_[Ils].*`, regex=TRUE) %>% ### <<<<---------------------------------------------------------- COOL NEW FUNCTION
  ungroup() %>% 
  mutate(group = str_replace_all(.variable, c("b_" = "",
                                              "sd_" = "",
                                              "__Intercept" = ""))) %>% 
  ggplot(aes(y = group, x = .value)) + 
  geom_vline(xintercept = 0, size = 0.8) +
  stat_halfeye(fill = temp_colour) +
  scale_y_discrete(labels = c("Population intercept", "Longevity", "Sample size")) +
  labs(x = "Posterior temperature coefficient", y = NULL) +
  theme_ridges(center_axis_labels = TRUE, grid = T, line_size = 0.3) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 12))

#_______________________________________________________________________________
### 5b. Variance effects

vareff <- temp_lon_biome %>%
  gather_draws(`sd_.*|sigma`, regex=TRUE) %>%  
  ungroup() %>% 
  mutate(group = str_replace_all(.variable, c("b_" = "",
                                              "sd_" = "",
                                              "__Intercept" = ""))) %>%
  ggplot(aes(y = group, x = .value)) + 
  geom_vline(xintercept = 0, size = 0.8) +
  stat_halfeye(fill = temp_colour) +
  scale_y_discrete(labels = c("Phylogenetic covariance", 
                              "Population-level variance", 
                              "Species-level variance")) +
  labs(x = "Posterior estimate", y = NULL) +
  theme_ridges(center_axis_labels = TRUE, grid = T, line_size = 0.3) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 12))

#_______________________________________________________________________________
### 5c. Biome effect

# Summary factor levels
fct_biomes <- temp_lon_biome %>%
  gather_draws(`b_biome.*`, regex = TRUE) %>%
  mean_qi() %>% 
  ungroup() %>% 
  mutate(group = str_replace(.variable, "b_biome", "")) %>% 
  arrange(.value) %>% 
  pull(group)

fct_biome_labs <- c("Temperate grassland", "Flooded grassland", 
                    "Mediterranean",  "Montane", "Tundra", 
                    "Tropical broadleaf", "Temperate broadleaf", 
                    "Desert", "Tropical coniferous", 
                    "Temperate coniferous",  "Tropical grassland", 
                    "Tropical moist")

# Point intervals from the posterior summary
# temp_lon_biome %>%
#   gather_draws(`b_biome.*`, regex = TRUE) %>%
#   median_qi(.width = c(.95, .66)) %>% 
#   ungroup() %>% 
#   mutate(group = factor(str_replace(.variable, "b_biome", ""),
#                         levels = fct_biomes)) %>%
#   ggplot(aes(y = group, x = .value, xmin = .lower, xmax = .upper)) +
#   geom_pointinterval() 

biomeeff <- temp_lon_biome %>%
  gather_draws(`b_biome.*`, regex=TRUE) %>%
  ungroup() %>% 
  mutate(group = factor(str_replace(.variable, "b_biome", ""),
                        levels = fct_biomes)) %>%
  ggplot(aes(y = group, x = .value)) +
  geom_vline(xintercept = 0, size = 0.8) +
  stat_halfeye(fill = temp_colour) +
  scale_y_discrete(labels = fct_biome_labs) +
  labs(x = "Posterior temperature coefficient", y = NULL) +
  theme_ridges(center_axis_labels = TRUE, grid = T, line_size = 0.3) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 12))

## Save 
ggsave((maineff / vareff) | biomeeff,
       filename = "plots/phylogenetic_regression/meta_regression_full/temp_effects_all.jpeg",
       width = 25, height = 23, units = "cm", dpi = 600)


##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 6. Posterior predictions ####

preddat <- expand_grid(longevity = seq(-2.1,2, length.out = 50),
                       biome = "Temperate broadleaf and mixed forests",
                       sample_size = mean(mam_temp$sample_size),
                       phylo = unique(mam_temp$phylo)) %>% 
  mutate(species = phylo)

post_pred <- brms::posterior_predict(temp_lon_biome, newdata = preddat) 

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
  add_fitted_draws(temp_lon_biome, n = 100) %>% 
  ungroup() %>% 
  group_by(longevity, .draw) %>% 
  summarise(mnpost = mean(.value))

ggplot(mam_temp, aes(x = longevity, y = coef_temp)) +
  geom_hline(yintercept = 0) +
  geom_point(aes(size = sample_size), alpha = 0.3) +
  geom_line(data = postsimrums, aes(y = mnpost, group = .draw), 
            alpha = 0.3, colour = temp_colour) +
  geom_line(data = postpred_dat, aes(y = postmn)) +
  scale_size_continuous(range = c(1,8), guide = F) +
  labs(x = "Standardised longevity", y = "Temperature coefficient") +
  coord_cartesian(ylim = c(-2.5,2.5)) +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank()) +
  ggsave(filename = "plots/phylogenetic_regression/meta_regression_full/longevity_effect_temp.jpeg",
         width = 15, height = 15, units = "cm", dpi = 800)

ggplot(postpred_dat, aes(x = longevity, y = postmn)) +
  geom_smooth(stat = "identity", 
              aes(ymin = lwrPI, ymax = uprPI),
              alpha = 0.7)
