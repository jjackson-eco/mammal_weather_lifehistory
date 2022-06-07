#################################################################
##                                                             ##
##           Global climate and population dynamics            ##
##                                                             ##
##       COMADRE life history data and GAM coefficients        ##
##                                                             ##
##                    Linear regression                        ##
##                                                             ##
##                     March 22nd 2022                         ##
##                                                             ##
#################################################################

## Investigating the association between life history parameters from COMADRE and
## absolute weather coefficients

rm(list = ls())
options(width = 100)

## Packages
library(brms)
library(tidybayes)
library(tidyverse)
library(patchwork)
library(viridis)
library(psych)

# colours
temp_colour <- "#d45371"
precip_colour <- "#30738e"

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 1. Load data ####

## Weather effects
load("data/pgr_weather/mnanom_5km_GAM.RData", verbose = TRUE)
glimpse(mnanom_5km_GAM)

## COMADRE Life history data
load("../rawdata/mam_comadre.RData", verbose = TRUE)

## Species names to link data
load("../rawdata/GBIF_species_names_mamUPDATE.RData", verbose = TRUE)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 2. Coefficient data ####

# Species-level data from comadre
mam_comadre_spp <- mam_comadre %>% 
  group_by(gbif.species.id) %>% 
  summarise(gbif.species = gbif.species[1],
            generation_time = mean(generation_time, na.rm = T),
            life_expectancy = mean(life_expectancy, na.rm = T),
            adult_survival = mean(adult_survival, na.rm = T)) %>% 
  ungroup()

# coefficient data merge
mam_coef <- mnanom_5km_GAM %>% 
  ungroup() %>%  ## <- The groups ruin the z-transformations
  left_join(x = ., y = dplyr::select(lpd_gbif, Binomial, gbif.species.id), 
            by = "Binomial") %>% 
  left_join(x = ., y = dplyr::select(mam_comadre_spp, -2),
            by = "gbif.species.id") %>% 
  mutate(species = gsub(" ", "_", gbif_species),
         phylo = species,
         # absolute values of z transformed coefficients
         abs_temp = abs(as.numeric(scale(coef_temp))),
         abs_precip = abs(as.numeric(scale(coef_precip))),
         # z transform absolute latitude for models
         lat = as.numeric(scale(abs(Latitude))),
         # observation-level term for residual variance (not sure if needed)
         OLRE = 1:n(),
         sample_size = as.numeric(scale(log(n_obs)))) %>% 
  dplyr::select(id = ID, id_block = ID_block, n_obs, sample_size, 
                order = Order, species, phylo, 
                biome, lat, generation_time, life_expectancy, adult_survival,
                abs_temp, abs_precip) %>% 
  drop_na(generation_time, life_expectancy, adult_survival, abs_precip)

# summary statistics
n_distinct(mam_coef$species)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3. Exploratory plots ####

mam_coef_spp <- mam_coef %>% 
  group_by(species) %>% 
  summarise(mn_temp = mean(abs_temp),
            se_temp = sd(abs_temp)/sqrt(n()),
            mn_precip = mean(abs_precip),
            se_precip = sd(abs_precip)/sqrt(n()),
            n_obs = n(),
            adult_survival = adult_survival[1],
            generation_time = generation_time[1],
            life_expectancy = life_expectancy[1])
#_______________________________________________________________________________
## Temperature
t_as <- ggplot(mam_coef_spp, aes(x = adult_survival, y = mn_temp)) +
  geom_errorbar(aes(ymax = mn_temp + se_temp,
                    ymin = mn_temp - se_temp),
                width = 0.01, colour = temp_colour) +
  geom_point(size = 3, colour = temp_colour) + 
  labs(x = "Adult survival", 
       y = "|Temperature effect on abundance|") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank())

t_gt <- ggplot(mam_coef_spp, aes(x = log10(generation_time), y = mn_temp)) +
  geom_errorbar(aes(ymax = mn_temp + se_temp,
                    ymin = mn_temp - se_temp),
                width = 0.01, colour = temp_colour) +
  geom_point(size = 3, colour = temp_colour) + 
  scale_x_log10() +
  labs(x = expression(paste(log[10], " Generation time")), 
       y = "|Temperature effect on abundance|") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank())

t_le <- ggplot(mam_coef_spp, aes(x = life_expectancy, y = mn_temp)) +
  geom_errorbar(aes(ymax = mn_temp + se_temp,
                    ymin = mn_temp - se_temp),
                width = 0.01, colour = temp_colour) +
  geom_point(size = 3, colour = temp_colour) + 
  scale_x_log10() +
  labs(x = expression(paste(log[10], " Life expectancy")), 
       y = "|Temperature effect on abundance|") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank())
#_______________________________________________________________________________
## Precipitation
p_as <- ggplot(mam_coef_spp, aes(x = adult_survival, y = mn_precip)) +
  geom_errorbar(aes(ymax = mn_precip + se_precip,
                    ymin = mn_precip - se_precip),
                width = 0.01, colour = precip_colour) +
  geom_point(size = 3, colour = precip_colour) + 
  labs(x = "Adult survival", 
       y = "|Precipitation effect on abundance|") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank())

p_gt <- ggplot(mam_coef_spp, aes(x = log10(generation_time), y = mn_precip)) +
  geom_errorbar(aes(ymax = mn_precip + se_precip,
                    ymin = mn_precip - se_precip),
                width = 0.01, colour = precip_colour) +
  geom_point(size = 3, colour = precip_colour) + 
  scale_x_log10() +
  labs(x = expression(paste(log[10], " Generation time")), 
       y = "|Precipitation effect on abundance|") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank())

p_le <- ggplot(mam_coef_spp, aes(x = life_expectancy, y = mn_precip)) +
  geom_errorbar(aes(ymax = mn_precip + se_precip,
                    ymin = mn_precip - se_precip),
                width = 0.01, colour = precip_colour) +
  geom_point(size = 3, colour = precip_colour) + 
  scale_x_log10() +
  labs(x = expression(paste(log[10], " Life expectancy")), 
       y = "|Precipitation effect on abundance|") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank())

(t_as + t_gt + t_le)/
  (p_as + p_gt + p_le)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 4. Link with life history data ####

load("data/lifehistory.RData", verbose = TRUE)

lh_dat_all <-  mnanom_5km_GAM %>% 
  ungroup() %>%  ## <- The groups ruin the z-transformations
  left_join(x = ., y = dplyr::select(lpd_gbif, Binomial, gbif.species.id), 
            by = "Binomial") %>% 
  left_join(x = ., y = dplyr::select(mam_comadre_spp, -2),
            by = "gbif.species.id") %>% 
  left_join(x = ., y = dplyr::select(lifehistory, c(3,7,8)),
            by = "gbif.species.id") %>% 
  dplyr::select(id = ID, gbif_species, gbif.species.id,
                litter, longevity, generation_time, 
                life_expectancy, adult_survival) %>% 
  drop_na(generation_time, life_expectancy, adult_survival, 
          litter, longevity) %>% 
  group_by(gbif_species) %>% 
  summarise_all(.funs = first) %>% 
  mutate_at(.vars = 6:8, .funs = function(x){log(x + 1)}) %>% 
  rename_all(.funs = str_to_sentence) %>% 
  rename(`Generation time` = Generation_time, 
         `Life expectancy` = Life_expectancy,
         `Adult survival` = Adult_survival)

jpeg(filename = "plots/lifehistory_raw/life_history_covariance_finaldata.jpeg",
     width = 15, height = 15, units = "cm", res = 1000)
pairs.panels(lh_dat_all[,c(4:8)], ellipses = FALSE,
             method = "pearson", lm = TRUE)
dev.off()

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 5. Models ####

# Model data
cdat <- mam_coef %>% 
  mutate(gent = as.numeric(scale(generation_time)),
         le = as.numeric(scale(life_expectancy)),
         as = adult_survival) %>% 
  dplyr::select(id, species, sample_size, gent, 
                le, as, abs_temp, abs_precip)

# Base model chain test
plot(density(rgamma(1000,shape = 2, scale = 0.6)))
plot(density(rexp(1000, rate = 10)))


#_______________________________________________________________________________
## Temperature
set.seed(666)
temp_base <- brm(
  abs_temp ~ 1 + sample_size + (1| species),
  data = cdat, family = Gamma(link = "log"),
  prior = c(
    prior(normal(0, 1), class =  Intercept),
    prior(normal(0, 1), class = b, coef = "sample_size"),
    prior(exponential(2), class = sd, group = "species"),
    prior(gamma(2,0.5), class = shape)),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 3, iter = 4000, warmup = 2000
)

temp_as <- brm(
  abs_temp ~ 1 + sample_size + as + (1| species),
  data = cdat, family = Gamma(link = "log"),
  prior = c(
    prior(normal(0, 1), class =  Intercept),
    prior(normal(0, 1), class = b),
    prior(exponential(2), class = sd, group = "species"),
    prior(gamma(2,0.5), class = shape)),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 3, iter = 4000, warmup = 2000
)

temp_gt <- brm(
  abs_temp ~ 1 + sample_size + gent + (1| species),
  data = cdat, family = Gamma(link = "log"),
  prior = c(
    prior(normal(0, 1), class =  Intercept),
    prior(normal(0, 1), class = b),
    prior(exponential(2), class = sd, group = "species"),
    prior(gamma(2,0.5), class = shape)),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 3, iter = 4000, warmup = 2000
)

temp_le <- brm(
  abs_temp ~ 1 + sample_size + le + (1| species),
  data = cdat, family = Gamma(link = "log"),
  prior = c(
    prior(normal(0, 1), class =  Intercept),
    prior(normal(0, 1), class = b),
    prior(exponential(2), class = sd, group = "species"),
    prior(gamma(2,0.5), class = shape)),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 3, iter = 4000, warmup = 2000
)

## Model comparisons
temp_base <- add_criterion(temp_base, criterion = c("loo","waic"))
temp_as <- add_criterion(temp_as, criterion = c("loo","waic"))
temp_gt <- add_criterion(temp_gt, criterion = c("loo","waic"))
temp_le <- add_criterion(temp_le, criterion = c("loo","waic"))

mod_comp_temp <- as.data.frame(loo_compare(temp_base, temp_as, temp_gt, temp_le))

#_______________________________________________________________________________
## Precipitation
set.seed(666)
precip_base <- brm(
  abs_precip ~ 1 + sample_size + (1| species),
  data = cdat, family = Gamma(link = "log"),
  prior = c(
    prior(normal(0, 1), class =  Intercept),
    prior(normal(0, 1), class = b, coef = "sample_size"),
    prior(exponential(2), class = sd, group = "species"),
    prior(gamma(2,0.5), class = shape)),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 3, iter = 4000, warmup = 2000
)

precip_as <- brm(
  abs_precip ~ 1 + sample_size + as + (1| species),
  data = cdat, family = Gamma(link = "log"),
  prior = c(
    prior(normal(0, 1), class =  Intercept),
    prior(normal(0, 1), class = b),
    prior(exponential(2), class = sd, group = "species"),
    prior(gamma(2,0.5), class = shape)),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 3, iter = 4000, warmup = 2000
)

precip_gt <- brm(
  abs_precip ~ 1 + sample_size + gent + (1| species),
  data = cdat, family = Gamma(link = "log"),
  prior = c(
    prior(normal(0, 1), class =  Intercept),
    prior(normal(0, 1), class = b),
    prior(exponential(2), class = sd, group = "species"),
    prior(gamma(2,0.5), class = shape)),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 3, iter = 4000, warmup = 2000
)

precip_le <- brm(
  abs_precip ~ 1 + sample_size + le + (1| species),
  data = cdat, family = Gamma(link = "log"),
  prior = c(
    prior(normal(0, 1), class =  Intercept),
    prior(normal(0, 1), class = b),
    prior(exponential(2), class = sd, group = "species"),
    prior(gamma(2,0.5), class = shape)),
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15),
  chains = 3, iter = 4000, warmup = 2000
)

## Model comparisons
precip_base <- add_criterion(precip_base, criterion = c("loo","waic"))
precip_as <- add_criterion(precip_as, criterion = c("loo","waic"))
precip_gt <- add_criterion(precip_gt, criterion = c("loo","waic"))
precip_le <- add_criterion(precip_le, criterion = c("loo","waic"))

mod_comp_precip <- as.data.frame(loo_compare(precip_base, precip_as, precip_gt, precip_le))

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 6. Model prediction plots ####

#_______________________________________________________________________________
## Data
pred_dat <- tibble(as = seq(0.4,1, length.out = 50),
                   le = seq(-0.9,7.7, length.out = 50),
                   gent = seq(-1.4,5.4, length.out = 50)) %>% 
  expand_grid(sample_size = mean(cdat$sample_size),
              species = unique(cdat$species))

#_______________________________________________________________________________
## Posterior predictions

# Temperature
temp_pred_as <-  brms::posterior_predict(temp_as, newdata = pred_dat, 
                                          type = "response")
temp_pred_gt <-  brms::posterior_predict(temp_gt, newdata = pred_dat, 
                                         type = "response") 
temp_pred_le <-  brms::posterior_predict(temp_le, newdata = pred_dat, 
                                         type = "response")

# Precipitation
precip_pred_as <-  brms::posterior_predict(precip_as, newdata = pred_dat, 
                                         type = "response")
precip_pred_gt <-  brms::posterior_predict(precip_gt, newdata = pred_dat, 
                                         type = "response") 
precip_pred_le <-  brms::posterior_predict(precip_le, newdata = pred_dat, 
                                         type = "response") 

#_______________________________________________________________________________
## Posterior summaries
posterior_summary_lh <- bind_rows(lapply(unique(pred_dat$as), function(x){
  
  cpos = which(pred_dat$as == x) # where is does data have that value, same for all life-history vars
  
  # posterior means - temp
  mn_tas = mean(temp_pred_as[,cpos])
  mn_tgt = mean(temp_pred_gt[,cpos])
  mn_tle = mean(temp_pred_le[,cpos])
  
  # posterior means - temp
  mn_pas = mean(precip_pred_as[,cpos])
  mn_pgt = mean(precip_pred_gt[,cpos])
  mn_ple = mean(precip_pred_le[,cpos])
  
  # prediction intervals - 80% quantiles - temp
  cQuant_tas = quantile(temp_pred_as[,cpos], c(0.1, 0.9))
  cQuant_tgt = quantile(temp_pred_gt[,cpos], c(0.1, 0.9))
  cQuant_tle = quantile(temp_pred_le[,cpos], c(0.1, 0.9))
  
  # prediction intervals - 80% quantiles - precip
  cQuant_pas = quantile(precip_pred_as[,cpos], c(0.1, 0.9))
  cQuant_pgt = quantile(precip_pred_gt[,cpos], c(0.1, 0.9))
  cQuant_ple = quantile(precip_pred_le[,cpos], c(0.1, 0.9))
  
  # return data
  return(tibble(as = x,
                gent = pred_dat[cpos[1],]$gent,
                le = pred_dat[cpos[1],]$le,
                mn_tas = mn_tas, upr_tas = cQuant_tas[2], lwr_tas = cQuant_tas[1],
                mn_tgt = mn_tgt, upr_tgt = cQuant_tgt[2], lwr_tgt = cQuant_tgt[1],
                mn_tle = mn_tle, upr_tle = cQuant_tle[2], lwr_tle = cQuant_tle[1],
                mn_pas = mn_pas, upr_pas = cQuant_pas[2], lwr_pas = cQuant_pas[1],
                mn_pgt = mn_pgt, upr_pgt = cQuant_pgt[2], lwr_pgt = cQuant_pgt[1],
                mn_ple = mn_ple, upr_ple = cQuant_ple[2], lwr_ple = cQuant_ple[1]))
}))


#_______________________________________________________________________________
## Posterior plots

# Species-level data scaling
cdat_spp <- mam_coef_spp %>% 
  mutate(generation_time = as.numeric(scale(generation_time)),
         life_expectancy = as.numeric(scale(life_expectancy)))

#_______________________________________________________________________________
## Temperature
t_as_p <- ggplot(cdat_spp, aes(x = adult_survival, y = mn_temp)) +
  geom_errorbar(aes(ymax = mn_temp + se_temp,
                    ymin = mn_temp - se_temp),
                width = 0.01, colour = temp_colour) +
  geom_point(size = 3, colour = temp_colour) + 
  geom_smooth(data = posterior_summary_lh,
              aes(x = as, y = mn_tas, ymax = upr_tas, ymin = lwr_tas),
              colour = "black", fill = temp_colour, alpha = 0.3,
              stat = "identity") +
  labs(x = "Adult survival", 
       y = "|Temperature effect on abundance|") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank())

t_gt_p <- ggplot(cdat_spp, aes(x = generation_time, y = mn_temp)) +
  geom_errorbar(aes(ymax = mn_temp + se_temp,
                    ymin = mn_temp - se_temp),
                width = 0.01, colour = temp_colour) +
  geom_point(size = 3, colour = temp_colour) + 
  geom_smooth(data = posterior_summary_lh,
              aes(x = gent, y = mn_tgt, ymax = upr_tgt, ymin = lwr_tgt),
              colour = "black", fill = temp_colour, alpha = 0.3,
              stat = "identity") +
  scale_x_continuous(limits = c(-1.5,4)) +
  labs(x = "Generation time (z scored)", 
       y = "|Temperature effect on abundance|") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank())

t_le_p <- ggplot(cdat_spp, aes(x = life_expectancy, y = mn_temp)) +
  geom_errorbar(aes(ymax = mn_temp + se_temp,
                    ymin = mn_temp - se_temp),
                width = 0.01, colour = temp_colour) +
  geom_point(size = 3, colour = temp_colour) + 
  geom_smooth(data = posterior_summary_lh,
              aes(x = le, y = mn_tle, ymax = upr_tle, ymin = lwr_tle),
              colour = "black", fill = temp_colour, alpha = 0.3,
              stat = "identity") +
  scale_x_continuous(limits = c(-1.1,4)) +
  labs(x = "Life expectancy (z scored)", 
       y = "|Temperature effect on abundance|") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank())

#_______________________________________________________________________________
## Precipitation
p_as_p <- ggplot(cdat_spp, aes(x = adult_survival, y = mn_precip)) +
  geom_errorbar(aes(ymax = mn_precip + se_precip,
                    ymin = mn_precip - se_precip),
                width = 0.01, colour = precip_colour) +
  geom_point(size = 3, colour = precip_colour) + 
  geom_smooth(data = posterior_summary_lh,
              aes(x = as, y = mn_tas, ymax = upr_tas, ymin = lwr_tas),
              colour = "black", fill = precip_colour, alpha = 0.3,
              stat = "identity") +
  labs(x = "Adult survival", 
       y = "|Precipitation effect on abundance|") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank())

p_gt_p <- ggplot(cdat_spp, aes(x = generation_time, y = mn_precip)) +
  geom_errorbar(aes(ymax = mn_precip + se_precip,
                    ymin = mn_precip - se_precip),
                width = 0.01, colour = precip_colour) +
  geom_point(size = 3, colour = precip_colour) + 
  geom_smooth(data = posterior_summary_lh,
              aes(x = gent, y = mn_tgt, ymax = upr_tgt, ymin = lwr_tgt),
              colour = "black", fill = precip_colour, alpha = 0.3,
              stat = "identity") +
  scale_x_continuous(limits = c(-1.5,4)) +
  labs(x = "Generation time (z scored)", 
       y = "|Precipitation effect on abundance|") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank())

p_le_p <- ggplot(cdat_spp, aes(x = life_expectancy, y = mn_precip)) +
  geom_errorbar(aes(ymax = mn_precip + se_precip,
                    ymin = mn_precip - se_precip),
                width = 0.01, colour = precip_colour) +
  geom_point(size = 3, colour = precip_colour) + 
  geom_smooth(data = posterior_summary_lh,
              aes(x = le, y = mn_tle, ymax = upr_tle, ymin = lwr_tle),
              colour = "black", fill = precip_colour, alpha = 0.3,
              stat = "identity") +
  scale_x_continuous(limits = c(-1.1,4)) +
  labs(x = "Life expectancy (z scored)", 
       y = "|Precipitation effect on abundance|") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank())


ggsave((t_as_p + t_gt_p + t_le_p)/
       (p_as_p + p_gt_p + p_le_p),
      filename = "plots/manuscript_figures/Supplementary figures/comadre_lh_predictions.jpeg",
      width = 30, height = 20, units = "cm", dpi = 1500)



