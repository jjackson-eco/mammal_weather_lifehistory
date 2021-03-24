####################################################
##                                                ##
##  Global climate and mammal population dynamics ##
##                                                ##
##              Manuscript figures                ##
##                                                ##
##                Jan 27th 2021                   ##
##                                                ##
####################################################

# GLobal map with studies -

rm(list = ls())
options(width = 100)

## Basic
library(tidyverse)
library(patchwork)
library(ggdist)
library(viridis)
library(cowplot)

## Bayesian
library(brms)
library(tidybayes)
library(ggridges)

## Phylogenetic
library(ape)
library(phangorn)
library(phytools)
library(caper)
library(ggtree)
library(treeio)

# Maps
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)
library(sf)
library(sp)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 1. Load data ####

#_______________
## Analysis data
load("data/mammal_analysis_data_GAM.RData", verbose = TRUE)

#_____________
## Best models
# Gaussian raw coefficients
load("results/gaussian_models/temp_biome_rawcoef.RData")
load("results/gaussian_models/precip_biome_rawcoef.RData")

# Gamma absolute coefficients
temp_lh_uni <- readRDS("results/UCloud_gamma_models/temp_lh_uni_gamma.rds")
precip_lh_uni <- readRDS("results/UCloud_gamma_models/precip_lh_uni_gamma.RDS")

#______________________
## Model selection data
load("results/gaussian_models/model_comparison_temp_rawcoef.RData")
load("results/gaussian_models/model_comparison_precip_rawcoef.RData")

load("results/UCloud_gamma_models/temperature_model_comparisons.RData")
temp_modcomp <- mod_comp
load("results/UCloud_gamma_models/precipitation_model_comparisons.RData")
precip_modcomp <- mod_comp
rm(mod_comp)

#______________________
## Posterior summaries from the gamma model
load("results/UCloud_gamma_models/posterior_summary_Gamma.RData", verbose = TRUE)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### Figure 1 - World map with studies ####

# Colour
point_colour <- viridis(20, option = "C")[14]

# The CRS - Robin projection for curved edges
myCRS <- CRS("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

# World map data
world_sf <- ne_countries(scale = "medium", returnclass = "sf")

# mammal coefficient spatial data
mam_coord <- mam_coef %>% 
  dplyr::select(id, Longitude, Latitude, n_obs, coef_temp, coef_precip) %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), 
           crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

map_plot <- ggplot(data = world_sf) +
  geom_sf(size = 0, fill = "lightgrey") + 
  geom_sf(data = mam_coord, aes(size = n_obs), 
          colour = point_colour, alpha = 0.3) +
  scale_size_continuous(range = c(0.1,3)) +
  coord_sf(crs = myCRS) +
  guides(size = guide_legend(title = "Record length (years)", 
                             override.aes = list(alpha = 1))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position = "top")

record_length <- mam_coef %>% 
  mutate(record_length = n_obs + 1) %>% 
  ggplot(aes(x = record_length)) +
  geom_histogram(bins = 13, fill = point_colour, colour = "black", size = 0.05) +
  scale_x_continuous(breaks = seq(10,40, by = 5)) +
  labs(x = "Record length (years)", y = "Number of records") +
  theme_bw(base_size = 6) +
  theme(panel.grid = element_blank())

order_freq <- ggplot(mam_coef, aes(x = order)) +
  geom_bar(fill = point_colour, colour = "black", size = 0.05) +
  labs(x = "", y = "Number of records") +
  coord_flip() +
  theme_bw(base_size = 6) +
  theme(panel.grid = element_blank())
  
# the final figure
figS1 <- map_plot / (record_length | order_freq ) +
  plot_layout(heights = c(4, 1))

ggsave(figS1, filename = "plots/manuscript_figures/figureS1.jpeg", dpi = 700,
       width = 18, height = 13, units = "cm")


##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### Figure 2 - Raw coefficient plots ####

##______________________________________________________________________________
## The mammal tree
mamMCC_coef$tip.label <- gsub("_", " ", mamMCC_coef$tip.label)

mamtree <- ggtree(mamMCC_coef, size = 0.4,) +
  geom_tiplab(size = 2, hjust = 1) +
  theme_tree() +
  scale_x_reverse(limits = c(200,0))

##______________________________________________________________________________
## The heatmaps
hmdat <- mam_coef %>% 
  mutate(label = gsub("_", " ", phylo)) %>% 
  group_by(phylo) %>% 
  mutate(record = 1:n()) %>% 
  ungroup() %>% 
  left_join(x = ., y = dplyr::select(mamtree$data, label, pos = y), ## Pull out the positions from the plotting tree data.frame
            by = "label") %>% 
  dplyr::select(phylo, pos, record, coef_temp, coef_precip) 

temp_hm <- ggplot(hmdat, aes(x = record, y = pos, fill = coef_temp)) + 
  geom_tile() +
  scale_fill_viridis_c(option = "C", name = "Temperature effect") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(breaks = 1:max(hmdat$pos), expand = c(0,0),
                     sec.axis = dup_axis()) +
  labs(x = "Record number", y = NULL) +
  guides(fill = guide_colorbar(barheight = 1.6, barwidth = 15, 
                               direction = "horizontal", 
                               label.position = "top", 
                               title.position = "top",
                               title.hjust = 0.5)) +
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "top")

precip_hm <- ggplot(hmdat, aes(x = record, y = pos, fill = coef_precip)) + 
  geom_tile() +
  scale_fill_viridis_c(option = "D", name = "Precipitation effect") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(breaks = 1:max(hmdat$pos), expand = c(0,0),
                     sec.axis = dup_axis()) +
  labs(x = "Record number", y = NULL) +
  guides(fill = guide_colorbar(barheight = 1.6, barwidth = 15, 
                               direction = "horizontal", 
                               label.position = "top", 
                               title.position = "top",
                               title.hjust = 0.5)) +
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "top")

ggsave(temp_hm + precip_hm + mamtree + plot_layout(widths = c(2,2,3)),
       filename = "plots/manuscript_figures/phylo_weather_effects.jpeg",
       width = 35, height = 33, units = "cm", dpi = 1500)

##______________________________________________________________________________
### Extracting the colours in photoshops colour picker
temp_colour <- "#990a80"
precip_colour <- "#287f79"

## The posterior plots
temp_intercept_posterior <- temp_biome %>%
  gather_draws(`b_Intercept`) %>% #tidybayes
  ungroup() %>% 
  ggplot(aes(y = .variable, x = .value)) + 
  stat_halfeye(fill = temp_colour) +
  geom_vline(xintercept = 0, size = 0.8) +
  scale_x_continuous(breaks = seq(-0.75,0.75,by = 0.25)) +
  scale_y_discrete(labels = expression(paste(bar(alpha)))) +
  labs(x = "Posterior estimate", y = NULL, tag = "a)") +
  theme_ridges(center_axis_labels = TRUE, grid = T, line_size = 0.3) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

precip_intercept_posterior <- precip_biome %>%
  gather_draws(`b_Intercept`) %>% #tidybayes
  ungroup() %>% 
  ggplot(aes(y = .variable, x = .value)) + 
  stat_halfeye(fill = precip_colour) +
  geom_vline(xintercept = 0, size = 0.8) +
  scale_x_continuous(breaks = seq(-0.75,0.75,by = 0.25)) +
  scale_y_discrete(labels = expression(paste(bar(alpha)))) +
  labs(x = "Posterior estimate", y = NULL, tag = "b)") +
  theme_ridges(center_axis_labels = TRUE, grid = T, line_size = 0.3) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

##______________________________________________________________________________
## Extract variance posteriors from the gaussian model
temp_sdpost <- temp_biome %>%
  gather_draws(`sd_.*`, regex = TRUE) %>% 
  ungroup() %>% 
  mutate(weathervar = "Temperature")

precip_sdpost <- precip_biome %>%
  gather_draws(`sd_.*`, regex = TRUE) %>% 
  ungroup() %>% 
  mutate(weathervar = "Precipitation")

## The plot
sdspp_plot <- bind_rows(temp_sdpost, precip_sdpost) %>% 
  mutate(weathervar = factor(weathervar, levels = c("Temperature", "Precipitation"))) %>% 
  filter(.variable == "sd_species__Intercept") %>% 
  ggplot(aes(y = .variable, x = .value, fill = weathervar)) + 
  stat_halfeye(size = 1.5,) +
  scale_fill_manual(values = c(temp_colour, precip_colour), guide = F) +
  scale_y_discrete(labels = expression(paste("Within species variance ", sigma[SPECIES]))) +
  scale_x_continuous(breaks = seq(0,0.6, by = 0.2), limits = c(0,0.6)) +
  facet_wrap(~ weathervar) +
  labs(x = NULL, y = NULL, tag = "c)") +
  theme_ridges(center_axis_labels = TRUE, grid = T, line_size = 0.3) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.background = element_blank())

sdphylo_plot <- bind_rows(temp_sdpost, precip_sdpost) %>% 
  mutate(weathervar = factor(weathervar, levels = c("Temperature", "Precipitation"))) %>% 
  filter(.variable == "sd_phylo__Intercept") %>% 
  ggplot(aes(y = .variable, x = .value, fill = weathervar)) + 
  stat_halfeye(size = 1.5,) +
  scale_fill_manual(values = c(temp_colour, precip_colour), guide = F) +
  scale_y_discrete(labels = expression(paste("Phylogenetic covariance ", sigma[PHYLO]^2))) +
  scale_x_continuous(breaks = seq(0,0.6, by = 0.2), limits = c(0,0.6)) +
  facet_wrap(~ weathervar) +
  labs(x = "Posterior estimate", y = NULL) +
  theme_ridges(center_axis_labels = TRUE, grid = T, line_size = 0.3) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.background = element_blank(), 
        strip.text = element_blank())

ggsave(sdspp_plot/sdphylo_plot, filename = "plots/manuscript_figures/figure1_inset.jpeg",
       width = 18.4, height = 8.5, units = "cm", dpi = 1500)

layout_fig2 <- "AABB###
                DDEEFFF
                DDEEFFF
                DDEEFFF
                DDEEFFF
                DDEEFFF
                DDEEFFF
                DDEEFFF
                DDEEFFF
                DDEEFFF
                DDEEFFF"

fig1 <- wrap_plots(A = temp_intercept_posterior, 
           B = precip_intercept_posterior,
           D = temp_hm, 
           E = precip_hm,
           `F` = mamtree,
           design = layout_fig2) 

ggsave(fig1, filename = "plots/manuscript_figures/figure1.jpeg",
       width = 35, height = 36.3, units = "cm", dpi = 1500)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### Figure 3 - Posteriors predictions of life-history ####

#_______________________________________________________________________________
## Posterior half eye

# posteriors
temp_lh_post <- temp_lh_uni %>%
  gather_draws(`b_l.*`, regex = TRUE) %>% 
  ungroup() 

precip_lh_post <- precip_lh_uni %>%
  gather_draws(`b_l.*`, regex = TRUE) %>% 
  ungroup() 

# plots
post_temp_lon <- temp_lh_post %>% 
  filter(.variable == "b_longevity") %>% 
  ggplot(aes(x = .value, y = .variable)) +
  stat_halfeye(fill = temp_colour, size = 2) +
  geom_vline(xintercept = 0, size = 0.8) + 
  scale_y_discrete(labels = expression(paste(beta[LONGEVITY]))) +
  scale_x_continuous(breaks = seq(-0.75, 0.5, by = 0.25)) +
  labs(x = "Posterior estimate", y = NULL) +
  theme_ridges(center_axis_labels = TRUE, grid = T, line_size = 0.3) +
  theme(axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        strip.background = element_blank(), 
        strip.text = element_blank())

post_temp_lit <- temp_lh_post %>% 
  filter(.variable == "b_litter") %>% 
  ggplot(aes(x = .value, y = .variable)) +
  stat_halfeye(fill = temp_colour, size = 2) +
  geom_vline(xintercept = 0, size = 0.8) + 
  scale_y_discrete(labels = expression(paste(beta[LITTER]))) +
  scale_x_continuous(breaks = seq(-0.75, 0.75, by = 0.25)) +
  labs(x = "Posterior estimate", y = NULL) +
  theme_ridges(center_axis_labels = TRUE, grid = T, line_size = 0.3) +
  theme(axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        strip.background = element_blank(), 
        strip.text = element_blank())

post_precip_lon <- precip_lh_post %>% 
  filter(.variable == "b_longevity") %>% 
  ggplot(aes(x = .value, y = .variable)) +
  stat_halfeye(fill = precip_colour, size = 2) +
  geom_vline(xintercept = 0, size = 0.8) + 
  scale_y_discrete(labels = expression(paste(beta[LONGEVITY]))) +
  scale_x_continuous(breaks = seq(-0.75, 0.5, by = 0.25)) +
  labs(x = "Posterior estimate", y = NULL) +
  theme_ridges(center_axis_labels = TRUE, grid = T, line_size = 0.3) +
  theme(axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        strip.background = element_blank(), 
        strip.text = element_blank())

post_precip_lit <- precip_lh_post %>% 
  filter(.variable == "b_litter") %>% 
  ggplot(aes(x = .value, y = .variable)) +
  stat_halfeye(fill = precip_colour, size = 2) +
  geom_vline(xintercept = 0, size = 0.8) + 
  scale_y_discrete(labels = expression(paste(beta[LITTER]))) +
  scale_x_continuous(breaks = seq(-0.75, 0.75, by = 0.25)) +
  labs(x = "Posterior estimate", y = NULL) +
  theme_ridges(center_axis_labels = TRUE, grid = T, line_size = 0.3) +
  theme(axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        strip.background = element_blank(), 
        strip.text = element_blank())

#_______________________________________________________________________________
## Binned data plots

lon_temp_bin <- ggplot(longevity_binned_data, aes(x = lon_bin, y = mntemp)) +
  geom_errorbar(aes(ymax = mntemp + setemp, ymin = mntemp - setemp),
                width = 0.03, colour = temp_colour) +
  geom_point(size = 4, colour = temp_colour) +
  geom_smooth(data = posterior_summary_longevity,
              aes(x = longevity, y = post_mn_temp, 
                  ymax = uprPI_temp, ymin = lwrPI_temp),
              colour = "black", fill = temp_colour, alpha = 0.3,
              stat = "identity") +
  coord_cartesian(ylim = c(0,3)) +
  scale_size_continuous(range = c(2,8), guide = F) +
  labs(x = "Standardised longevity", y = "|Temperature effect on abundance|",
       tag = "a)") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank()) 

lit_temp_bin <- ggplot(litter_binned_data, aes(x = lit_bin, y = mntemp)) +
  geom_errorbar(aes(ymax = mntemp + setemp, ymin = mntemp - setemp),
                width = 0.03, colour = temp_colour) +
  geom_point(size = 4, colour = temp_colour) +
  geom_smooth(data = posterior_summary_litter,
              aes(x = litter, y = post_mn_temp, 
                  ymax = uprPI_temp, ymin = lwrPI_temp),
              colour = "black", fill = temp_colour, alpha = 0.3,
              stat = "identity") +
  scale_size_continuous(range = c(2,8), guide = F) +
  labs(x = "Standardised litter size", y = "|Temperature effect on abundance|",
       tag = "b)") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank()) 

lon_precip_bin <- ggplot(longevity_binned_data, aes(x = lon_bin, y = mnprecip)) +
  geom_errorbar(aes(ymax = mnprecip + seprecip, ymin = mnprecip - seprecip),
                width = 0.03, colour = precip_colour) +
  geom_point(size = 4, colour = precip_colour) +
  geom_smooth(data = posterior_summary_longevity,
              aes(x = longevity, y = post_mn_precip, 
                  ymax = uprPI_precip, ymin = lwrPI_precip),
              colour = "black", fill = precip_colour, alpha = 0.3,
              stat = "identity") +
  coord_cartesian(ylim = c(0,4)) +
  scale_size_continuous(range = c(2,8), guide = F) +
  labs(x = "Standardised longevity", y = "|Precipitation effect on abundance|",
       tag = "c)") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank()) 

lit_precip_bin <- ggplot(litter_binned_data, aes(x = lit_bin, y = mnprecip)) +
  geom_errorbar(aes(ymax = mnprecip + seprecip, ymin = mnprecip - seprecip),
                width = 0.03, colour = precip_colour) +
  geom_point(size = 4, colour = precip_colour) +
  geom_smooth(data = posterior_summary_litter,
              aes(x = litter, y = post_mn_precip, 
                  ymax = uprPI_precip, ymin = lwrPI_precip),
              colour = "black", fill = precip_colour, alpha = 0.3,
              stat = "identity") +
  coord_cartesian(ylim = c(0,3.5)) +
  scale_size_continuous(range = c(2,8), guide = F) +
  labs(x = "Standardised litter size", y = "|Precipitation effect on abundance|",
       tag = "d)") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank()) 


#_______________________________________________________________________________
## the figure
fig2a <- lon_temp_bin +
  inset_element(post_temp_lon, 0.3, 0.7, 1, 1, align_to = "panel")

fig2b <- lit_temp_bin +
  inset_element(post_temp_lit, 0, 0.7, 0.6, 1, align_to = "panel")

fig2c <- lon_precip_bin +
  inset_element(post_precip_lon, 0.3, 0.7, 1, 1, align_to = "panel")

fig2d <- lit_precip_bin +
  inset_element(post_precip_lit, 0, 0.7, 0.6, 1, align_to = "panel")

# the final figure
(fig2a | fig2b)/
  (fig2c | fig2d)

ggsave((fig2a | fig2b)/
         (fig2c | fig2d),
       file = "plots/manuscript_figures/figure2.jpeg",
       width = 20, height = 20, units = "cm", dpi = 600)


##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### Results text ####

# raw data
load("../rawdata/mam.RData", verbose = T)

# number of observations
summarise(mam_coef, mn_n = mean(n_obs + 1), md_n = median(n_obs + 1))

# the number for each order
table(mam_coef$order)

# records per species
mam_coef %>% 
  group_by(species) %>% 
  summarise(n = n()) %>% 
  summarise(min_n = min(n), max_n = max(n),
            mean_n = mean(n), md_n = median(n))


# coefficient responses
mean(mam_coef$coef_temp_raw)
sd(mam_coef$coef_temp_raw)

mean(mam_coef$coef_precip_raw, na.rm = T)
sd(mam_coef$coef_precip_raw, na.rm = T)

quantile(mam_coef$coef_temp_raw, probs = c(0.025, 0.975))
quantile(na.omit(mam_coef$coef_precip_raw), probs = c(0.025, 0.975))
    
# extreme records
filter(mam_coef, coef_temp_raw > 1 | coef_temp_raw < -1)
filter(mam_coef, coef_precip_raw > 1 | coef_precip_raw < -1)

# gaussian model
temp_biome
precip_biome

filter(mam_coef, species == "Myodes_glareolus") %>% 
  summarise(min_t = min(coef_temp_raw), max_t = max(coef_temp_raw),
            min_p = min(coef_precip_raw, na.rm = TRUE), 
            max_p = max(coef_precip_raw, na.rm = TRUE))

filter(mam_coef, order == "Rodentia") %>% 
  summarise(min_t = min(coef_temp_raw), max_t = max(coef_temp_raw),
            min_p = min(coef_precip_raw, na.rm = TRUE), 
            max_p = max(coef_precip_raw, na.rm = TRUE))


# gamma
temp_lh_uni
precip_lh_uni

temp_modcomp
precip_modcomp

# raw life-history
load("../rawdata/mam_dski.RData")

# minimum and maximum longevity
filter(mam_coef, longevity == min(longevity))$species
filter(mam_coef, longevity == max(longevity))$species

# 80 years
filter(mam_dski, gbif.species == "Loxodonta africana") %>% 
  as.data.frame()

# 10 months
filter(mam_dski, gbif.species == "Akodon azarae") %>% 
  as.data.frame()

summarise(posterior_summary_longevity, min_t = min(post_mn_temp),
          max_t = max(post_mn_temp), min_p = min(post_mn_precip),
          max_p = max(post_mn_precip),
          tfold = max_t/min_t, pfold = max_p/min_p)

# minimum and maximum litter size
filter(mam_coef, litter== min(litter))$species
filter(mam_coef, litter == max(litter))$species

# 1 baby
filter(mam_dski, gbif.species == "Rhinolophus ferrumequinum") %>% 
  as.data.frame()

# 13 babys
filter(mam_dski, gbif.species == "Thylamys elegans") %>% 
  as.data.frame()

summarise(posterior_summary_litter, min_t = min(post_mn_temp),
          max_t = max(post_mn_temp), min_p = min(post_mn_precip),
          max_p = max(post_mn_precip),
          tfold = max_t/min_t, pfold = max_p/min_p)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### Figure S14 - Posteriors predictions of sample size effects ####

# prediction data
preddat_n <- expand_grid(longevity = mean(mam_coef$longevity),
                         litter = mean(mam_coef$litter),
                         bodymass = mean(mam_coef$bodymass),
                         sample_size = seq(-1.3,2.6, length.out = 50),
                         phylo = unique(mam_coef$phylo)) %>% 
  mutate(species = phylo)

# posterior predictions
temp_pred_n <-  brms::posterior_predict(temp_lh_uni, newdata = preddat_n, 
                                          type = "response") 
precip_pred_n <-  brms::posterior_predict(precip_lh_uni, 
                                            newdata = preddat_n, 
                                            # account for NAs in the precipitation coefficients
                                            allow_new_levels = TRUE,
                                            type = "response") 



# summary data for each sample size
posterior_summary_n <- bind_rows(lapply(unique(preddat_n$sample_size), function(x){
  
  cpos = which(preddat_n$sample_size == x) # vector of positions corresponding to where the values are equal to x
  
  # posterior mean 
  post_mn_temp = mean(temp_pred_n[,cpos])
  post_mn_precip = mean(precip_pred_n[,cpos])
  
  
  # return data
  return(tibble(sample_size = x, 
                post_mn_temp = post_mn_temp, 
                post_mn_precip = post_mn_precip))
}))

temp_n <- ggplot(mam_coef, aes(x = sample_size, y = abs_temp)) +
  geom_point(size = 3, colour = temp_colour, alpha = 0.3) +
  geom_line(data = posterior_summary_n, aes(y = post_mn_temp),
            colour = "black", size = 1.5) +
  coord_cartesian(ylim = c(0,4)) +
  scale_x_continuous(breaks = seq(min(mam_coef$sample_size),max(mam_coef$sample_size), length = 10),
                     labels = round(seq(9,35, length = 10), 1)) +
  labs(x = "Record length (years)", y = "|Temperature effect on abundance|",
       tag = "a)") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank()) 

precip_n <- ggplot(mam_coef, aes(x = sample_size, y = abs_precip)) +
  geom_point(size = 3, colour = precip_colour, alpha = 0.3) +
  geom_line(data = posterior_summary_n, aes(y = post_mn_precip),
            colour = "black", size = 1.5) +
  coord_cartesian(ylim = c(0,4)) +
  scale_x_continuous(breaks = seq(min(mam_coef$sample_size),max(mam_coef$sample_size), length = 10),
                     labels = round(seq(9,35, length = 10), 1)) +
  labs(x = "Record length (years)", y = "|Precipitation effect on abundance|",
       tag = "b)") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank()) 

ggsave(temp_n + precip_n, 
       filename = "plots/meta_regression/sample_size_posterior_prediction.jpeg",
       width = 25, height = 14, units = "cm", dpi = 500)

