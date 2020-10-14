#########################################################
##                                                     ##
##        Global climate and population dynamics       ##
##                                                     ##
##  Mammal Phylogeny and linear weather effects plots  ##
##                                                     ##
##                  August 21st 2020                   ##
##                                                     ##
#########################################################
rm(list = ls())
options(width = 100)

library(tidyverse)
library(ape)
library(phangorn)
library(phytools)
library(caper)
library(ggtree)
library(treeio)
library(cowplot)

##__________________________________________________________________________________________________
#### 1. Loading data ####

## Phylogenetic
load("../Phlogenetics/mam_lpd_pruned.RData", verbose = TRUE)

## Weather effects
load("data/pgr_weather/mnanom_5km.RData", verbose = TRUE)

##__________________________________________________________________________________________________
#### 2. Merging and dropping tips ####

## Getting weather effects for each checked species name
mam_res_names <- mnanom_5km %>% 
  ungroup() %>% 
  dplyr::select(ID, Order, coef_temp, coef_precip) %>% 
  left_join(x = ., 
            y = dplyr::select(LPD_tree_update, ID, checked_speciesname),
            by = "ID") %>% 
  mutate(checked_speciesname = gsub(" ", "_", checked_speciesname)) %>% 
  group_by(checked_speciesname) %>% 
  summarise(Order = Order[1],
            coef_temp = mean(coef_temp), 
            coef_precip = mean(coef_precip),
            .groups = "drop") %>% 
  # restrict to reasonable effect sizes for plotting purposes
  filter(coef_temp <= 2 & coef_temp >= -2 & coef_precip <= 2 & coef_precip >= -2) %>%
  # converting to a factor of effect sizes for plotting
  mutate(coef_temp_fac = cut(coef_temp, breaks = c(-2, seq(-0.5,0.5,length = 7), 2),
                             labels = c("-2 to -0.5", "-0.5 to -0.33", "-0.33 to -0.17", "-0.17 to 0",
                                        "0 to 0.17", "0.17 to 0.33", "0.33 to 0.5", "0.5 to 2")),
         coef_precip_fac = cut(coef_precip, breaks = c(-2, seq(-0.5,0.5,length = 7), 2),
                               labels = c("-2 to -0.5", "-0.5 to -0.33", "-0.33 to -0.17", "-0.17 to 0",
                                          "0 to 0.17", "0.17 to 0.33", "0.33 to 0.5", "0.5 to 2")),
         genus_label = gsub("_.*", "", checked_speciesname))
  
## Dropping tips
mamMCC_weatheff <- keep.tip(mamMCC_pruned, mam_res_names$checked_speciesname)

## Adding in weather effects and Order columns by converting to a tbl_tree and treedata object
tree_plotdat <- as_tibble(mamMCC_weatheff) %>% 
  full_join(x = ., y = mam_res_names, by = c("label" = "checked_speciesname")) %>% 
  as.treedata()

##__________________________________________________________________________________________________
#### 3. Phylogenetic plotting ####

##________________________________________
# Temperature

# inset of the effect sizes
temp_hist <- ggplot(mam_res_names, aes(x = coef_temp_fac, fill = coef_temp_fac)) + 
  geom_bar(colour = "black", show.legend = F, size = 0.1) +
  scale_fill_brewer(palette = "Spectral",direction = -1, na.value = "black") +
  labs(x = "Temperature coefficient", y = "Frequency") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 8, angle = 20, vjust = 0.7),
        panel.background = element_blank())

# the tree
temp_tree <- ggtree(tree_plotdat, aes(colour = coef_temp_fac), layout = "circular", size = 0.6) +
    geom_tippoint(aes(colour = coef_temp_fac), size = 5,show.legend = FALSE, offset = 2) +
    geom_tiplab(aes(label = genus_label), size = 3, colour = "black", offset = 5) +
    scale_color_brewer(palette = "Spectral",direction = -1, na.value = "black", guide = F)

ggdraw() +
  draw_plot(temp_tree) +
  draw_plot(temp_hist, x = 0.4, y = 0.32, width = 0.39, height = 0.18) +
  ggsave("plots/weather_pop_growth/mam_temp_tree.jpeg",
         dpi = 400,height = 11, width = 11, units = "in")

##________________________________________
# Precipitation

# inset of the effect sizes
precip_hist <- ggplot(mam_res_names, aes(x = coef_precip_fac, fill = coef_precip_fac)) + 
  geom_bar(colour = "black", show.legend = F, size = 0.1) +
  scale_fill_brewer(palette = "Spectral",direction = -1, na.value = "black") +
  labs(x = "Precipitation coefficient", y = "Frequency") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 8, angle = 20, vjust = 0.7),
        panel.background = element_blank())

# the tree
precip_tree <- ggtree(tree_plotdat, aes(colour = coef_precip_fac), layout = "circular", size = 0.6) +
  geom_tippoint(aes(colour = coef_precip_fac), size = 5,show.legend = FALSE, offset = 2) +
  geom_tiplab(aes(label = genus_label), size = 3, colour = "black", offset = 5) +
  scale_color_brewer(palette = "Spectral",direction = -1, na.value = "black", guide = F) 

ggdraw() +
  draw_plot(precip_tree) +
  draw_plot(precip_hist, x = 0.4, y = 0.32, width = 0.39, height = 0.18) +
  ggsave("plots/weather_pop_growth/mam_precip_tree.pdf",
         height = 11, width = 11, units = "in")

rm(list=ls())
options(width = 100)

##__________________________________________________________________________________________________
#### 4. Coefficient plot ####

dat <- tibble(int = -1:1, slp = 1:-1, 
              col = letters[1:3])

ggplot(dat) +
  geom_abline(aes(intercept = int, slope = slp, colour = col), size = 2) +
  coord_cartesian(xlim = c(0,2), ylim = c(-1.5,1.5)) +
  scale_colour_manual(values = c("#d53e4f", "#e6f598", "#3288bd"),
                      labels = c("−2 to −0.5", "−0.17 to 0", "0.5 to 2"),
                      name = "Temperature\ncoefficient") +
  labs(x = "Temperature anomaly", y = "Abundance") +
  theme_bw(base_size = 17) +
  theme(legend.position = c(0.7,0.85),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        legend.background = element_blank()) +
  ggsave(filename = "plots/weather_pop_growth/coefficient_diagram.jpeg", dpi = 400,
         width = 5, height = 5, units = "in")

