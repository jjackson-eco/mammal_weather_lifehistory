####################################################
##                                                ##
##     Global climate and population dynamics     ##
##                                                ##
##     Life-history weather effect exploration    ##
##                                                ##
##                Oct 28th 2020                   ##
##                                                ##
####################################################

## Exploration of potential links between weather effects and life-history

rm(list = ls())
options(width = 100)

library(tidyverse)
library(viridis)
library(patchwork)

##__________________________________________________________________________________________________
#### 1. Load data ####

# mammal coefficient data - Mean anomaly in a 5km buffer radius
load("data/pgr_weather/mnanom_5km.RData")
glimpse(mnanom_5km)

# life-history data
load("data/lifehistory.RData")
glimpse(lifehistory)

# species names 
load("../../LPI_rawdata/GBIF_species_names.RData", verbose = T)
rm(dski_gbif, mamMCC_gbif)

##__________________________________________________________________________________________________
#### 2. Merging data ####

## 2a. Species names for mammal data
mamcoef <- mnanom_5km %>% 
  ungroup() %>% 
  left_join(x = ., y = lpd_gbif, by = "Binomial")

# any not in gbif? NOPE :)
length(which(is.na(mamcoef$gbif.species.id)))


## 2b. Merge with life-history data
coef_lh <- mamcoef %>% 
  left_join(x = ., 
            y = dplyr::select(lifehistory, gbif.species.id, IUCNstatus, 
                              litter, longevity, bodymass), 
            by = "gbif.species.id")
glimpse(coef_lh)


# How many without life history info - It's good news
coef_lh %>% 
  summarise(percent.no.lit = (length(which(is.na(litter))) / n())*100,
            percent.no.lon = (length(which(is.na(longevity))) / n())*100,
            percent.no.bod = (length(which(is.na(bodymass))) / n())*100)

##__________________________________________________________________________________________________
#### 3. Raw weather coefficient life-history plots ####

# set up
tempcol <- viridis(20, option = "plasma")[13]
preccol <- viridis(20, option = "plasma")[4]

mytheme <- theme_bw(base_size = 15) +
  theme(panel.grid = element_blank())

#_________________________________________________
## 3a. temperature
tbm <- ggplot(coef_lh, aes(x = bodymass, y = coef_temp)) +
  geom_point(colour = tempcol, alpha = 0.6, size = 3) +
  labs(x = NULL, y = "Temperature coefficient") +
  mytheme

tln <- ggplot(coef_lh, aes(x = longevity, y = coef_temp)) +
  geom_point(colour = tempcol, alpha = 0.6, size = 3) +
  labs(x = NULL, y = NULL) +
  mytheme

tlt <- ggplot(coef_lh, aes(x = litter, y = coef_temp)) +
  geom_point(colour = tempcol, alpha = 0.6, size = 3) +
  labs(x = NULL, y = NULL) +
  mytheme

tbm + tln + tlt

#_________________________________________________
## 3b. precipitation
pbm <- ggplot(coef_lh, aes(x = bodymass, y = coef_precip)) +
  geom_point(colour = preccol, alpha = 0.6, size = 3) +
  labs(x = "Standardised adult bodymass", y = "Precipitation coefficient") +
  mytheme

pln <- ggplot(coef_lh, aes(x = longevity, y = coef_precip)) +
  geom_point(colour = preccol, alpha = 0.6, size = 3) +
  labs(x = "Standardised maximum longevity", y = NULL) +
  mytheme

plt <- ggplot(coef_lh, aes(x = litter, y = coef_precip)) +
  geom_point(colour = preccol, alpha = 0.6, size = 3) +
  labs(x = "Standardised litter size", y = NULL) +
  mytheme

ggsave((tbm | tln | tlt) /
       (pbm | pln | plt),
       filename = "plots/weather_pop_growth/life_history_weathercoef.jpeg", 
       width = 28, height = 18, units = "cm", dpi = 500)


##__________________________________________________________________________________________________
#### 4. Density dependence and Linear trend + life-history ####

ddcol <- viridis(20, option = "viridis")[13]
trcol <- viridis(20, option = "viridis")[4]

#_________________________________________________
## 4a. Density dependence
ddbm <- ggplot(coef_lh, aes(x = bodymass, y = coef_abun)) +
  geom_point(colour = ddcol, alpha = 0.6, size = 3) +
  labs(x = NULL, y = "Density dependence coefficient") +
  mytheme

ddln <- ggplot(coef_lh, aes(x = longevity, y = coef_abun)) +
  geom_point(colour = ddcol, alpha = 0.6, size = 3) +
  labs(x = NULL, y = NULL) +
  mytheme

ddlt <- ggplot(coef_lh, aes(x = litter, y = coef_abun)) +
  geom_point(colour = ddcol, alpha = 0.6, size = 3) +
  labs(x = NULL, y = NULL) +
  mytheme

#_________________________________________________
## 4b. Linear trend
trbm <- ggplot(coef_lh, aes(x = bodymass, y = coef_trend)) +
  geom_point(colour = trcol, alpha = 0.6, size = 3) +
  labs(x = "Standardised adult bodymass", y = "Linear trend coefficient") +
  mytheme

trln <- ggplot(coef_lh, aes(x = longevity, y = coef_trend)) +
  geom_point(colour = trcol, alpha = 0.6, size = 3) +
  labs(x = "Standardised maximum longevity", y = NULL) +
  mytheme

trlt <- ggplot(coef_lh, aes(x = litter, y = coef_trend)) +
  geom_point(colour = trcol, alpha = 0.6, size = 3) +
  labs(x = "Standardised litter size", y = NULL) +
  mytheme

ggsave((ddbm | ddln | ddlt) /
       (trbm | trln | trlt),
       filename = "plots/weather_pop_growth/life_history_abundance_coef.jpeg", 
       width = 28, height = 18, units = "cm", dpi = 500)

##__________________________________________________________________________________________________
#### 5. Binned weather coefficient life-history plots ####

#_________________________________________________
## 5a. Bodymass
summary(coef_lh$bodymass)
bm_dat <- coef_lh %>% 
  filter(is.na(bodymass) == FALSE) %>% 
  mutate(bm_bin = as.numeric(as.character(cut(bodymass, 
                      breaks = seq(-1.6, 3.2, by = 0.2), 
                      labels = seq(-1.6, 3.0, by = 0.2))))) %>% 
  group_by(bm_bin) %>% 
  summarise(mntemp = mean(coef_temp),
            setemp = sd(coef_temp)/sqrt(n()),
            mnprecip = mean(coef_precip, na.rm = T),
            seprecip = sd(coef_precip, na.rm = T)/sqrt(n()),
            mndd = mean(coef_abun),
            sedd = sd(coef_abun)/sqrt(n()),
            mntr = mean(coef_trend),
            setr = sd(coef_trend)/sqrt(n())) 

bmb_temp <- ggplot(bm_dat, aes(x = bm_bin, y = mntemp)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymax = mntemp + setemp, ymin = mntemp - setemp),
                width = 0.02, colour = tempcol) +
  geom_point(colour = tempcol, size = 3) +
  labs(x = "Standardised adult bodymass", y = "Temperature coefficient") +
  mytheme

bmb_precip <- ggplot(bm_dat, aes(x = bm_bin, y = mnprecip)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymax = mnprecip + seprecip, ymin = mnprecip - seprecip),
                width = 0.02, colour = preccol) +
  geom_point(colour = preccol, size = 3) +
  labs(x = "Standardised adult bodymass", y = "Precipitation coefficient") +
  mytheme

bmb_dd <- ggplot(bm_dat, aes(x = bm_bin, y = mndd)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymax = mndd + sedd, ymin = mndd - sedd),
                width = 0.05, colour = ddcol) +
  geom_point(colour = ddcol, size = 3) +
  labs(x = "Standardised adult bodymass", y = "Density dependence coefficient") +
  mytheme

bmb_tr <- ggplot(bm_dat, aes(x = bm_bin, y = mntr)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymax = mntr + setr, ymin = mntr - setr),
                width = 0.05, colour = trcol) +
  geom_point(colour = trcol, size = 3) +
  labs(x = "Standardised adult bodymass", y = "Linear trend coefficient") +
  mytheme

bmb_temp | bmb_precip | bmb_dd | bmb_tr

#_________________________________________________
## 5b. Longevity
summary(coef_lh$longevity)
ln_dat <- coef_lh %>% 
  filter(is.na(longevity) == FALSE) %>% 
  mutate(ln_bin = as.numeric(as.character(cut(longevity, 
                                              breaks = seq(-2.2, 2.2, by = 0.2), 
                                              labels = seq(-2.2, 2.0, by = 0.2))))) %>% 
  group_by(ln_bin) %>% 
  summarise(mntemp = mean(coef_temp),
            setemp = sd(coef_temp)/sqrt(n()),
            mnprecip = mean(coef_precip, na.rm = T),
            seprecip = sd(coef_precip, na.rm = T)/sqrt(n()),
            mndd = mean(coef_abun),
            sedd = sd(coef_abun)/sqrt(n()),
            mntr = mean(coef_trend),
            setr = sd(coef_trend)/sqrt(n())) 

lnb_temp <- ggplot(ln_dat, aes(x = ln_bin, y = mntemp)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymax = mntemp + setemp, ymin = mntemp - setemp),
                width = 0.02, colour = tempcol) +
  geom_point(colour = tempcol, size = 3) +
  labs(x = "Standardised maximum longevity", y = "Temperature coefficient") +
  mytheme

lnb_precip <- ggplot(ln_dat, aes(x = ln_bin, y = mnprecip)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymax = mnprecip + seprecip, ymin = mnprecip - seprecip),
                width = 0.02, colour = preccol) +
  geom_point(colour = preccol, size = 3) +
  labs(x = "Standardised maximum longevity", y = "Precipitation coefficient") +
  mytheme

lnb_dd <- ggplot(ln_dat, aes(x = ln_bin, y = mndd)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymax = mndd + sedd, ymin = mndd - sedd),
                width = 0.05, colour = ddcol) +
  geom_point(colour = ddcol, size = 3) +
  labs(x = "Standardised maximum longevity", y = "Density dependence coefficient") +
  mytheme

lnb_tr <- ggplot(ln_dat, aes(x = ln_bin, y = mntr)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymax = mntr + setr, ymin = mntr - setr),
                width = 0.05, colour = trcol) +
  geom_point(colour = trcol, size = 3) +
  labs(x = "Standardised maximum longevity", y = "Linear trend coefficient") +
  mytheme

lnb_temp | lnb_precip | lnb_dd | lnb_tr

#_________________________________________________
## 5c. Litter size
summary(coef_lh$litter)
lt_dat <- coef_lh %>% 
  filter(is.na(litter) == FALSE) %>% 
  mutate(lt_bin = as.numeric(as.character(cut(litter, 
                                              breaks = seq(-1.2, 2.8, by = 0.2), 
                                              labels = seq(-1.2, 2.6, by = 0.2))))) %>% 
  group_by(lt_bin) %>% 
  summarise(mntemp = mean(coef_temp),
            setemp = sd(coef_temp)/sqrt(n()),
            mnprecip = mean(coef_precip, na.rm = T),
            seprecip = sd(coef_precip, na.rm = T)/sqrt(n()),
            mndd = mean(coef_abun),
            sedd = sd(coef_abun)/sqrt(n()),
            mntr = mean(coef_trend),
            setr = sd(coef_trend)/sqrt(n())) 

ltb_temp <- ggplot(lt_dat, aes(x = lt_bin, y = mntemp)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymax = mntemp + setemp, ymin = mntemp - setemp),
                width = 0.02, colour = tempcol) +
  geom_point(colour = tempcol, size = 3) +
  labs(x = "Standardised litter size", y = "Temperature coefficient") +
  mytheme

ltb_precip <- ggplot(lt_dat, aes(x = lt_bin, y = mnprecip)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymax = mnprecip + seprecip, ymin = mnprecip - seprecip),
                width = 0.02, colour = preccol) +
  geom_point(colour = preccol, size = 3) +
  labs(x = "Standardised litter size", y = "Precipitation coefficient") +
  mytheme

ltb_dd <- ggplot(lt_dat, aes(x = lt_bin, y = mndd)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymax = mndd + sedd, ymin = mndd - sedd),
                width = 0.05, colour = ddcol) +
  geom_point(colour = ddcol, size = 3) +
  labs(x = "Standardised litter size", y = "Density dependence coefficient") +
  mytheme

ltb_tr <- ggplot(lt_dat, aes(x = lt_bin, y = mntr)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymax = mntr + setr, ymin = mntr - setr),
                width = 0.05, colour = trcol) +
  geom_point(colour = trcol, size = 3) +
  labs(x = "Standardised litter size", y = "Linear trend coefficient") +
  mytheme

ltb_temp | ltb_precip | ltb_dd | ltb_tr

#_________________________________________________
## 5d. Big plot

ggsave(
  (bmb_temp | bmb_precip | bmb_dd | bmb_tr) /
  (lnb_temp | lnb_precip | lnb_dd | lnb_tr) /
  (ltb_temp | ltb_precip | ltb_dd | ltb_tr),
  filename = "plots/weather_pop_growth/life_history_coefficients_binned.jpeg",
  width = 39, height = 35, units = "cm", dpi = 600
)






