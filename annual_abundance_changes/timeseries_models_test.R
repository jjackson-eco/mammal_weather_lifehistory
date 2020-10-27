#####################################################
##                                                 ##
##     Global climate and population dynamics      ##
##                                                 ##
##           Testing timeseries models             ##
##                                                 ##
##                May 26th 2020                    ##
##                                                 ##
#####################################################
rm(list = ls())
options(width = 100)

library(tidyverse)
library(grid)
library(gridExtra)

##__________________________________________________________________________________________________
#### 1. Load data ####

load("../rawdata/mam_IDblocks.RData")
glimpse(mam_IDblocks)

##__________________________________________________________________________________________________
#### 2. Calculating temporal autocorrelation across studies ####

# Applying acf over all studies and returning data
mam_acf <- mam_IDblocks %>% 
  group_by(ID_block) %>% 
  group_modify(~{
    cab = .$ln_abundance
    
    cacf = acf(cab, lag.max = 4, plot = F)
    
    tibble(dplyr::select(., c(3:8, 13, 14))[1,],
           lag = cacf$lag, cor = cacf$acf, 
           ci = qnorm(0.975)/sqrt(length(cab)),
           n_obs = length(cab))
  }) 

ci_median <- qnorm(0.975)/sqrt(14) # confidence limit for the median record with 10 obs

# Having a look across studies
acf_plot <- ggplot(mam_acf, aes(x = lag, y = abs(cor), group = ID_block, colour = n_obs)) +
  geom_hline(yintercept = ci_median, colour = "blue", linetype = "dashed") +
  #geom_hline(yintercept = -ci_median, colour = "blue", linetype = "dashed") +
  geom_hline(yintercept = 0) +
  geom_line(alpha = 0.3)  +
  scale_color_gradient(low = "cadetblue1", high = "midnightblue", 
                       name = "Length of\nrecord\n(years)") +
  labs(x = "Lag length (years)", 
       y = "Absolute correlation coefficient in ln abundance") +
  theme_bw(base_size = 14)

# Number of significant lags across studies for each lag
propsig_plot <- mam_acf %>% 
  group_by(lag) %>% 
  summarise(n.sig = length(which(abs(cor) >= ci_median)),
            prop.sig = n.sig/n()) %>% 
  ggplot(aes(x = lag, y = prop.sig)) +
  geom_col(fill = "lightblue", colour = "black", size = 0.2) +
  labs(x = "Lag", y = "Proportion of significant autocorrelation values") +
  theme_bw(base_size = 13)

ggsave(grid.arrange(acf_plot, propsig_plot, ncol = 1),
       filename = "plots/annual_abundance/timeseries_model_testing/temporal_autocorrelation.jpeg",
       width = 10, height = 11, units = "in", dpi = 400)

##__________________________________________________________________________________________________
#### 3. First order AR model vs. White Noise models across the records and blocks ####

# Fitting AR models for all ID_blocks and returning model selection data
mam_ARmodels <- mam_IDblocks %>% 
  # Have to exclude these studies that have stange problem that prevents ML being fit - exp growth???
  filter(ID_block != "11175_1" & ID_block != "10734_1" & ID_block != "4803_1") %>% 
  group_by(ID_block) %>% 
  group_modify(~{
    #convert abundance to a timeseries
    cts = ts(.$ln_abundance)
    
    # fit autoregressive model and NULL model of white noise
    mod_AR = arima(cts, order = c(1,0,0), method = "ML")
    mod_WN = arima(cts, order = c(0,0,0), method = "ML")
    
    tibble(dplyr::select(., c(3:8, 13, 14))[1,],
           coef_AR = coef(mod_AR)[1],
           AIC_AR = AIC(mod_AR), AIC_WN = AIC(mod_WN),
           n_obs = length(cts))
  }) %>% ungroup() %>% 
  mutate(dAIC = AIC_AR - AIC_WN,
         n_gr = n_obs - (n_obs %% 5))

hist(mam_ARmodels$dAIC)

for(i in unique(mam_IDblocks$ID_block)){
  print(i)
  if(i != "11175_1" & i != "10734_1" & i != "4803_1"){
  cdat = filter(mam_IDblocks, ID_block == i)
  cts = ts(cdat$ln_abundance)
  mod_AR = arima(cts, order = c(1,0,0), method = "ML")}
}

ggplot(mam_ARmodels, aes(y = dAIC, x = factor(n_gr), 
                         colour = factor(n_gr), 
                         fill = factor(n_gr))) +
  geom_hline(yintercept = -2) +
  geom_violin(alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.7) +
  scale_colour_viridis_d(end = 0.8, guide = F,
                         aesthetics = c("colour", "fill")) +
  labs(x = "Years of record (grouped by 5-year bins)", y = "Change in AIC from AR(1) to WN") +
  annotate("text", x = 1, y = -90, hjust = 0,
           label = paste(round(length(which(mam_ARmodels$dAIC <= -2))/nrow(mam_ARmodels) * 100, 1),
                         "% of timeseries records\nhave AIC difference <= -2")) +
  theme_bw(base_size = 12) +
  ggsave("plots/annual_abundance/timeseries_model_testing/AIC_AR.jpeg",
         width = 8, height = 5, units = "in",
         dpi = 400)

# 50% of studies have support from the predictive performance of an AR(1) model
length(which(mam_ARmodels$dAIC <= -2))/nrow(mam_ARmodels)



