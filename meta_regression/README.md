# Meta regression exploring patterns of weather effects on abundance change in the terrestrial mammals

#### 2021-03-17
#### John Jackson

---

This markdown is intended as an accompaniment to the scripts contained within the directory `meta_regression/`, to walk through the Bayesian meta-regression framework carried out in this study. We carried out two forms of meta-regression model during these analysis:

1. Gaussian meta-regression on weather coefficient values to estimate consistent patterns across the mammals.
2. Gamma meta-regressions on absolute weather coefficients to explore relationships to life-history.


Here we generate the key results and findings presented in the manuscript. Please refer to the scripts mentioned in each section of the markdown for full details on each section.

There are 4 main sections and scripts:

## 1. Prior Predictive Simulation

<details>
  <summary>Click here to expand</summary>

### `testing_and_prior_predictive_simulation/prior_predictive_simulation.R`

Before fitting our Bayesian meta-regression models, we need to develop effective priors. We opted to use conservative, regularising priors following McElreath 2020, which gave estimates within the parameter space observed in the raw data. This was achieved through prior predictive simulation (PPS). Here, we compare the estimates and expectation of priors to the limits of observed data to inform the priors. In addition, priors were further tuned to improve the efficiency/accuracy of Markov chains during model selection analyses.

In all cases, we chose conservative regularising priors to reflect the high number of parameters in phylogenetically or spatially controlled models. 

Here we walk through simulations carried out to inform priors for the global intercept (mean weather coefficient), the beta coefficients for differences (i.e. biome), the life-history slope effects and the random effects variance terms (i.e. phylogenetic covariance and species variance). We present priors of increasing regularisation.

### Intercept terms

Here we used normal priors to describe the intercept of population responses across records. For all priors we used a mean of 0 as we had no prior expectation regarding the direction of weather effects. Then, we used three simulated priors to inform the priors used in the study:

1. Weak prior - Normal(0,10)
2. Medium regularising prior - Normal(0,2)
3. Regularising prior - Normal(0,0.5)

Here we compare the intercept priors to the observed coefficient bounds:

<img src="../plots/meta_regression/prior_predictive_simulation/weather_coefficient_pps.jpeg" width="600" />

### Coefficient difference beta terms

For the beta coefficients describing differences in grouping variables, we looked at pairwise differences in all observed coefficients to inform the parameter space for the prior. Again, we used normal priors and explored the same prior parameters. Here are the simulated differences in coefficients expected by the prior to an intercept of 0.

<img src="../plots/meta_regression/prior_predictive_simulation/coefficient_differences_beta.jpeg" width="600" />

### Life-history effect simulations

For the slope terms describing the effect of life-history on population responses to weather, we used an intercept of 0 once more and then simulated beta slope terms using the same normal priors explored previously. We then predicted weather effects from simulated life-history values between -2 and 2. These plots present the predictions of weather effects for each of the normal priors. The solid and dashed black lines are the maximum observed coefficients for temperature and precipitation, respectively.

<img src="../plots/meta_regression/prior_predictive_simulation/life_history_effect_pps.jpeg" width="800" />

### Random effect variance terms

We used exponential priors when considering variance terms relating to the mixed effects in the meta-regression, which mainly were used for phylogenetic covariance and species variance. Exponential terms were beneficial here as they are non-zero and flexible for exploring large variances. Here, we explored the priors of varying exponential rates from 0.5-20, and their consequent distributions of variance terms. Smaller rates give weaker priors with a wider range of variance terms. In our case, particularly for phylogenetic covariance, we do not expect variance terms > 1. We present the resulting distributions from exponential priors of varying rates. The solid lines indicate a variance of 1.

<img src="../plots/meta_regression/prior_predictive_simulation/random_effect_variance_pps.jpeg" width="600" />

</details>

In all cases, more regularising, conservative priors were much more representative of observed restrictions (i.e. maxima and minima) of the raw data. Furthermore, we only presented isolated priors, without exploring the consequences of adding a greater number of parameters e.g. random effects that would further restrict the coefficients obtained.

Thus, in all subsequent bayesian models, we used 

## 2. Other annual weather variables and scales
<details>
  <summary>Click here to expand</summary>

### `GAM_weather_pop_growth_ALLweathervars.R`

Now we want to repeat the same GAM modelling framework but expand to calculate coefficients for all of our annual weather variables and spatial scales. We begin in very much the same way, but don't exclude any of the spatial scales or weather variables.

```
##__________________________________________________________________________________________________
#### 1. Load data ####

# mammal data
load("../rawdata/mammal.RData")
glimpse(mammal)

# annual weather anomaly - focus on just the mean anomaly in this script at a 5km range
mam_chelsa_annual <- readRDS("data/mam_chelsa_annual.RDS") %>% 
  dplyr::select(-c(4:6))
glimpse(mam_chelsa_annual)

##__________________________________________________________________________________________________
#### 2. Joining data ####

mammal_weather <- mammal %>% 
  left_join(., y = mam_chelsa_annual, by = c("ID", "year"))

```

To estimate weather effects for each record, we iterate through weather variables and spatial scales for each, fit a GAM model that also incorporates trend and temporal autocorrelation (as above), and extract the weather effects. 

```
##__________________________________________________________________________________________________
#### 3. GAMs for each variable and scale for each record ####

# 3a. set up iteration data
# Ignoring number of odd days vars for now - they follow a zero inflated pattern
iter_dat <- expand_grid(ID_block = unique(mammal_weather$ID_block),
                               scale = unique(mammal_weather$scale),
                               weather_var = colnames(mammal_weather)[24:39])

# 3b. weather coefficients for each variable
pgr_weather_res <- bind_rows(lapply(X = 1:nrow(iter_dat), function(x){
  
  crow = iter_dat[x,]
  
  # current data
  cdat = mammal_weather %>% 
    filter(ID_block == crow$ID_block, scale == crow$scale) %>% 
    dplyr::select(ID_block, year, ln_abundance,
                  weather_val = crow$weather_var,
                  pop_growth_rate)
  
  # record info
  rec_info = mammal_weather %>% 
    filter(ID_block == crow$ID_block, scale == crow$scale) %>% 
    dplyr::select(2:17) %>% 
    slice(1)
  
  # model
  if(length(which(is.na(cdat$weather_val) == T)) > 0){modcoef = rep(NA,4)}
  else{mod_weather = gamm(pop_growth_rate ~ 
                            s(year, bs = "tp", k = 5) + weather_val,
                          data = cdat, 
                          family = gaussian,
                          correlation = corARMA(form = ~ year, p = 1),
                          method = "REML")
       modcoef = coef(mod_weather$gam)}
  
  # returning data
  cat('\r',"Your Job is",round((x/nrow(iter_dat))*100, 0),"% Complete       ")
  return(tibble(crow, coef_weather = modcoef[2], 
                rec_info))
}))
  
# 3c. Adding in weather variable labels
pgr_weather_res <- pgr_weather_res %>% 
  mutate(weather_var_lab = stringr::str_to_sentence(gsub("_", " ", weather_var))) %>% 
  mutate(weather_var_lab = gsub("emp", "emperature", weather_var_lab),
         weather_var_lab = gsub("recip", "recipitation", weather_var_lab))
```

This gives us weather coefficients for each variable and scale of our 494 records. Assuming first that all spatial scales are ~identical in their effect size, here we plot the density of the weather coefficient for each of the weather variables. This shows that coefficients of weather effects are largely very small across records. However, there are some cases with large weather coefficients and some distributions that suggest there may be patterns.

<img src="../plots/weather_pop_growth/coef_weather_vars.jpeg" width="600" />

We can also have a look at how the weather coefficients we obtained are different based on the buffer radius or spatial scale that was chosen. Below we can see a pairs.panel plot that displays the correlations in all weather coefficients between the scales. You can see that they are virtually identical.

<img src="../plots/weather_pop_growth/scale_weather_coef.jpeg" width="600" />

### Variance in weather `GAM_weather_variance_pop_growth.R`

An important feature of responses to the environment is that species may be more responsive to variance in the environment rathe than the central tendency i.e. a mean annual temperature anomaly of 0 may not reflect the fact that in reality there were big fluctuations in the monthly weather variables. Therefore, we repeat the process highlighted in step 1 using the 5km buffer radius, but this time using weather variance.

If we refer to `weather_variables/annual_weather_variables.R` in the root github repository (i.e. section 2 of the main workflow), we calculate the weather variance of the raw monthly temperature and precipitation variables such that:

```
# Weather variance
temp_variance = var(temp)
precip_variance = var(precip)
```
Apart from this, the framework for estimating the weather effects on population growth rates is identical using the GAM model. Please see `GAM_weather_variance_pop_growth.R` for the full details.

</details>

## 3. Model comparison
<details>
  <summary>Click here to expand</summary>

### `weather_popgrowth_method_comparison.R`

To test the implications of our choice of GAMs for modelling the underlying patterns of weather effects on population growth rates, we performed a method comparison for different ways of assessing the weather effects. We evaluated a set of five candidate models for estimating the weather effects, which were as follows:

1. A fully naive linear model only including the weather effect.
2. A linear model accounting for the trend in population growth rate.
3. A linear model incorporating trend and the previous years abundance.
4. A glmmTMB model with an autoregressive term for year of order AR(1).
5. A GAM with a corse smoothing term for the year trend and an autoregressive term for year of order AR(1)

These models were specified generally as follows:

```
    #_______________________________________________
    # Temperature
    
    # Model 1. Fully naive model - simple linear regression
    temp_naive = lm(pop_growth_rate ~ mean_temp_anomaly, data = .)
    coef_temp_naive = coef(temp_naive)[2]
    
    # Model 2. Linear regression accounting for trend 
    temp_lintr = lm(pop_growth_rate ~ mean_temp_anomaly + year, data = .)
    coef_temp_lintr = coef(temp_lintr)[2]
    
    # Model 3. Linear regression with trend and past abundance
    temp_linear = lm(pop_growth_rate ~ mean_temp_anomaly + ln_abundance + year, data = .)
    coef_temp_linear = coef(temp_linear)[2]
    
    # Model 4. glmmTMB with autoregression(1) for year
    temp_TMB = glmmTMB(pop_growth_rate ~ mean_temp_anomaly + ar1(as.factor(year_f) + 0 | ID), 
                       family = 'gaussian', data = .)
    coef_temp_TMB = as.numeric(coef(temp_TMB)$cond$ID[length(coef(temp_TMB)$cond$ID)])
    
    # Model 5. GAMM with a coarse year smoothing term and an autoregression (1) correlation structure 
    temp_gamm = gamm(pop_growth_rate ~ s(year, bs = "tp", k = 5) + mean_temp_anomaly,
                    data = ., family = gaussian,
                    correlation = corARMA(form = ~ year, p = 1),
                    method = "REML")
    coef_temp_gamm = coef(temp_gamm$gam)[2]
    
    #_______________________________________________
    # Precipitation
    if(length(which(is.na(.$mean_precip_anomaly) == T)) == 0){
      
      # Model 1. Fully naive model - simple linear regression
      precip_naive = lm(pop_growth_rate ~ mean_precip_anomaly, data = .)
      coef_precip_naive = coef(precip_naive)[2]
      
      # Model 2. Linear regression accounting for trend 
      precip_lintr = lm(pop_growth_rate ~ mean_precip_anomaly + year, data = .)
      coef_precip_lintr = coef(precip_lintr)[2]
      
      # Model 3. Linear regression with trend and past abundance
      precip_linear = lm(pop_growth_rate ~ mean_precip_anomaly + ln_abundance + year, data = .)
      coef_precip_linear = coef(precip_linear)[2]
      
      # Model 4. glmmTMB with autoregression(1) for year
      precip_TMB = glmmTMB(pop_growth_rate ~ mean_precip_anomaly + ar1(as.factor(year_f) + 0 | ID), 
                         family = 'gaussian', data = .)
      coef_precip_TMB = as.numeric(coef(precip_TMB)$cond$ID[length(coef(precip_TMB)$cond$ID)])
      
      # Model 5. GAMM with a coarse year smoothing term and an autoregression (1) correlation structure 
      precip_gamm = gamm(pop_growth_rate ~ s(year, bs = "tp", k = 5) + mean_precip_anomaly,
                       data = ., family = gaussian,
                       correlation = corARMA(form = ~ year, p = 1),
                       method = "REML")
      coef_precip_gamm = coef(precip_gamm$gam)[2]}
```
And we also then explore the relationships between the different methods of calculating weather effects. Here we look at the pairwise correlation plot for all of the coefficients for temperature:

<img src="../plots/weather_pop_growth/temp_effect_comparison.jpeg" width="700" />

And for precipitation:

<img src="../plots/weather_pop_growth/precip_effect_comparison.jpeg" width="700" />

You can see that in particular for the GAM approach taken in this study, it was a good representation of the weather coefficients for both temperature and precipitation relative to other methods. And, as with our density dependence simulation in the previous section, we see that weather coefficients (environmental effects) are highly correlated even when we don't account for temporal autocorrelation or trend. The only exception to these reassuring findings is the glmmTMB method, which is not as well correlated. However these coefficients are still well related to the GAM coefficients.

</details>

## 4. Hypothesis exploration plots
<details>
  <summary>Click here to expand</summary>

### `hypothesis_exploration/`

With a weather effect for each record, we can start to explore the hypotheses of the study by looking at these coefficients across different taxanomic groups, ecological biomes, latitudes, and with respect to life-history variables. Please refer to the scripts in the `hypothesis_exploration/`.

### Spatial variables

Here we look at the distribution of the weather coefficients with respect to the biome and the latitude, both of which are often important in macro-ecological patterns. We would predict that generally, as the climate is more stable in tropical regions, the magnitude population responses to weather at low latitudes and tropical biomes is lower, with more extreme population changes in regions where weather is more changeable. However, exactly because the climate is more stable, we may also expect the opposite.

<img src="../plots/weather_pop_growth/coef_biome_mnanom_5km_GAM.jpeg" width="700" />
<img src="../plots/weather_pop_growth/coef_lat_mnanom_5km_GAM.jpeg" width="700" />

It does certainly look like there are some biomes with more extreme population responses to the weather. Furthermore there seems to be  wider spread of population responses at the most extreme latitudes.

### Evolutionary history

We also could predict that different taxonomic groups, and shared evolutionary history may be responsible for an organisms response to the weather, due to shared adaptation to changes in the environment. Here you have the coefficients distributions for each order of mammals in the study.

<img src="../plots/weather_pop_growth/coef_order_mnanom_5km_GAM.jpeg" width="700" />

We can also look at how these weather coefficients are distributed around the phylogenetic tree for the mammals. Do we see covariance in population responses to weather between closely related species?

<img src="../plots/weather_pop_growth/mam_temp_tree.jpeg" width="800" />
<img src="../plots/weather_pop_growth/mam_precip_tree.jpeg" width="800" />

### Life-history

We also may expect that demographic traits, and traits related to an organisms mode of life, their life-history, has a part to play in the response to changes in the weather. We can look at the temperature and precipitation coefficients with respect to our three key life-history variables maximum longevity, litter size and adult body mass. 

<img src="../plots/weather_pop_growth/life_history_weathercoef.jpeg" width="800" />

These coefficients or effects sizes form the basis of our meta-regression approach across taxa, where we will explore these coefficient patterns in detail.

</details>
