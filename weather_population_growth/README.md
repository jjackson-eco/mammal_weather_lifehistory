# Weather effects on population growth rates for the terrestrial mammals

#### 2020-06-09
#### John Jackson

---

This markdown is intended as an accompaniment to the scripts contained within the directory `weather_population_growth`, to walk through the process of calculating the effect of annual weather on population growth rates for each timeseries record in the LPD for the terrestrial mammals. This serves as the first step in a two-step meta regression approach to explore global weather effects on population growth. Please refer to the scripts mentioned in each section of the markdown for full details on each section.

There are 2 main sections and scripts:

## 1. Annual mean weather anomalies
<details>
  <summary>Click here to expand</summary>

### `annual_mean_anomaly.R`

First, we will walk through the process for calculating weather effects using the mean annual weather anomaly for a 5km buffer radius around the study site to demonstrate the process before expanding this out to look across different radius sizes and for different weather variables. We need to join the annual chelsa anomaly data with our population growth data first:

```
##__________________________________________________________________________________________________
#### 1. Load data ####

# mammal data
load("../rawdata/mammal.RData")
glimpse(mammal)

# annual weather anomaly - focus on just the mean anomaly in this script at a 5km range
mam_chelsa_annual <- readRDS("data/mam_chelsa_annual.RDS") %>% 
  filter(scale == "scale_5km") %>% 
  dplyr::select(ID,year, weather_scale = scale, mean_temp_anomaly, mean_precip_anomaly)
glimpse(mam_chelsa_annual)

##__________________________________________________________________________________________________
#### 2. Joining data ####

mammal_weather <- mammal %>% 
  left_join(., y = mam_chelsa_annual, by = c("ID", "year"))
```

Now using the raw data here, we can begin to explore some of the hypotheses of our study. We want to see if there are consistent weather effects on population growth rates across taxa or spatial scales. We can plot out the population growth rate in each year with respect to the taxanomic Order of the corresponding species and the ecological realm.

<img src="../plots/weather_pop_growth/species_anomal_pgr.jpeg" width="700" />
<img src="../plots/weather_pop_growth/realm_anomal_pgr.jpeg" width="700" />

There doesn't seem to be an obvious pattern across realms or taxa, which is perhaps what we expect.

### Calculating weather effects on population growth rate

Now, we want to look at this hypothesis explictly using the timeseries data from each study, whilst accounting for density dependence and any temporal trends in the data. We estimate weather effects on population growth rate for each record using a linear model. Here, the population growth rate *r* at time *t* is given by

<img src="../plots/weather_pop_growth/model_equation.png" width="500" />

where beta 0 is the intercept, omega (*w*) gives the weather variable at time *t* with coefficient beta 1, X gives the *ln* abundance at time *t* with coefficient beta 2, and y gives the year at time *t* with coefficient beta 3. Thus, this linear model estimates the effect of weather, density dependence and the trend on population growth rates, respectively. In this script, the weather variable is the annual mean temperature and precipitation anomaly at a 5km buffer radius. We calculated linear models and extracted the beta coefficients as follows:

```
pgr_weather <- mammal_weather %>% 
  group_by(ID_block) %>% 
  group_modify(~{
    
    # Temperature
    mod_temp = lm(pop_growth_rate ~ mean_temp_anomaly + ln_abundance + year, data = .)
    
    # Precipitation + dealing with NA values
    if(length(which(is.na(.$mean_precip_anomaly) == T)) == 0){
    mod_precip = lm(pop_growth_rate ~ mean_precip_anomaly + ln_abundance + year, data = .)
    coef_precipmod = mod_precip$coefficients}
    else{coef_precipmod = rep(NA,4)}
    
    tibble(.[1,],
           coef_temp = mod_temp$coefficients[2],
           coef_precip = coef_precipmod[2],
           coef_abun = mod_temp$coefficients[3], 
           coef_trend = mod_temp$coefficients[4],
           
           coef_abun2 = coef_precipmod[3], 
           coef_trend2 = coef_precipmod[4],
           n_obs = nrow(.))
  }) 
```

Note that we extract coefficients for the abundance and the year for both temperature and precipitation. These are largely very similar, but vary based on the differences in the weather effect.

Now we have model coefficients for each of the 502 10> year records for the terrestrial mammals. We can now look at comparative patterns in these coefficients. First, the overall density distributions of each of the coefficients across the records

<img src="../plots/weather_pop_growth/overall_coefficients_mnanom_5km.jpeg" width="700" />

As with the raw data, we can see again from this that there doesn't seem to be a consistent pattern of weather effects on population growth rates for either precipitation or temperature. The same goes for trend effects. This is indicative of population-specific (or other level i.e phylogentic or spatial) responses. 

We can also start to explore the hypotheses of the study by looking at these coefficients across different taxanomic groups, ecological realms and latitude.

<img src="../plots/weather_pop_growth/coef_order_mnanom_5km.jpeg" width="700" />
<img src="../plots/weather_pop_growth/coef_realm_mnanom_5km.jpeg" width="700" />
<img src="../plots/weather_pop_growth/coef_lat_mnanom_5km.jpeg" width="700" />

These coefficients or effects sizes form the basis of our meta-regression approach across taxa.

</details>

## 2. Other annual weather variables and scales
<details>
  <summary>Click here to expand</summary>

### `annual_weather_variables.R`

Now we want to repeat the same linear modelling framework but expand to calculate coefficients for all of our annual weather variables and spatial scales. We begin in very much the same way, but don't exclude any of the spatial scales or weather variables.

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

To estimate weather effects for each record, we iterate through weather variables and spatial scales for each, fit a linear model that also incorporates density dependence and trend effects (as above), and extract the weather effects. 

```
##__________________________________________________________________________________________________
#### 3. Linear models for each variable and scale for each record ####

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
  else{mod_weather = lm(pop_growth_rate ~ weather_val + ln_abundance + year, data = cdat)
       modcoef = coefficients(mod_weather)}
  
  # returning data
  cat('\r',"Your Job is",round((x/nrow(iter_dat))*100, 0),"% Complete       ")
  return(tibble(crow, coef_weather = modcoef[2], 
                coef_abun = modcoef[3], coef_trend = modcoef[4],
                rec_info))
}))
  
# 3c. Adding in weather variable labels
pgr_weather_res <- pgr_weather_res %>% 
  mutate(weather_var_lab = stringr::str_to_sentence(gsub("_", " ", weather_var))) %>% 
  mutate(weather_var_lab = gsub("emp", "emperature", weather_var_lab),
         weather_var_lab = gsub("recip", "recipitation", weather_var_lab))
```

This gives us weather coefficients for each variable and scale of our 502 records. Assuming first that all spatial scales are ~identical in their effect size, here we plot the density of the weather coefficient for each of the weather variables.

<img src="../plots/weather_pop_growth/coef_weather_vars.jpeg" width="800" />

</details>