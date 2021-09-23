# Annual abundance data to population growth rates
#### 2021-09-23
#### John Jackson

---

This markdown is intended as an accompaniment to the scripts contained within the directory `annual_abundance_changes`, to walk through the process of managing annual abundance data from the Living Planet Database for the terrestrial mammals, calculating annual population growth rates and exploring how best to account for temporal autocorrelation in abundance. Please refer to the scripts mentioned in each section of the markdown for full details on each section. We initially planned to use a detrending method to account for temporal trends in the abundance time series, but this was replaced with directly accounting for temporal autocorrelation in the models. There are several scripts here that were just exploratory i.e. detrending scripts,

There are 4 main sections and scripts:

## 1. Exploring and accounting for gaps in the abundance timeseries
<details>
  <summary>Click here to expand</summary>

### `timeseries_gap_exploration.R`

For our record question of how weather influences annual population changes in vertebrates, one potential problem with the vertebrate abundance data from the Living Planet Database is that there are gaps in the timeseries. The number of these gaps and the way we deal with them is important. This first section is intended to explore the pervasiveness of these gaps across our records, and deal with them in an appropriate way for further analysis.

We first have to restrict the data to only include records that have sufficient data with which to explore annual changes abundance in relation to CHELSA weather data. We only include observations that overlap with the CHELSA data i.e. between 1979-2013, and those that have 5 or more years (first pass) of abundance data:

```
mam <- mam %>% 
  filter(year >= 1979 & year <= 2013) %>% # important to do this first
  group_by(ID) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  filter(n >= 5) %>% 
  dplyr::select(-n)
```

Now we present the example of a population timeseries with 5 observations of population density from the Fossa (Cryptoprocta ferox), which is endemic to Madagascar. There are 5 observations between 2008 and 2013.

<img src="../plots/annual_abundance/timeseries_gaps/fossa_timeseries.jpeg" width="600"/>

However, we have a gap in the data in 2009. We could interpolate this value when we detrend the timeseries, but since we are investigating changes in annual abundance, we are actually interested in annual deviation in abundance. Therefore, a better strategy is to map these gaps across all of our records and investigate if there is a way of splitting the timeseries up when we investigate the effect of weather.

We make use of the differences in the $years column to split each record's timeseries in to blocks. We compute some summary statistics for each record `mam_gaps`, including the number of these blocks, the proportion of the data that is a timeseries with 1-year transitions, and the longest continuous block in the record. In `mam_blocks`, we record all the blocks for each record and keep the raw data:

```
mam_gaps <- mam %>% 
  group_by(ID) %>% 
  group_modify(~{
    cyears = .$year
    diff_cyears = diff(cyears)
    cumsum_blocks = cumsum(c(1, diff_cyears != 1))
    
    summarise(., Binomial = Binomial[1],
              study_length = length(cyears),
              no_consecutive_blocks = n_distinct(cumsum_blocks),
              prop_1year_transitions = sum(diff_cyears == 1)/ length(diff_cyears),
              longest_block = max(table(cumsum_blocks)))
  }) %>% 
  ungroup()

mam_blocks <- mam %>% 
  group_by(ID) %>%
  mutate(block = cumsum(c(1, diff(year) != 1)),
         max_block = max(block)) %>% 
  ungroup() %>% 
  dplyr::select(ID, Binomial, Order, ln_abundance, 
                year, block, max_block) %>% 
  left_join(x = ., y = dplyr::select(mam_gaps, -c(Binomial, no_consecutive_blocks)),
            by = "ID") %>% 
  arrange(desc(longest_block)) %>% 
  mutate(ID = factor(ID, levels = unique(.$ID)))
```
### Time series gap summaries

Encouragingly, the majority of the data is occurring in continuous blocks with 1-year transitions. Here we have the distribution of the proportion of each record that is occurring in 1-year transitions. We can see that the majority have all their data as 1-year transitions. Furthermore >76% (870) of the records have more than 2 thirds of their data in 1-year transitions.

<img src="../plots/annual_abundance/timeseries_gaps/Proportion_1year_transitions.jpeg" width="700" />

However, this doesn't quite give us the full picture because we also need to know how many blocks each timeseries occurs in. Here we plot the number of blocks that each record occurs in against the number of years in its longest block. The colour denotes the proportion of the timeseries occurring in 1-year transitions.

<img src="../plots/annual_abundance/timeseries_gaps/Consecutive_blocks.jpeg" width="700" />

So it does appear that there are some records that are primarily in timeseries with 1-year transitions (lighter colours), but do occur over quite a few blocks of observations. We can also plot these blocks of observations as timelines, where we see the years of data for each record ID. I have split these up based on the number of blocks that the timeseries occurs in for ease. Here first you have the records that are just occurring in 1 consecutive block. Points and lines indicate where there is data for each record ID (row).

<img src="../plots/annual_abundance/timeseries_gaps/record_timelines/1_consecutive_block_timeline.jpeg" width="700" />

These are the 'gold standard' records that occur solely in one consecutive chain of annual observations (with more than 10 years of data). However, the picture becomes a little bit more complex when we look at records that occur in a greater number of blocks. Here you can see the records that occur in 5 blocks.

<img src="../plots/annual_abundance/timeseries_gaps/record_timelines/5_consecutive_block_timeline.jpeg" width="700" />

We can see here that there are scenarios where there are longer consecutive blocks of observations, with smaller satellite blocks that have fewer observations. Furthermore, we can see in records with more blocks, there are situations where there are several separate blocks of 1-year transitions, but that many blocks have less than 5 observations.

### Data selection

We are selecting data based on the sizes of the blocks for each record - We only want to retain blocks within a record that have 10 or more consecutive annual observations.

```
# IDs and blocks that we want to keep - 502 out of 2756 ID-block combinations
ID_block_keep <- mam_blocks %>% 
  mutate(ID = as.numeric(as.character(ID))) %>% 
  group_by(ID, block) %>% 
  summarise(ID_block = paste0(ID[1],"_",block[1]),
            block_keep = if_else(n() >= 10, 1, 0)) %>%
  ungroup() %>% 
  filter(block_keep == 1)

# Restricting the dataset
mam_IDblocks <- mam %>% 
  group_by(ID) %>%
  mutate(block = cumsum(c(1, diff(year) != 1)),
         ID_block = paste0(ID[1],"_",block)) %>% 
  ungroup() %>% 
  filter(ID_block %in% ID_block_keep$ID_block == T) %>% 
  select(1,21,22,2:9,11:20)

```

This gives us a final dataset with which to assess how weather affects changes in abundance, stored in `mammal`. Each Study has consecutive blocks within the study with at least ten years of data within them. Here is a comparison to the initial raw data:

![](../plots/annual_abundance/data_summary.jpeg)

This equates to 38% of the initial observations, 20% of the initial records, and 31% of the species.

</details>

---

## 2. Investigating temporal autocorrelation
<details>
  <summary>Click here to expand</summary>

### `timeseries_models_test.R`

Because we are working with abundance timeseries, accounting for and exploring the temporal autocorrelation in the timeseries is crucial for the study. We can treat each block in each record explicitly as a timeseries object, investigate whether we have evidence for temporal autocorrelation, and then fit an Autoregressive timeseries model.

First, we apply the `acf` function to each of the record blocks to investigate temporal autocorrelation. We only explored up to a lag of 4 years here, because there was insufficient data for most studies to explore a greater number of lags. The acf function computes the correlation for a vector at all specified lags. We then have to use `qnorm` at the 95% confidence limit with the sample size to give the confidence limit, or level of significance. This is done both for each vector, and as a general confidence limit for the median record length (14 years). Absolute correlation values above this limit are significant.

```
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

ci_median <- qnorm(0.975)/sqrt(14) # confidence limit for the median record with 14 obs
```

We can then plot out the autocorrelation for each of the studies, and the proportion that are significant at each time lag. Here you have the absolute correlation coefficient at each lag for each record, and the dashed line is the significance confidence level for the median study length. Then you have the proportion of the studies at each lag that are significant based on their record length.

<img src="../plots/annual_abundance/timeseries_model_testing/temporal_autocorrelation.jpeg" width="700"/>

So it seems like we have good evidence from a reasonable proportion of records (>25%) that there is temporal autocorrelation at a lag of one year. However, beyond a single year of lag, the absolute correlation and proportion of significant records drops substantially. This is most likely due to the fact we are dealing with annual data, and the restrictive sample size of each record.

Therefore, we conclude it is important to incorporate temporal autocorrelation of lag 1. We can take this further by fitting explicit order 1 autoregressive models, AR(1), to each ~stationary timeseries (have to exclude a couple of studies which seem to exhibit exponential growth). We can then compare the predictive performance (AIC here) of the AR(1) models to a white noise (WN) model, which just assumes random noise with a normal distribution at a specified mean and variance.

```
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
```

The AIC here gives us the predictive performance of the AR(1) and WN models for each record. Lower AIC values indicate improved predictive performance, and we consider a difference of >=2 to indicate improved predictive performance. Here are the results across the studies, grouped by the length of the study in 5-year bins. You have the difference in AIC for each study (points) as well as a density violin across each bin. The black horizontal line is a difference of 2.

<img src="../plots/annual_abundance/timeseries_model_testing/AIC_AR.jpeg" width="700"/>

We can see that just under 50% of the records have statistical support for an AR(1) model compared to white noise. This suggests that it is important to incorporate lag 1 autocorrelation in abundance when exploring annual weather effects.

</details>

---

## 3. Simulating the effect of density dependence on estimating weather effects
<details>
  <summary>Click here to expand</summary>

### `density_dependence_simulation.R`

Following on from the investigation of temporal autocorrelation in abundance measures in the LPD data, a pervasive feature of population dynamics that manifests in temporally autocorrelated abundance is density dependence. Investigating density dependence directly is beyond the scope of this study. Instead however, we wanted to investigate whether the presence of density dependence in the system had the potential to bias any estimates of weather effects on abundance changes.

Here, we used a simulation following the advice of Dr. Dylan Childs to test whether the presence of density dependence might bias our estimation of weather effects on abundance changes. The aim here is to simulate a number (`n_sims = 1000`) of timeseries that have both a given strength of density dependence (`bDD`), and an effect of the environment on abundance changes (`bENV`). Given these parameters and some variance/noise, we simulate timeseries and then look at whether retrofitting a linear model to the timeseries enables us to accurately estimate the environmental effect. In other words, does the presence of density dependence prevent us from estimating the weather effects accurately.

The simulation has the following starting parameters:

```
t_step = 5       # timestep i.e. the length of each timeseries (here years) - testing with 5 to start with

p <- list()      # parameters for our models
p$b0 <- 1        # Initial value to feed to the density dependence
p$bDD <- 0.5     # coefficient for density dependence
p$bENV <- 0.2    # coefficient for the environment
p$sigma <- 1     # noise variance
p$sigmaENV <- 1  # environmental variance

x0 <- 2.0        # starting abundance

n_sims <- 1000   # number of simulations
```

Then we want to create functions that i) generate our timeseries/environmental data using the parameters defined above, ii) runs a linear model to estimate the effect of the environment and density dependence, and the iii) run the simulation and pull out the necessary coefficients from the retrofit model.

```
## Iteration function to create the two components of the model
iter <- function(p, x0, t_step) {
  x    <- numeric(t_step + 1)
  env  <- rnorm(t_step + 1, 0, p$sigmaENV)
  x[1] <- x0
  noise <- rnorm(t_step, 0, p$sigma)
  for (t  in seq_len(t_step)) {
    x[t+1] <- p$b0 + (1 - p$bDD) * x[t] + p$bENV * env[t] + noise[t]
  }
  list("x" = x, "env" = env)
}

## Fit the model function 
fit_mod <- function (p, x, beta) {
  unname(coef(lm(diff(x$x) ~ x$x[-(length(x$x))] + x$env[-(length(x$env))])))[beta]
}

## Function to run the simulation and pull out the necessary values
run_sim <- function(p, x0, t_step, n_sims, beta) {
  slopes <- numeric(n_sims)
  for (i in seq_len(n_sims)){
    x <- iter(p, x0, t_step)
    slopes[i] <- fit_mod(p, x, beta = beta)
  }
  slopes
}
```

Then, we ran this simulation for varying lengths of timeseries (5,10,20,50,100 and 1000 years) and looked at the coefficients generated in the model relative to the parameters put in. You can see the results here for estimated the environmental effect (`bENV`)

<img src="../plots/annual_abundance/Accounting_for_autocorrelation_testing/Environment_coefficient_simulation.jpeg" width="700"/>

and for the density dependence effect (`bDD`)

<img src="../plots/annual_abundance/Accounting_for_autocorrelation_testing/DD_coefficient_simulation.jpeg" width="700"/>

Where the dashed line is the simulation retrofit coefficient, and the solid line is the actual model parameter. You can see from the upper figure that for the environmental effect, even for short timeseries the estimation of an environmental effect is not hindered by the presence of density dependence. This is good news for our estimation of weather effects on population growth rates for the LPD data. However, this is not the case for density dependence estimation, which is really affected by the length of the timeseries. Very hard to accurately estimate density dependence from small timeseries.

Together with the results of section 2, this suggests that accounting for temporal autocorrelation in abundance (with lag 1) should be sufficient to estimate weather effects on abundance changes.

</details>

---

## 4. Calculating annual population growth rates
<details>
  <summary>Click here to expand</summary>

### `annual_population_growth_rate.R`

Using the abundance data that has been split in to consecutive blocks, we will calculate per-capita annual population growth rates, which we will then relate to annual weather data. This is calculated as r = N~t+1~/N~t~ , where N is the abundance on the ln scale at time t. We store the results in the `mammal` data frame for future reference.

```
# This will remove one year from each ID_block
mammal <- mam_IDblocks %>% 
  mutate(raw_abundance2 = raw_abundance + 1) %>% # For the 0 observations, adjusting all to non-zeros
  group_by(ID_block) %>% 
  group_modify(~{
    t0 <- .$raw_abundance2[-(length(.$raw_abundance2))] # get rid the last obs
    t1 <- .$raw_abundance2[-1]                        # get rid of the first obs
    
    mutate(., pop_growth_rate = c(log(t1/t0),NA))
  }) %>% 
  ungroup() %>% 
  filter(is.na(pop_growth_rate) == F)
```

This per-capita growth rate gives us a response variable with which to explore weather effects. The per-capita growth rates are distributed as follows:

<img src="../plots/annual_abundance/annual_population_growth_rates/pop_growth_rate_histogram.jpeg" width="600"/>

We can also explore the relationship between population growth rate and abundance, which gives us an indication of how density dependence may be acting across our abundance timeseries. Each point here is one year (t) of each record.

<img src="../plots/annual_abundance/annual_population_growth_rates/density_dependence_mam.jpeg" width="700"/>

We can see evidence of a slight negative trend between abundance and population growth, indicative of negative density dependence.

### Odd population growth rates and 0's in the raw data

On inspection of the histogram of our per-capita population growth rates, while the vast majority of the data falls within a reasonable range, we can see that there are a number of odd population growth rates, at both tails of the distribution. In order to explore this in a little bit more detail, lets have a look at the density of these per-capita population growth rates.

<img src="../plots/annual_abundance/annual_population_growth_rates/population_growth_rate_distribution.jpeg" width="700"/>

We can see that the vast vast majority of the data is centered around a population growth rate of 1 indicating no annual change. The dashed lines here are the 1% and 99% quantiles of the population growth rates, indicating that 98% of the data is occurring between a population growth rate of **0.54 and 1.9**. This makes sense - in a single year we perhaps do not expect that many timeseries observations display a change in abundance great than a factor of 2, which if true would indicate very very rapid population changes.

We can see that there are a number of points outside of this range however, and also some annual population changes well outside of this range that deserve some exploration. So, we can pull out all of the timeseries observations where there these per-capita population growth rates lying outside of this range:

```
# The quantiles
mamquant <- tibble(quant = quantile(mammal$pop_growth_rate, probs = c(0.01,0.99)))

# Pull out the ID blocks below-above the 1%-99% quantiles 
unusual_pop_growth <- mammal %>% 
  filter(pop_growth_rate > mamquant$quant[2] | 
           pop_growth_rate < mamquant$quant[1]) 

# and the raw abundance changes
unusual_dat <- mam_IDblocks %>% 
  filter(ID_block %in% unusual_pop_growth$ID_block == TRUE)
```

And then we can have a look at the timeseries where there are these rapid population changes. Below is the timeseries of ln Abundance for all records that have at least one population growth rate outside of the normal range. Each panel is a record, and the black line is their full timeseries. The coloured points are the years with odd population growth rate observations, where the colour of the point is the population growth rate.

<img src="../plots/annual_abundance/annual_population_growth_rates/unusual_abundance_changes.jpeg" width="900"/>

For quite a few of these records, we can see that the rapid increases/decreases are real and for the populations, and probably just because of quite rapid fluctuations in the population from year to year. A good example of a record like this is **18740_1**. However, there are several suspicious examples. It seems there are some records where the ln abundance drops to 1, which indicates a **raw abundance of 0**. Here, because of the log transformation, this means there are huge huge jumps in the population when it drops to 0. Furthermore, we are unsure whether these are *true* 0s, or instead reflect sampling variation i.e. years where individuals were not recorded and should have missing data but have been recorded as 0 by accident. A good example of a key candidate for this is record **9888_1**, where there are biannual drops in the population to 0. This is highly unlikely, especially because this is a population of lions, for which abundance is expected to be more stable that 5-fold changes in abundance in a single year.

The key seems to be in abundance timeseries where there are 0s. Indeed, when we restrict to look at only records with 0s in the raw abundance data, the very same populations come up. Here are the abundance timeseries for records with raw 0s:

<img src="../plots/annual_abundance/annual_population_growth_rates/raw0_abundance.jpeg" width="1000"/>

The 0 abundances map very well to the extreme population growth rates. Now we need to find a reliable way of removing suspicious records, because unfortunately detailed information about the sampling regime is not given in the meta data `mam_meta`. 

The really extreme population growth rates seem to be coming from records that have a large proportion of 0s, and those with long continuous periods with raw 0s. Lets look at the proportion of raw abundance data that is 0 in all of these records, with the 75% quantile showed with the vertical line:

<img src="../plots/annual_abundance/annual_population_growth_rates/raw0_proportion.jpeg" width="800"/>

Similar suspects are present, including **9888_1**. So we will use this criteria to restrict the raw population growth data for our models. More specifically, **we removed raw 0s from records where there was >30% of the record with 0s**. So, we will remove any observations with a 0. But then, we need to retain our criteria of >10 years of population data for our analysis. This means that the majority of these suspect records need to be deleted, except for **9887_1**, for which the first half only needs to be deleted.

```
# Full Record Unreliable
full_unreliable <- c("10096_1", "17787_1", "23155_2", "23156_2", 
                     "5835_1", "8103_1", "9863_1", "9888_1")

# Part of record unreliable - early part of record 9887_1
years_unreliable_9887_1 <- (mam_IDblocks %>% filter(ID_block == "9887_1" 
                                                    & raw_abundance == 0))$year

## Restrict population growth data 
mammal <- mammal %>% 
  filter(!((ID_block %in% full_unreliable == TRUE) |
           (ID_block == "9887_1" & year %in% years_unreliable_9887_1 == TRUE)))
```

This gives `mammal`, our prepared data for analysis with annual per-capita population growth rates from ln Abundance

</details>





