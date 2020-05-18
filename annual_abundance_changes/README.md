# Detrending annual abundance change data
#### 2020-05-06 
#### John Jackson

---

This mardown is intended as an accompaniment to the scripts contained within the directory `annual_abundance_changes`, to walk through the process of detrending annual abundance data from the Living Planet Database for the terrestrial mammals. Please refer to the scripts mentioned in each section of the markdown for full details on each section.

There are 2 main sections and scripts:

## 1. Exploring and accounting for gaps in the abundance timeseries
<details>
  <summary>Click here to expand</summary>

### `timeseries_gap_exploration.R`

For our record question of how weather influencs annual population changes in vertebrates, one potential problem with the vertebrate abundance data from the Living Planet Database is that there are gaps in the timeseries. The number of these gaps and the way we deal with them is important. This first section is intended to explore the pervasiveness of these gaps across our records, and deal with them in an appropriate way for further analysis.

We first have to restrict the data to only include records that have sufficient data with which to explore annual changes abundance in relation to CHELSA weather data. We only include observations that overlap with the CHELSA data i.e. between 1979-2013, and those that have 5 or more years of abundance data:

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

<img src="../plots/annual_abundance/fossa_timeseries.jpeg" width="600"/>

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
  dplyr::select(ID, Binomial, Order, scaled_abundance, 
                year, block, max_block) %>% 
  left_join(x = ., y = dplyr::select(mam_gaps, -c(Binomial, no_consecutive_blocks)),
            by = "ID") %>% 
  arrange(desc(longest_block)) %>% 
  mutate(ID = factor(ID, levels = unique(.$ID)))
```
### Time series gap summaries

Encouragingly, the majority of the data is occuring in continuous blocks with 1-year transitions. Here we have the distribution of the proportion of each record that is occuring in 1-year transitions. We can see that the majority have all their data as 1-year transitions. Furthermore >76% (870) of the records have more than 2 thirds of their data in 1-year transitions.

<img src="../plots/annual_abundance/Proportion_1year_transitions.jpeg" width="700" />

However, this doesn't quite give us the full picture because we also need to know how many blocks each timeseries occurs in. Here we plot the number of blocks that each record occurs in against the number of years in its longest block. The colour denotes the proportion of the timeseries occuring in 1-year transitions.

<img src="../plots/annual_abundance/Consecutive_blocks.jpeg" width="700" />

So it does appear that there are some records that are primarily in timeseries with 1-year transitions (lighter colours), but do occur over quite a few blocks of observations. We can also plot these blocks of observations as timelines, where we see the years of data for each record ID. I have split these up based on the number of blocks that the timeseries occurs in for ease. Here first you have the records that are just occuring in 1 consecutive block. Points and lines indicate where there is data for each record ID (row).

<img src="../plots/annual_abundance/record_timelines/1_consecutive_block_timeline.jpeg" width="700" />

These are the 'gold standard' records that occur solely in one consecutive chain of annual observations (with more than 5 years of data). However, the picture becomes a little bit more complex when we look at records that occur in a greater number of blocks. Here you can see the records that occur in 5 blocks.

<img src="../plots/annual_abundance/record_timelines/5_consecutive_block_timeline.jpeg" width="700" />

We can see here that there are scenarios where there are longer consecutive blocks of observations, with smaller satellite blocks that have fewer observations. Furthermore, we can see in records with more blocks, there are situations where there are several separate blocks of 1-year transitions, but that many blocks have less than 5 observations.

### Data selection

We are selecting data based on the sizes of the blocks for each record - We only want to retain blocks within a record that have 5 or more consecutive annual observations.

```
# IDs and blocks that we want to keep - 901 out of 2756 ID-block combinations
ID_block_keep <- mam_blocks %>% 
  mutate(ID = as.numeric(as.character(ID))) %>% 
  group_by(ID, block) %>% 
  summarise(ID_block = paste0(ID[1],"_",block[1]),
            block_keep = if_else(n() >= 5, 1, 0)) %>%
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

This gives us a final dataset with which to assess how weather affects changes in abundance, stored in `mammal`. Each Study has atleast 5 years of raw data, and consecutive blocks within the study with at least five years of data within them. Here is a comparison to the initial raw data:

![](../plots/annual_abundance/data_summary.jpeg)

This equates to ~51% of the initial observations, 33% of the initial records, and 48% of the species.

</details>


---

## 2. Detrending abundance timseries and calculating population growth rates
<details>
  <summary>Click here to expand</summary>

### `detrend_population_growth_rate.R`

In this section, using the abundance data that has been split in to consecutive blocks, we will detrend each consecutive block to extract the residual abundance, and then calculate per-captita population growth rates from these detrended abundances. The rationale for detrending the raw abundance data was to focus on annual population changes that were not as influenced by underlying trends driving population dynamics such as habitat loss or other human disturbances. By detrending, we are focussing on residual annual population changes, which we expect can be driven by changes in weather.

Here we perform a simple linear, vertical detrend of the scaled abundance data. This is visualised (and compared to orthogonal detrending) in `testing_detrending.R`. If we simulate a timeseries with observations of abundance (y) at different timepoints (x), providing it is justified to use a linear model to capture the timeseries trend, which we assume here, we can fit a linear trend to the data (left). Then, to detrend here we take the vertical deviations from the fitted line, or the residuals from the linear model. Vertical residuals are justified here because there isn't error in our independent variable (year), or at least far far less than there is in abundance. You can see an example of this below, with a simulated timeseries and then green dashed line segments to visualise the vertical residuals used here (right).

<img src="../plots/annual_abundance/detrending_linear.jpeg" width="700" />

We are repeating this detrending for all consecutive blocks of each record ID separately. This is step 2 from `detrend_population_growth_rate.R`, detrending for each block_ID in `mam_IDblocks`. Importantly however for population growth rate calculations, we are then scaling the residual abundances calculated to be centered around 10 <- IS THIS JUSTIFIED?????

```
# extracting the residuals after fitting a linear trend to each timeseries
mam_detrend <- mam_IDblocks %>% 
  group_by(ID_block) %>% 
  group_modify(~{
    mod = lm(scaled_abundance ~ year, data = .) 
    
    resid_ab = mod$residuals + 10 # Centre around 10 for sensible population growth rate calculations
    
    mutate(., residual_abundance = resid_ab,
           coef = mod$coefficients[2])
  }) %>% 
  ungroup()
```

</details>





