#####################################################
##                                                 ##
##     Global climate and population dynamics      ##
##                                                 ##
##   Annual Abundance Timeseries Gap Exploration   ##
##                                                 ##
##               April 23rd 2020                   ##
##                                                 ##
#####################################################
rm(list = ls())
options(width = 100)

library(tidyverse)
library(gridExtra)
library(grid)

## Load in the raw mammal data
load("../rawdata/mam.RData", verbose = T)

mam_raw <- mam # for summary later
summarise(mam_raw, n = n(), n_records = n_distinct(ID))

##__________________________________________________________________________________________________
#### 1. Set up data ####

# Keeping only years from CHELSA and only 5 years of study beyond that.
mam <- mam %>% 
  filter(year >= 1979 & year <= 2013) %>% # important to do this first
  group_by(ID) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  filter(n >= 5) %>% 
  dplyr::select(-n)

n_distinct(mam$ID) # Initial number of records

# Check for duplicated years
# dupyer <- mam %>% 
#   group_by(ID) %>% 
#   summarise(Binomial = Binomial[1],
#             duplicate_years = n() -n_distinct(year))

##__________________________________________________________________________________________________
#### 2. Testing - Fossa  ####

test <- filter(mam, ID == 25396)

ggplot(test, aes(x = year, y = scaled_abundance)) + 
  geom_point(size = 2) +
  labs(x = "Year", y = "Scaled abundance") +
  theme_bw(base_size = 12) +
  ggsave("plots/annual_abundance/fossa_timeseries.jpeg",
         width = 5, height = 3.5, units = "in", dpi = 400)

tyear <- test$year
diff_tyear <- diff(tyear)
cumsum_blocks <- cumsum(c(1, diff_tyear != 1))
cumsum_blocks_tab <- table(cumsum_blocks)

no_consecutive_blocks <- n_distinct(cumsum_blocks)
prop_1year_transitions <- sum(diff_tyear == 1)/ length(diff_tyear)

longest_block <- max(cumsum_blocks_tab)

# # interesting function
# consecutive_blocks <-  length(split(tyear, cumsum(c(1, diff(tyear) != 1))))

##__________________________________________________________________________________________________
#### 3. Exploring gaps for all records ####

mam_gaps <- mam %>% 
  group_by(ID) %>% 
  group_modify(~{
    cyears = .$year
    diff_cyears = diff(cyears)
    cumsum_blocks = cumsum(c(1, diff_cyears != 1))
    
    summarise(., Binomial = Binomial[1],
              record_length = length(cyears),
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

##__________________________________________________________________________________________________
#### 4. Exploration plots ####

#____________________________________________________________
## 4a. Proportion of the records in 1-year transitions

## This is good news - most records are predominantly in 1 year blocks
prophist <- ggplot(mam_gaps, aes(x = prop_1year_transitions)) +
  geom_histogram(bins = 10, fill = "lightblue", 
                 colour = "black", size = 0.1) +
  labs(x = "Proportion of record in 1-year transitions",
       y = "Number of records") +
  theme_bw(base_size = 11)

propbx <- ggplot(mam_gaps, aes(y = prop_1year_transitions)) +
  geom_boxplot(fill = "lightblue", 
                colour = "black") +
  scale_x_continuous(limits = c(-0.75,0.75)) +
  labs(x = " ", 
       y = "Proportion of record in 1-year transitions") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(grid.arrange(prophist, propbx, ncol = 2, widths = c(2,1)),
       filename = "plots/annual_abundance/Proportion_1year_transitions.jpeg",
        width = 6, height = 4, units = "in", dpi = 400)

# 25th percentile is 0.66667 i.e. 2/3rds of the record in 1-year transitions
quantile(mam_gaps$prop_1year_transitions,probs = c(0.025,0.25,0.5,0.75,0.975))

# >76% (870) of the records have more than 2 thirds of its data in 1-year transitions
mam_gaps %>%
  filter(prop_1year_transitions >= 2/3) %>% 
  summarise(tot = n(), full = nrow(mam_gaps),
            prop = tot/full) 

#____________________________________________________________
## 4b. Consecutive blocks

ggplot(mam_gaps, aes(x = no_consecutive_blocks, y = longest_block, 
                     colour = prop_1year_transitions)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_viridis_c(begin = 0.1,end = 0.9,
                        name = "Proportion\nof record in\n1-year transitions") +
  scale_x_continuous(breaks = 1:10, labels = 1:10) +
  scale_y_continuous(breaks = seq(0,35, by = 5), labels = seq(0,35, by = 5)) +
  labs(x = "Number of blocks", y = "Longest block (years)") +
  guides(colour = guide_colorbar(barheight = 8, barwidth = 1.4)) +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank()) +
  ggsave(filename = "plots/annual_abundance/Consecutive_blocks.jpeg",
         width = 6, height = 4, units = "in", dpi = 400)
  
ggplot(filter(mam_gaps, prop_1year_transitions >= 2/3), 
       aes(x = no_consecutive_blocks, y = longest_block, 
                     colour = prop_1year_transitions)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_viridis_c(begin = 0.1,end = 0.9,
                        name = "Proportion\nof record in\n1-year transitions") +
  scale_x_continuous(breaks = 1:10, labels = 1:10) +
  scale_y_continuous(breaks = seq(0,35, by = 5), labels = seq(0,35, by = 5)) +
  labs(x = "Number of blocks", y = "Longest block (years)") +
  guides(colour = guide_colorbar(barheight = 8, barwidth = 1.4)) +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank()) +
  ggsave(filename = "plots/annual_abundance/Consecutive_blocks_highproportion.jpeg",
         width = 6, height = 4, units = "in", dpi = 400)

#____________________________________________________________
## 4c. record timelines

for(i in 1:10){
  
  cdat = filter(mam_blocks, max_block == i)
  
  if(i == 1){tit = paste0("Records with ", i, " consecutive block")}
    else{tit = paste0("Records with ", i, " consecutive blocks")}
  
  ggplot(cdat, aes(x = ID, y = year, 
                         colour = factor(block))) + 
    geom_line(size = i/2, show.legend = FALSE) + 
    geom_point(size = i/2) +
    scale_colour_viridis_d(begin = 0.1, end = 0.9, guide = F) +
    scale_y_continuous(breaks = seq(1980,2015,by = 5),
                       labels = c(1980, "", 1990, "", 
                                  2000, "", 2010, "")) +
    labs(x = "Record ID", y = "Year",
         title = tit) +
    coord_flip() +
    theme_bw(base_size = 28) +
    theme(panel.grid = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          strip.background = element_blank()) +
    ggsave(filename = paste0("plots/annual_abundance/record_timelines/", i,
                             "_consecutive_block_timeline.jpeg"),
           width = 210, height = 297, units = "mm", dpi = 400)
}

##__________________________________________________________________________________________________
#### 5. Data restriction and save ####

#____________________________________________________________
## 5a. Restricting by the blocks
# Want only blocks from records that have 5 or more observations

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

# This leaves us with 844 records or 33% of the initial 2539 records, 
# and 10,489 of the 20,379 raw observations, equating to about 51%

# Summarising this in a nice table
mam_datasum <- data.frame(Dataset = c("Raw data", "Study data"),
                          Observations  = c(nrow(mam_raw), 
                                            nrow(mam_IDblocks)),
                          Records = c(n_distinct(mam_raw$ID), 
                                      n_distinct(mam_IDblocks$ID)),
                          Species = c(n_distinct(mam_raw$Binomial),
                                      n_distinct(mam_IDblocks$Binomial)))

General_sum <- tableGrob(mam_datasum, rows = NULL, theme = ttheme_minimal(base_size = 16))
ggsave(grid.arrange(General_sum), 
       filename = "plots/annual_abundance/data_summary.jpeg",
       width = 7, height = 3, units = "in",
       dpi = 400)

#____________________________________________________________
## 5b. Saving the data

save(mam_IDblocks, file = "../rawdata/mam_IDblocks.RData")


