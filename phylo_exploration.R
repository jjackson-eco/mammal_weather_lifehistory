#####################################################
##                                                 ##
##     Global climate and population dynamics      ##
##                                                 ##
##            Mammal phylogeny data                ##
##                                                 ##
##             October 14th 2020                   ##
##                                                 ##
#####################################################
rm(list = ls())
options(width = 100)

library(tidyverse)
library(ape)
library(phangorn)
library(phytools)
library(caper)

##__________________________________________________________________________________________________
#### 1. Load data ####

load("../Phlogenetics/mam_lpd_pruned.RData", verbose = T)

## Have a look at the data
str(mamMCC_pruned)
glimpse(LPD_tree_update)

##__________________________________________________________________________________________________
#### 2. Simple plot ####

jpeg(file = "plots/mam_LPD_tree.jpeg", width = 7, height = 7, units = "in", res = 800)
plot(mamMCC_pruned, cex = 0.2, type = "f", label.offset = 0.3)
dev.off()

