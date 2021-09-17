####################################################
##                                                ##
##     Global climate and population dynamics     ##
##                                                ##
##      Weather 10 Popgrowth timeseries data      ##
##                                                ##
##                    Sept 2021                   ##
##                                                ##
####################################################

rm(list = ls())
options(width = 100)

### State-space model approach based on 10 timeseries records, the longest 5 and 5 with 9 years of data.
### Aim: 1) Can we easily add the weather covariates into the SSM? 
###      2) Are the weather coefficients in long and/or short time series similar to the estimates from GAMM models? (growth ~ weather) 

library(jagsUI) 
library(mgcv)  

##****************************************
## Look at the data (after loading data)
##****************************************

load("../rawdata/mammal_weather_CL.RData")
head(mammal_weather_CL)

##**********************************************************
## State-space model - Plotting the observed and estimated N
## for each population separately (see pdf)
## No covariate included
##**********************************************************

IDs <- unique(mammal_weather_CL$ID)

record <- mammal_weather_CL[mammal_weather_CL$ID==IDs[1],]

# Specify model in BUGS language
cat(file="ssm.jags", "
  model {

  # Priors and constraints
  logN[1] ~ dnorm(5.6, 0.01)           # Prior for initial population size
  mean.r ~ dnorm(0, 0.001)             # Prior for mean growth rate
  sigma.proc ~ dunif(0, 1)             # Prior for sd of state process
  sigma2.proc <- pow(sigma.proc, 2)
  tau.proc <- pow(sigma.proc, -2)
  sigma.obs ~ dunif(0, 1)              # Prior for sd of observation process
  sigma2.obs <- pow(sigma.obs, 2)
  tau.obs <- pow(sigma.obs, -2)

  # Likelihood
  # State process
  for (t in 1:(T-1)){
    r[t] ~ dnorm(mean.r, tau.proc)
    logN[t+1] <- logN[t] + r[t]
  }
  # Observation process
  for (t in 1:T) {
    y[t] ~ dnorm(logN[t], tau.obs)
  }

  # Population sizes on real scale
  for (t in 1:T) {
    N[t] <- exp(logN[t])
  }
  }
")

# Bundle data
jags.data <- list(y=log(record$raw_abundance_no0), T=length(record$year))

# Initial values
inits <- function(){list(sigma.proc=runif(1, 0, 1), mean.r=rnorm(1), 
                         sigma.obs=runif(1, 0, 1), logN=c(rnorm(1, 5.6, 0.1), 
                                                          rep(NA, (length(record$year)-1))))}

# Parameters monitored
parameters <- c("r", "mean.r", "sigma2.obs", "sigma2.proc", "N")

# MCMC settings
ni <- 200000
nt <- 6
nb <- 100000
nc <- 3
na <- 5000

# Call JAGS from R (BRT 3 min)
record.ssm <- jags(jags.data, inits, parameters, "ssm.jags", 
                   n.chains=nc, n.thin=nt, n.iter=ni, 
                   n.burnin=nb, n.adapt=na)

# Draw figure
fitted <- lower <- upper <- numeric()
year <- min(record$year):max(record$year)
n.years <- length(record$year) 
for (i in 1:n.years){
    fitted[i] <- mean(record.ssm$sims.list$N[,i])
    lower[i] <- quantile(record.ssm$sims.list$N[,i], 0.025)
    upper[i] <- quantile(record.ssm$sims.list$N[,i], 0.975)}
m1 <- min(c(fitted, record$raw_abundance_no0, lower), na.rm=TRUE)
m2 <- max(c(fitted, record$raw_abundance_no0, upper), na.rm=TRUE)
  
par(mar=c(4.5, 4, 1, 1))
plot(0, 0, ylim=c(m1, m2), xlim=c(1, n.years), ylab="Population size", xlab="Year", col="black", type="l", lwd=2, axes=FALSE, frame=FALSE)
axis(2, las=1)
axis(1, at=1:n.years, labels=year)
polygon(x=c(1:n.years, n.years:1), y=c(lower, rev(upper)), col="gray90", border="gray90")
points(record$raw_abundance_no0, type="l", col="black", lwd=2)
points(fitted, type="l", col="blue", lwd=2)
legend(x=1, y=max(record$raw_abundance_no0), legend=c("Counts", "Estimates"), lty=c(1, 1), lwd=c(2, 2), col=c("black", "blue"), bty="n", cex=1)


##********************************************************
## State-space model - Weather coefficients ('beta')
## and comparison with GAMM coefficients for 
## temp and prec anomalies
##
## All coefficients are saved in "ssmresults" 
## for all populations, using a loop
##*******************************************************

ssmresults <- NULL
IDs <- unique(mammal_weather_CL$ID)

for (i in 1:length(IDs)){
  record <- mammal_weather_CL[mammal_weather_CL$ID==IDs[i],]
  
  #############################
  ####     TEMPERATURE     ####
  #############################
  
  ################
  ###SSM
  
  # Specify model in BUGS language
  cat(file="ssmtemp.jags", "
  model {
  
  # Priors and constraints
  logN[1] ~ dnorm(5.6, 0.01)           # Prior for initial population size
  b0 ~ dnorm(0, 0.001)                 # Prior for mean growth rate
  beta ~ dnorm(0, 0.001)               # Prior for slope parameter
  sigma.proc ~ dunif(0, 1)             # Prior for sd of state process
  sigma2.proc <- pow(sigma.proc, 2)
  tau.proc <- pow(sigma.proc, -2)
  sigma.obs ~ dunif(0, 1)              # Prior for sd of observation process
  sigma2.obs <- pow(sigma.obs, 2)
  tau.obs <- pow(sigma.obs, -2)
  
  # Likelihood
  # State process
  for (t in 1:(T-1)){
    r[t] <- b0 + beta * x[t] + epsilon[t]   # Linear model for the population growth rate
    epsilon[t] ~ dnorm(0, tau.proc)         # Random noise of the population growth rate
    logN[t+1] <- logN[t] + r[t]
  }
  # Observation process
  for (t in 1:T) {
    y[t] ~ dnorm(logN[t], tau.obs)
  }
  
  # Population sizes on real scale
  for (t in 1:T) {
    N[t] <- exp(logN[t])
  }
  }
  ")
  
  # Bundle data
  jags.data <- list(y=log(record$raw_abundance_no0), T=length(record$year), x=record$mean_temp_anomaly)
  
  # Initial values
  inits <- function(){list(sigma.proc=runif(1, 0, 1), sigma.obs=runif(1, 0, 1), b0=rnorm(1),logN=c(rnorm(1, 5.6, 0.1), rep(NA, (length(record$year)-1))))}
  
  # Parameters monitored
  parameters <- c("r",  "sigma2.obs", "sigma2.proc", "N", "b0", "beta")
  
  
  # MCMC settings
  ni <- 200000
  nt <- 6
  nb <- 100000
  nc <- 3
  na <- 5000
  
  # Call JAGS from R (BRT 3 min)
  ssmtemptot <- jags(jags.data, inits, parameters, "ssmtemp.jags", 
                     n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, n.adapt=na)
  #print(ssmtemptot, digits=3)
  tempbetaSSM <- ssmtemptot$mean$beta
  
  
  ################
  ###GAMM
  mod_temp = gamm(pop_growth_rate ~ s(year, bs = "tp", k = 5) + mean_temp_anomaly,
                  data=record,family = gaussian,
                  correlation = corARMA(form = ~ year, p = 1),
                  method = "REML")
  tempbetaGAMM = coef(mod_temp$gam)
  
  
  
  
  #############################
  ####    PRECIPITATION    ####
  #############################
  
  ################
  ###SSM
  
  # Specify model in BUGS language
  cat(file="ssmprec.jags", "
  model {

  # Priors and constraints
  logN[1] ~ dnorm(5.6, 0.01)           # Prior for initial population size
  b0 ~ dnorm(0, 0.001)                 # Prior for mean growth rate
  #beta ~ dnorm(0, 0.001)               # Prior for slope parameter
  beta ~ dunif(-10, 10)               # Prior for slope parameter

  sigma.proc ~ dunif(0, 1)             # Prior for sd of state process
  sigma2.proc <- pow(sigma.proc, 2)
  tau.proc <- pow(sigma.proc, -2)
  sigma.obs ~ dunif(0, 1)              # Prior for sd of observation process
  sigma2.obs <- pow(sigma.obs, 2)
  tau.obs <- pow(sigma.obs, -2)

  # Likelihood
  # State process
  for (t in 1:(T-1)){
    r[t] <- b0 + beta * x[t] + epsilon[t]   # Linear model for the population growth rate
    epsilon[t] ~ dnorm(0, tau.proc)         # Random noise of the population growth rate
    logN[t+1] <- logN[t] + r[t]
  }
  # Observation process
  for (t in 1:T) {
    y[t] ~ dnorm(logN[t], tau.obs)
  }

  # Population sizes on real scale
  for (t in 1:T) {
    N[t] <- exp(logN[t])
  }
  }
  ")
  
  
  
  # Bundle data
  jags.data <- list(y=log(record$raw_abundance_no0), T=length(record$year), x=record$mean_precip_anomaly)
  
  # Initial values
  inits <- function(){list(sigma.proc=runif(1, 0, 1),  sigma.obs=runif(1, 0, 1), b0 = rnorm(1), logN=c(rnorm(1, 5.6, 0.1), rep(NA, (length(record$year)-1))))}
  
  # Parameters monitored
  parameters <- c("r", "sigma2.obs", "sigma2.proc", "N", "b0", "beta")
  
  # Call JAGS from R (BRT 3 min)
  ssmprectot <- jags(jags.data, inits, parameters, "ssmprec.jags", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, n.adapt=na)
  precbetaSSM <- ssmprectot$mean$beta
  
  
  
  ################
  ###GAMM
  # Precipitation + dealing with NA values
  if(length(which(is.na(record$mean_precip_anomaly) == T)) == 0){
    mod_precip = gamm(pop_growth_rate ~ s(year, bs = "tp", k = 5) + mean_precip_anomaly,
                      data = record, family = gaussian,
                      correlation = corARMA(form = ~ year, p = 1),
                      method = "REML")
    precbetaGAMM = coef(mod_precip$gam)}else{precbetaGAMM = rep(NA,15)}     # Arbitrary long NA vector
  
  
  ### Gather the data
  nyear <-length(record$year)
  essai <- cbind(record[1,1:7], nyear, tempbetaSSM, tempbetaGAMM[2], precbetaSSM, precbetaGAMM[2])
  rownames(essai) <- NULL
  ssmresults <- rbind(ssmresults,essai)
  
}

ssmresults













##***************************************************************
## State-space model - With a different prior for 'beta'
## A uniform distribution instead of a normal distribution
## and see how priors influenced the results
##***************************************************************

ssmresults_uni <- NULL
IDs <- unique(mammal_weather_CL$ID)

for (i in 1:length(IDs)){
  record <- mammal_weather_CL[mammal_weather_CL$ID==IDs[i],]
  
  #############################
  ####     TEMPERATURE     ####
  #############################
  
  ################
  ###SSM
  
  # Specify model in BUGS language
  cat(file="ssmtemp.jags", "
model {

# Priors and constraints
logN[1] ~ dnorm(5.6, 0.01)           # Prior for initial population size
b0 ~ dnorm(0, 0.001)                 # Prior for mean growth rate
beta ~ dunif(-10, 10)                # Prior for slope parameter
sigma.proc ~ dunif(0, 1)             # Prior for sd of state process
sigma2.proc <- pow(sigma.proc, 2)
tau.proc <- pow(sigma.proc, -2)
sigma.obs ~ dunif(0, 1)              # Prior for sd of observation process
sigma2.obs <- pow(sigma.obs, 2)
tau.obs <- pow(sigma.obs, -2)

# Likelihood
# State process
for (t in 1:(T-1)){
  r[t] <- b0 + beta * x[t] + epsilon[t]   # Linear model for the population growth rate
  epsilon[t] ~ dnorm(0, tau.proc)         # Random noise of the population growth rate
  logN[t+1] <- logN[t] + r[t]
}
# Observation process
for (t in 1:T) {
  y[t] ~ dnorm(logN[t], tau.obs)
}

# Population sizes on real scale
for (t in 1:T) {
  N[t] <- exp(logN[t])
}
}
")
  
  # Bundle data
  jags.data <- list(y=log(record$raw_abundance_no0), T=length(record$year), x=record$mean_temp_anomaly)
  
  # Initial values
  inits <- function(){list(sigma.proc=runif(1, 0, 1), sigma.obs=runif(1, 0, 1), b0=rnorm(1),logN=c(rnorm(1, 5.6, 0.1), rep(NA, (length(record$year)-1))))}
  
  # Parameters monitored
  parameters <- c("r",  "sigma2.obs", "sigma2.proc", "N", "b0", "beta")
  
  
  # MCMC settings
  ni <- 200000
  nt <- 6
  nb <- 100000
  nc <- 3
  na <- 5000
  
  # Call JAGS from R (BRT 3 min)
  ssmtemptot <- jags(jags.data, inits, parameters, "ssmtemp.jags", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, n.adapt=na)
  #print(ssmtemptot, digits=3)
  tempbetaSSM <- ssmtemptot$mean$beta
  
  
  ################
  ###GAMM
  mod_temp = gamm(pop_growth_rate ~ s(year, bs = "tp", k = 5) + mean_temp_anomaly,
                  data=record,family = gaussian,
                  correlation = corARMA(form = ~ year, p = 1),
                  method = "REML")
  tempbetaGAMM = coef(mod_temp$gam)
  
  
  
  
  #############################
  ####    PRECIPITATION    ####
  #############################
  
  ################
  ###SSM
  
  # Specify model in BUGS language
  cat(file="ssmprec.jags", "
  model {

  # Priors and constraints
  logN[1] ~ dnorm(5.6, 0.01)           # Prior for initial population size
  b0 ~ dnorm(0, 0.001)                 # Prior for mean growth rate
  beta ~ dunif(-10, 10)                # Prior for slope parameter

  sigma.proc ~ dunif(0, 1)             # Prior for sd of state process
  sigma2.proc <- pow(sigma.proc, 2)
  tau.proc <- pow(sigma.proc, -2)
  sigma.obs ~ dunif(0, 1)              # Prior for sd of observation process
  sigma2.obs <- pow(sigma.obs, 2)
  tau.obs <- pow(sigma.obs, -2)

  # Likelihood
  # State process
  for (t in 1:(T-1)){
    r[t] <- b0 + beta * x[t] + epsilon[t]   # Linear model for the population growth rate
    epsilon[t] ~ dnorm(0, tau.proc)         # Random noise of the population growth rate
    logN[t+1] <- logN[t] + r[t]
  }
  # Observation process
  for (t in 1:T) {
    y[t] ~ dnorm(logN[t], tau.obs)
  }

  # Population sizes on real scale
  for (t in 1:T) {
    N[t] <- exp(logN[t])
  }
  }
  ")
  
  
  
  # Bundle data
  jags.data <- list(y=log(record$raw_abundance_no0), T=length(record$year), x=record$mean_precip_anomaly)
  
  # Initial values
  inits <- function(){list(sigma.proc=runif(1, 0, 1),  sigma.obs=runif(1, 0, 1), b0 = rnorm(1), logN=c(rnorm(1, 5.6, 0.1), rep(NA, (length(record$year)-1))))}
  
  # Parameters monitored
  parameters <- c("r", "sigma2.obs", "sigma2.proc", "N", "b0", "beta")
  
  # Call JAGS from R (BRT 3 min)
  ssmprectot <- jags(jags.data, inits, parameters, "ssmprec.jags", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, n.adapt=na)
  precbetaSSM <- ssmprectot$mean$beta
  
  
  
  ################
  ###GAMM
  # Precipitation + dealing with NA values
  if(length(which(is.na(record$mean_precip_anomaly) == T)) == 0){
    mod_precip = gamm(pop_growth_rate ~ s(year, bs = "tp", k = 5) + mean_precip_anomaly,
                      data = record, family = gaussian,
                      correlation = corARMA(form = ~ year, p = 1),
                      method = "REML")
    precbetaGAMM = coef(mod_precip$gam)}else{precbetaGAMM = rep(NA,15)}     # Arbitrary long NA vector
  
  
  ### Gather the data
  nyear <-length(record$year)
  essai <- cbind(record[1,1:7], nyear, tempbetaSSM, tempbetaGAMM[2], precbetaSSM, precbetaGAMM[2])
  rownames(essai) <- NULL
  ssmresults_uni <- rbind(ssmresults_uni,essai)
  
}

ssmresults_uni



save(mammal_weather_CL, ssmresults, ssmresults_uni, file= "SSMresults.Rdata")







###ssmresults_uni
ssmresults_uni
#      ID         Family               Binomial           gbif_species                                                        Biome  Latitude  Longitude nyear tempbetaSSM tempbetaGAMM[2] precbetaSSM precbetaGAMM[2]
#1    102        Bovidae    Hippotragus_equinus    Hippotragus equinus Tropical and subtropical grasslands, savannas and shrublands -23.83333   31.50000     9  -1.5397097      -1.7980647  0.61698205      0.68961838
#2   3412        Canidae          Canis_latrans          Canis latrans                Temperate grasslands, savannas and shrublands  48.10611  -95.22806     9   2.6758600       1.6724629 -0.22457627      0.15757768
#3   3413        Canidae          Canis_latrans          Canis latrans                        Temperate broadleaf and mixed forests  47.95500  -94.78639     9  -0.1345889      -0.3553621 -2.03368776     -0.33689653
#4   3414        Canidae          Canis_latrans          Canis latrans                        Temperate broadleaf and mixed forests  47.31667  -92.54750     9  -0.6572392      -2.4498224  1.40350195     -0.21480807
#5   3415        Canidae          Canis_latrans          Canis latrans                Temperate grasslands, savannas and shrublands  45.33611  -94.40472     9  -0.2044879      -1.4640354 -0.10879361     -0.34475118
#6  18764      Leporidae       Lepus_americanus       Lepus americanus                                         Boreal forests/taiga  61.02759 -138.41045    34   0.4354323       0.2458739 -0.01761835     -0.02976037
#7  18887     Cricetidae Peromyscus_maniculatus Peromyscus maniculatus                                         Boreal forests/taiga  61.02759 -138.41045    34  -0.2164236      -0.1662930 -0.25571488     -0.24333106
#8  19103 Antilocapridae  Antilocapra_americana  Antilocapra americana                Temperate grasslands, savannas and shrublands  47.63095 -100.45888    34   0.6464810       0.8111799 -0.27253255     -0.33598661
#9  19783       Cervidae      Rangifer_tarandus      Rangifer tarandus                                         Boreal forests/taiga  48.62557  -56.42212    34  -0.1252751      -0.1636454  0.05642041      0.02840967
#10 22039       Cervidae Odocoileus_virginianus Odocoileus virginianus                        Temperate broadleaf and mixed forests  46.56531  -66.46191    34   0.9311289       0.9327810  0.11543555      0.01756898




ssmresults
#      ID         Family               Binomial           gbif_species                                                        Biome  Latitude  Longitude nyear tempbetaSSM tempbetaGAMM[2] precbetaSSM precbetaGAMM[2]
#1    102        Bovidae    Hippotragus_equinus    Hippotragus equinus Tropical and subtropical grasslands, savannas and shrublands -23.83333   31.50000     9 -1.60466631      -1.7980647  0.62643026      0.68961838
#2   3412        Canidae          Canis_latrans          Canis latrans                Temperate grasslands, savannas and shrublands  48.10611  -95.22806     9  2.76371560       1.6724629 -0.16317128      0.15757768
#3   3413        Canidae          Canis_latrans          Canis latrans                        Temperate broadleaf and mixed forests  47.95500  -94.78639     9 -0.07384233      -0.3553621 -2.08451767     -0.33689653
#4   3414        Canidae          Canis_latrans          Canis latrans                        Temperate broadleaf and mixed forests  47.31667  -92.54750     9 -0.62386433      -2.4498224  1.42849884     -0.21480807
#5   3415        Canidae          Canis_latrans          Canis latrans                Temperate grasslands, savannas and shrublands  45.33611  -94.40472     9 -0.22191731      -1.4640354 -0.12545840     -0.34475118
#6  18764      Leporidae       Lepus_americanus       Lepus americanus                                         Boreal forests/taiga  61.02759 -138.41045    34  0.45917941       0.2458739 -0.01279401     -0.02976037
#7  18887     Cricetidae Peromyscus_maniculatus Peromyscus maniculatus                                         Boreal forests/taiga  61.02759 -138.41045    34 -0.22803522      -0.1662930 -0.24643092     -0.24333106
#8  19103 Antilocapridae  Antilocapra_americana  Antilocapra americana                Temperate grasslands, savannas and shrublands  47.63095 -100.45888    34  0.66552194       0.8111799 -0.27100883     -0.33598661
#9  19783       Cervidae      Rangifer_tarandus      Rangifer tarandus                                         Boreal forests/taiga  48.62557  -56.42212    34 -0.09298552      -0.1636454  0.05520931      0.02840967
#10 22039       Cervidae Odocoileus_virginianus Odocoileus virginianus                        Temperate broadleaf and mixed forests  46.56531  -66.46191    34  0.94859184       0.9327810  0.10733196      0.01756898

















##********************************************************
## Look at some specific weather coefficients
## Exploration
##*******************************************************

plot(ssmresults$tempbetaSSM,ssmresults$`tempbetaGAMM[2]`)
abline(coef = c(0,1))

plot(ssmresults$precbetaSSM,ssmresults$`precbetaGAMM[2]`)
abline(coef = c(0,1))


#record 4: precipitation  - very different prec weather coeffcients betwenn GAMM and SSM model
##SSM beta= 1.42313603  --  GAMM beta = -0.21480807
record <- mammal_weather_CL[mammal_weather_CL$ID==IDs[4],]
library(ggplot2)
ggplot(record, aes(mean_precip_anomaly,pop_growth_rate))+geom_point()+stat_smooth(method = "lm", col = "red")

#loog at convergence and posterior distribution based on "ssmprectot" from record (ID==IDs[4])
par(mfrow = c(2,2))  ;   traceplot(ssmprectot)     
hist(ssmprectot$sims.list$beta)
hist(ssmprectot$sims.list$b0)

# Summarize posteriors
print(ssmprectot, digits=3)
#              mean     sd    2.5%    50%  97.5% overlap0     f  Rhat n.eff
#r[1]         0.567  0.373  -0.282  0.615  1.280     TRUE 0.934 1.000 49998
#r[2]        -0.206  0.307  -0.724 -0.248  0.453     TRUE 0.763 1.000 49998
#r[3]         0.440  0.289  -0.080  0.450  0.924     TRUE 0.952 1.000 15145
#r[4]         0.058  0.220  -0.380  0.042  0.527     TRUE 0.622 1.000 49998
#r[5]         0.187  0.206  -0.246  0.196  0.602     TRUE 0.862 1.000 23198
#r[6]         0.161  0.272  -0.361  0.128  0.753     TRUE 0.775 1.000 49998
#r[7]        -0.265  0.296  -0.799 -0.312  0.365     TRUE 0.822 1.000 49998
#r[8]         0.349  0.256  -0.200  0.376  0.825     TRUE 0.918 1.000 23915
#sigma2.obs   0.095  0.128   0.000  0.052  0.468    FALSE 1.000 1.000 15120
#sigma2.proc  0.162  0.174   0.001  0.108  0.683    FALSE 1.000 1.000 15866
#N[1]         7.872  2.459   4.515  7.369 13.898    FALSE 1.000 1.000 49998
#N[2]        13.874  4.090   7.322 13.702 23.186    FALSE 1.000 1.000 49998
#N[3]        11.173  2.772   7.504 10.565 17.768    FALSE 1.000 1.000 23789
#N[4]        17.233  3.436  10.947 17.339 24.386    FALSE 1.000 1.000 49998
#N[5]        18.308  3.921  11.303 18.394 26.718    FALSE 1.000 1.000 36898
#N[6]        22.163  4.993  13.049 22.428 31.999    FALSE 1.000 1.000 49998
#N[7]        25.972  5.763  16.061 25.754 38.978    FALSE 1.000 1.000 17141
#N[8]        20.085  5.504  12.911 18.626 33.463    FALSE 1.000 1.001 10760
#N[9]        28.434  7.549  16.844 27.343 46.515    FALSE 1.000 1.001 49998
#b0           0.137  0.150  -0.180  0.139  0.448     TRUE 0.860 1.000 26038
#beta         1.409  1.476  -1.676  1.443  4.315     TRUE 0.857 1.000 49998 #large sd 
#deviance    -6.070 17.404 -52.497 -1.417 15.724     TRUE 0.547 1.001  5359





#record 3: precipitation SSM beta= -2.10564141     GAMM beta = -0.33689653
record <- mammal_weather_CL[mammal_weather_CL$ID==IDs[3],]
ggplot(record, aes(mean_temp_anomaly,pop_growth_rate))+geom_point()+stat_smooth(method = "lm", col = "red")

ggplot(record, aes(mean_precip_anomaly,pop_growth_rate))+geom_point()+stat_smooth(method = "lm", col = "red")
par(mfrow = c(2,2))  ;   traceplot(ssmprectot)     
hist(ssmprectot$sims.list$beta)

# Summarize posteriors
print(ssmprectot, digits=3)
