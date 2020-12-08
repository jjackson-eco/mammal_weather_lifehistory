####################################################
##                                                ##
##     Global climate and population dynamics     ##
##                                                ##
##       Density dependence simulation + env      ##
##                                                ##
##                Dec 8th 2020                    ##
##                                                ##
####################################################

## Adapted from advice from Dylan Z Childs

rm(list=ls())
options(width = 100)

### Setting up the initial parameters
t_step = 5

p <- list()
p$b0 <- 1        # Initial value to feed to the density dependence
p$bDD <- 0.5     # coefficient for density dependence
p$bENV <- 0.2    # coefficient for the environment
p$sigma <- 1     # noise variance
p$sigmaENV <- 1  # environmental variance

x0 <- 2.0        # starting abundance

n_sims <- 1000

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

## Function to run the simulation and pull out the necessary work
run_sim <- function(p, x0, t_step, n_sims, beta) {
  slopes <- numeric(n_sims)
  for (i in seq_len(n_sims)){
    x <- iter(p, x0, t_step)
    slopes[i] <- fit_mod(p, x, beta = beta)
  }
  slopes
}
 
## Simulations to extract the weather effect
set.seed(666)
jpeg(filename = "plots/annual_abundance/Accounting_for_autocorrelation_testing/Environment_coefficient_simulation.jpeg",
     width = 20, height = 13, units = "cm",res = 400)
par(mfrow = c(2,3))
for(i in c(5,10,20,50,100,1000)){
  csim = run_sim(p, x0, t_step = i, n_sims = n_sims, beta = 3)
  mn_sim = mean(csim)
  hist(csim, xlab = "Environment coefficient",
       main = paste0(i, " years: mean diff = ", 
                     round(p$bENV - mn_sim, 3)))
  abline(v = p$bENV, lty = 1, lwd = 2,
         col = "darkgreen")
  abline(v = mn_sim, lty = 2, lwd = 2, col = "green")
}
dev.off()

## Simulations to extract the density dependence effect
set.seed(666)
jpeg(filename = "plots/annual_abundance/Accounting_for_autocorrelation_testing/DD_coefficient_simulation.jpeg",
     width = 20, height = 13, units = "cm",res = 400)
par(mfrow = c(2,3))
for(i in c(5,10,20,50,100,1000)){
  csim = run_sim(p, x0, t_step = i, n_sims = n_sims, beta = 2)
  mn_sim = mean(csim)
  hist(csim, xlab = "Density dependence coefficient",
       main = paste0(i, " years: mean diff = ", 
                     round(-p$bDD - mn_sim, 3)))
  abline(v = -p$bDD, lty = 1, lwd = 2, 
         col = "darkred")
  abline(v = mn_sim, lty = 2, lwd = 2, col = "red")
}
dev.off()


