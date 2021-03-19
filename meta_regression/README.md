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

Thus, in all subsequent Bayesian models in model selection, we used regularising priors, i.e. normal priors with standard deviation < 1 and exponential priors with rates > 5. Please see meta-regression scripts for specific details on each prior.

## 2. General form of the Gaussian meta-regression
<details>
  <summary>Click here to expand</summary>

### `phylo_temp_GAM.R`
### `phylo_precip_GAM.R`

Now we are going to present the general framework that was used to fit the meta-regressions in models and that were used to explore consistent patterns across the mammals using Guassian models. In the general models that were explored for consistent patterns we incorporated species variance and also phylogenetic covariance. However, we also explored spatial autocorrelation (see section 4).

```

```

</details>

## 3. General form of the Gamma meta-regression
<details>
  <summary>Click here to expand</summary>

</details>

## 4. Spatial autocorrelation
<details>
  <summary>Click here to expand</summary>

</details>
