
<!-- README.md is generated from README.Rmd. Please edit that file -->

<img src="man/figures/logo.png" height="139" alt="coevolve Logo"/>[<img src="https://raw.githubusercontent.com/stan-dev/logos/master/logo_tm.png" align="right" height="139" alt="Stan Logo"/>](https://mc-stan.org/)

# coevolve

<!-- badges: start -->

[![R-CMD-check](https://github.com/ScottClaessens/coevolve/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ScottClaessens/coevolve/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Overview

The **coevolve** package allows the user to fit Bayesian dynamic
coevolutionary phylogenetic models in Stan. These models can be used to
estimate how variables have coevolved over evolutionary time and to
assess causal directionality (X → Y vs. Y → X) and contingencies (X,
then Y) in evolution.

While existing methods only allow pairs of binary traits to coevolve
(e.g.,
[BayesTraits](https://www.evolution.reading.ac.uk/BayesTraitsV4.1.2/BayesTraitsV4.1.2.html)),
the **coevolve** package allows users to include multiple traits of
different data types, including binary, ordinal, count, and continuous
traits.

## Installation

To use the **coevolve** package, you must first install the `cmdstanr`
package (see full installation instructions here:
<https://mc-stan.org/cmdstanr/>).

``` r
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
```

You can then install the development version of **coevolve** with:

``` r
# install.packages("devtools")
devtools::install_github("ScottClaessens/coevolve")
```

## How to use coevolve

``` r
library(coevolve)
```

As an example, we analyse the coevolution of political and religious
authority in 97 Austronesian societies. These data were compiled and
analysed in [Sheehan et
al. (2023)](https://www.nature.com/articles/s41562-022-01471-y). Both
variables are four-level ordinal variables reflecting increasing levels
of authority. We use a phylogeny of Austronesian languages to assess
patterns of coevolution.

``` r
fit <-
  coev_fit(
    data = authority$data,
    variables = list(
      political_authority = "ordered_logistic",
      religious_authority = "ordered_logistic"
    ),
    id = "language",
    tree = authority$phylogeny,
    # manually set prior
    prior = list(A_offdiag = "normal(0, 2)"),
    # arguments for cmdstanr
    parallel_chains = 4,
    iter_sampling = 500,
    refresh = 0,
    seed = 1
  )
#> Running MCMC with 4 parallel chains...
#> 
#> Chain 3 finished in 910.2 seconds.
#> Chain 1 finished in 910.8 seconds.
#> Chain 4 finished in 917.5 seconds.
#> Chain 2 finished in 955.5 seconds.
#> 
#> All 4 chains finished successfully.
#> Mean chain execution time: 923.5 seconds.
#> Total execution time: 955.7 seconds.
#> Warning: 30 of 2000 (2.0%) transitions ended with a divergence.
#> See https://mc-stan.org/misc/warnings for details.
```

The results can be investigated using:

``` r
summary(fit)
#> Variables: political_authority = ordered_logistic 
#>            religious_authority = ordered_logistic 
#>      Data: authority$data (Number of observations: 97)
#> Phylogeny: authority$phylogeny (Number of trees: 1)
#>     Draws: 4 chains, each with iter = 500; warmup = 1000; thin = 1
#>            total post-warmup draws = 2000
#> 
#> Autoregressive selection effects:
#>                     Estimate Est.Error  2.5% 97.5% Rhat Bulk_ESS Tail_ESS
#> political_authority    -0.66      0.51 -1.90 -0.03 1.01      590      996
#> religious_authority    -0.80      0.59 -2.15 -0.03 1.00     1315     1007
#> 
#> Cross selection effects:
#>                                           Estimate Est.Error  2.5% 97.5% Rhat
#> political_authority ⟶ religious_authority     2.28      1.03  0.28  4.36 1.01
#> religious_authority ⟶ political_authority     1.76      1.10 -0.37  4.02 1.00
#>                                           Bulk_ESS Tail_ESS
#> political_authority ⟶ religious_authority      543      327
#> religious_authority ⟶ political_authority      500      801
#> 
#> Drift parameters:
#>                                              Estimate Est.Error  2.5% 97.5%
#> sd(political_authority)                          2.03      0.82  0.35  3.59
#> sd(religious_authority)                          1.26      0.78  0.06  2.85
#> cor(political_authority,religious_authority)     0.27      0.31 -0.39  0.76
#>                                              Rhat Bulk_ESS Tail_ESS
#> sd(political_authority)                      1.00      416      475
#> sd(religious_authority)                      1.01      313      779
#> cor(political_authority,religious_authority) 1.00     1235     1378
#> 
#> Continuous time intercept parameters:
#>                     Estimate Est.Error  2.5% 97.5% Rhat Bulk_ESS Tail_ESS
#> political_authority     0.21      0.95 -1.65  2.08 1.00     2042     1491
#> religious_authority     0.20      0.94 -1.55  2.06 1.00     2337     1365
#> 
#> Ordinal cutpoint parameters:
#>                        Estimate Est.Error  2.5% 97.5% Rhat Bulk_ESS Tail_ESS
#> political_authority[1]    -1.33      0.89 -3.04  0.53 1.00     1034      241
#> political_authority[2]    -0.56      0.87 -2.19  1.35 1.00     1023      241
#> political_authority[3]     1.66      0.89 -0.04  3.44 1.00      979      253
#> religious_authority[1]    -1.49      0.98 -3.42  0.54 1.00      739      226
#> religious_authority[2]    -0.80      0.95 -2.64  1.19 1.00      752      206
#> religious_authority[3]     1.62      0.98 -0.29  3.71 1.00     1049      209
#> Warning: There were 30 divergent transitions after warmup.
#> http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
```

The summary provides general information about the model and details on
the posterior draws for the model parameters. In particular, the output
shows the autoregressive selection effects (i.e., the effect of a
variable on itself in the future), the cross selection effects (i.e.,
the effect of a variable on another variable in the future), the amount
of drift, continuous time intercept parameters for the schocastic
differential equation, and cutpoints for the ordinal variables.

While this summary output is useful as a first glance, it is difficult
to interpret these parameters directly to infer directions of
coevolution. Another approach is to “intervene” in the system. We can
hold variables of interest at their average values and then increase one
variable by a standardised amount to see how this affects the optimal
trait value for another variable.

The `coev_plot_delta_theta()` function allows us to visualise
$\Delta\theta_{z}$ for all variable pairs in the model.
$\Delta\theta_{z}$ is defined as the change in the optimal trait value
of one variable which results from a one median absolute deviation
increase in another variable.

``` r
coev_plot_delta_theta(fit, limits = c(-5, 50))
#> Warning: Removed 119 rows containing non-finite outside the scale range
#> (`stat_density()`).
```

<img src="man/figures/README-authority-delta-theta-1.png" width="60%" style="display: block; margin: auto;" />

This plot suggests that both variables influence one another in their
coevolution. A standardised increase in political authority results in
an increase in the optimal trait value for religious authority, and vice
versa. In other words, these two variables reciprocally coevolve over
evolutionary time.

## Citing coevolve

When using the **coevolve** package, please cite the following papers:

- Ringen, E., Martin, J. S., & Jaeggi, A. (2021). Novel phylogenetic
  methods reveal that resource-use intensification drives the evolution
  of “complex” societies. *EcoEvoRXiv*.
  <https://doi.org/10.32942/osf.io/wfp95>
- Sheehan, O., Watts, J., Gray, R. D., Bulbulia, J., Claessens, S.,
  Ringen, E. J., & Atkinson, Q. D. (2023). Coevolution of religious and
  political authority in Austronesian societies. *Nature Human
  Behaviour*, *7*(1), 38-45.
  <https://doi.org/10.1038/s41562-022-01471-y>
