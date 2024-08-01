
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
#> Chain 4 finished in 4944.0 seconds.
#> Chain 3 finished in 5209.1 seconds.
#> Chain 2 finished in 5267.9 seconds.
#> Chain 1 finished in 5364.4 seconds.
#> 
#> All 4 chains finished successfully.
#> Mean chain execution time: 5196.4 seconds.
#> Total execution time: 5364.8 seconds.
#> Warning: 14 of 2000 (1.0%) transitions ended with a divergence.
#> See https://mc-stan.org/misc/warnings for details.
```

The results can be investigated using:

``` r
summary(fit)
#> Variables: political_authority = ordered_logistic 
#>            religious_authority = ordered_logistic 
#>      Data: authority$data (Number of observations: 97)
#>     Draws: 4 chains, each with iter = 500; warmup = 1000; thin = 1
#>            total post-warmup draws = 2000
#> 
#> Autoregressive selection effects:
#>                     Estimate Est.Error  2.5% 97.5% Rhat Bulk_ESS Tail_ESS
#> political_authority    -0.64      0.51 -1.90 -0.02 1.00     1720     1149
#> religious_authority    -0.67      0.52 -1.97 -0.03 1.00     1334     1183
#> 
#> Cross selection effects:
#>                                           Estimate Est.Error 2.5% 97.5% Rhat
#> political_authority ⟶ religious_authority     2.85      1.08 0.82  5.12 1.00
#> religious_authority ⟶ political_authority     2.56      1.06 0.54  4.73 1.00
#>                                           Bulk_ESS Tail_ESS
#> political_authority ⟶ religious_authority      780      889
#> religious_authority ⟶ political_authority      846     1032
#> 
#> Drift scale parameters:
#>                     Estimate Est.Error 2.5% 97.5% Rhat Bulk_ESS Tail_ESS
#> political_authority     1.47      0.78 0.11  3.14 1.00      811      733
#> religious_authority     1.20      0.74 0.06  2.81 1.00      675      661
#> 
#> Continuous time intercept parameters:
#>                     Estimate Est.Error  2.5% 97.5% Rhat Bulk_ESS Tail_ESS
#> political_authority     0.11      0.95 -1.74  2.00 1.01     3925     1345
#> religious_authority     0.14      0.94 -1.74  2.00 1.00     3321     1421
#> 
#> Ordinal cutpoint parameters:
#>                        Estimate Est.Error  2.5% 97.5% Rhat Bulk_ESS Tail_ESS
#> political_authority[1]    -1.00      0.85 -2.61  0.66 1.00     2427     1263
#> political_authority[2]    -0.31      0.84 -1.88  1.32 1.00     2406     1263
#> political_authority[3]     1.64      0.86 -0.01  3.30 1.00     2171     1532
#> religious_authority[1]    -1.33      0.88 -3.03  0.39 1.00     2391     1467
#> religious_authority[2]    -0.69      0.88 -2.44  1.03 1.00     2508     1317
#> religious_authority[3]     1.59      0.91 -0.15  3.42 1.00     2336     1398
#> Warning: There were 14 divergent transitions after warmup.
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
coev_plot_delta_theta(fit)
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
