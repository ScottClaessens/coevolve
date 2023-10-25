
<!-- README.md is generated from README.Rmd. Please edit that file -->

# coevolve

<!-- badges: start -->
<!-- badges: end -->

## Overview

The **coevolve** package allows the user to fit Bayesian dynamic
coevolutionary models in Stan.

## Installation

To use the **coevolve** package, you must first install the `cmdstanr`
package (see full installation instructions here:
<https://mc-stan.org/cmdstanr/>).

``` r
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/",
                 getOption("repos")))
```

You can then install the development version of `coevolve` with:

``` r
# install.packages("devtools")
devtools::install_github("ScottClaessens/coevolve")
```

## Example

We can simulate a phylogenetic tree with some data to see how the
package works. First, simulate a phylogenetic tree.

``` r
# set the random seed
set.seed(1)
# number of taxa
n <- 20
# random tree
tree <- ape::rtree(n)
```

Then, simulate data for a binary trait, an ordinal trait, and a count
trait.

``` r
# simulate data
d <- 
  data.frame(
    # id to match dataset to tree tips
    id = tree$tip.label,
    # simulate variables
    x = rbinom(n, 1, 0.5),
    y = ordered(sample(1:3, size = n, replace = TRUE)),
    z = rpois(n, 2)
  )

head(d)
#>    id x y z
#> 1 t10 1 1 2
#> 2 t14 1 1 2
#> 3 t20 0 2 1
#> 4  t7 1 1 3
#> 5  t9 1 1 2
#> 6 t15 0 1 2
```

We can then fit our Bayesian dynamic coevolutionary model in `cmdstanr`
with the `coev_fit()` function. We declare all variables and set the
response distributions for binary, ordinal, and count variables as
`bernoulli_logit`, `ordered_logistic`, and `poisson_log` respectively.

``` r
# load the coevolve package
library(coevolve)

# fit model
m <-
  coev_fit(
    data = d,
    variables = list(
      x = "bernoulli_logit",
      y = "ordered_logistic",
      z = "poisson_log"
    ),
    id = "id",
    tree = tree,
    # additional arguments for cmdstanr
    parallel_chains = 4,
    refresh = 0,
    show_messages = FALSE,
    seed = 1
  )
#> Running MCMC with 4 parallel chains...
#> 
#> Chain 4 finished in 1091.8 seconds.
#> Chain 1 finished in 1093.4 seconds.
#> Chain 3 finished in 1126.8 seconds.
#> Chain 2 finished in 1244.7 seconds.
#> 
#> All 4 chains finished successfully.
#> Mean chain execution time: 1139.2 seconds.
#> Total execution time: 1245.3 seconds.
#> Warning: 14 of 4000 (0.0%) transitions ended with a divergence.
#> See https://mc-stan.org/misc/warnings for details.

summary(m)
#> Variables: x = bernoulli_logit 
#>            y = ordered_logistic 
#>            z = poisson_log 
#>      Data: d (Number of observations: 20)
#>     Draws: 4 chains, each with iter = 1000; warmup = 1000; thin = 1
#>            total post-warmup draws = 4000
#> 
#> Autoregressive selection effects:
#>    mean median   sd  mad    q5  q95 rhat ess_bulk ess_tail
#> x  0.38   0.42 1.03 1.07 -1.33 1.98 1.00     4632     2943
#> y -0.07  -0.06 0.95 0.91 -1.66 1.47 1.00     3824     2972
#> z  0.01   0.07 0.84 0.83 -1.49 1.27 1.00     3621     3038
#> 
#> Cross selection effects:
#>        mean median   sd  mad    q5  q95 rhat ess_bulk ess_tail
#> x ⟶ y -0.07  -0.07 1.03 1.01 -1.74 1.64 1.00     5053     3145
#> x ⟶ z  0.05   0.05 0.97 0.95 -1.53 1.61 1.00     3279     2068
#> y ⟶ x -0.05  -0.06 0.95 0.93 -1.63 1.49 1.00     4514     2684
#> y ⟶ z  0.08   0.09 0.92 0.95 -1.42 1.53 1.00     3150     3068
#> z ⟶ x -0.09  -0.10 0.95 0.93 -1.60 1.52 1.00     4559     2250
#> z ⟶ y  0.28   0.28 0.97 0.98 -1.35 1.89 1.00     4933     3129
#> 
#> Drift scale parameters:
#>   mean median   sd  mad   q5  q95 rhat ess_bulk ess_tail
#> x 0.76   0.66 0.57 0.56 0.06 1.84 1.00     2520     1812
#> y 0.86   0.75 0.62 0.64 0.07 2.02 1.01     1944     1740
#> z 0.51   0.44 0.37 0.37 0.04 1.22 1.00     1923     1769
#> 
#> Note: Not all model parameters are displayed in this summary.
#> Warning: There were 14 divergent transitions after warmup.
#> http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
```

From the cross selection effects, we infer that the three traits do not
influence one another in their evolution.

We can also plot the cross selection effects from this model using the
`coev_plot_cross()` function.

``` r
coev_plot_cross(m)
```

<img src="man/figures/README-plot_cross-1.png" width="60%" style="display: block; margin: auto;" />

## Citing coevolve

When using the coevolve package, please cite the following papers:

- Ringen, E., Martin, J. S., & Jaeggi, A. (2021). Novel phylogenetic
  methods reveal that resource-use intensification drives the evolution
  of “complex” societies. *EcoEvoRXiv*.
  <https://doi.org/10.32942/osf.io/wfp95>
- Sheehan, O., Watts, J., Gray, R. D., Bulbulia, J., Claessens, S.,
  Ringen, E. J., & Atkinson, Q. D. (2023). Coevolution of religious and
  political authority in Austronesian societies. *Nature Human
  Behaviour*, *7*(1), 38-45.
  <https://doi.org/10.1038/s41562-022-01471-y>
