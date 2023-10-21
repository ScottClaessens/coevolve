
<!-- README.md is generated from README.Rmd. Please edit that file -->

# coevolve

<!-- badges: start -->
<!-- badges: end -->

## Overview

The **coevolve** package allows the user to fit Bayesian dynamic
coevolutionary models in Stan.

## Installation

You can install the development version of coevolve from
[GitHub](https://github.com/) with:

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
n <- 50
# random tree
tree <- ape::rtree(n)
```

Then, simulate data for two binary traits. We allow these traits to
evolve independently on the tree, but they do not influence one another
in their evolution. In other words, they do not coevolve.

``` r
# simulate data
d <- 
  data.frame(
    # id to match dataset to tree tips
    id = tree$tip.label,
    # simulate binary variables (0/1 integers)
    # that evolve independently on the tree
    x = as.integer(ape::rTraitDisc(tree)) - 1L,
    y = as.integer(ape::rTraitDisc(tree)) - 1L
  )

head(d)
#>    id x y
#> 1 t42 0 1
#> 2 t38 0 1
#> 3 t47 0 1
#> 4 t20 0 1
#> 5 t28 0 0
#> 6 t48 0 0
```

We can then fit our Bayesian dynamic coevolutionary model in `cmdstanr`
with the `coev_fit()` function. We declare both variables and set the
response distribution `bernoulli_logit`.

``` r
# load the coevolve package
library(coevolve)

# fit model
m <-
  coev_fit(
    data = d,
    variables = list(
      x = "bernoulli_logit",
      y = "bernoulli_logit"
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
#> Chain 2 finished in 422.8 seconds.
#> Chain 1 finished in 475.2 seconds.
#> Chain 3 finished in 490.7 seconds.
#> Chain 4 finished in 492.8 seconds.
#> 
#> All 4 chains finished successfully.
#> Mean chain execution time: 470.3 seconds.
#> Total execution time: 493.1 seconds.
#> Warning: 2 of 4000 (0.0%) transitions ended with a divergence.
#> See https://mc-stan.org/misc/warnings for details.

m
#> Variables: x = bernoulli_logit 
#>            y = bernoulli_logit 
#>      Data: d (Number of observations: 50)
#>     Draws: 4 chains, each with iter = 1000; warmup = 1000; thin = 1
#>            total post-warmup draws = 4000
#> 
#> Autoregressive selection effects:
#>    mean median   sd  mad    q5   q95 rhat ess_bulk ess_tail
#> x -1.17  -1.14 0.69 0.68 -2.37 -0.12 1.00  6171.48  2752.86
#> y -0.70  -0.68 0.82 0.80 -2.06  0.60 1.00  4817.27  2714.89
#> 
#> Cross selection effects:
#>       mean median   sd  mad    q5  q95 rhat ess_bulk ess_tail
#> x ⟶ y 0.78   0.77 0.62 0.59 -0.20 1.82 1.00  4309.08  2677.62
#> y ⟶ x 1.26   1.25 0.88 0.86 -0.17 2.70 1.00  3584.96  2557.74
#> 
#> Drift scale parameters:
#>   mean median   sd  mad   q5  q95 rhat ess_bulk ess_tail
#> x 0.71   0.59 0.54 0.51 0.06 1.76 1.00  2575.25  2194.75
#> y 2.04   2.03 0.79 0.76 0.68 3.36 1.00  2149.39  1495.79
```

From the cross selection effects, we correctly infer that the two traits
do not influence one another in their evolution.

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
