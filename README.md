
<!-- README.md is generated from README.Rmd. Please edit that file -->

# coevolve <img src="man/figures/logo.png" align="right" height="139" alt="" />

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
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
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
n <- 10
# random tree
tree <- ape::rcoal(n)
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
#>   id x y z
#> 1 t6 0 2 2
#> 2 t4 0 2 4
#> 3 t9 0 2 4
#> 4 t8 0 1 1
#> 5 t2 0 2 3
#> 6 t7 1 2 5
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
    seed = 1
  )
#> Running MCMC with 4 parallel chains...
#> 
#> Chain 1 finished in 572.3 seconds.
#> Chain 2 finished in 578.9 seconds.
#> Chain 4 finished in 653.3 seconds.
#> Chain 3 finished in 679.8 seconds.
#> 
#> All 4 chains finished successfully.
#> Mean chain execution time: 621.1 seconds.
#> Total execution time: 680.2 seconds.
```

The results can be investigated using:

``` r
summary(m)
#> Variables: x = bernoulli_logit 
#>            y = ordered_logistic 
#>            z = poisson_log 
#>      Data: d (Number of observations: 10)
#>     Draws: 4 chains, each with iter = 1000; warmup = 1000; thin = 1
#>            total post-warmup draws = 4000
#> 
#> Autoregressive selection effects:
#>   Estimate Est.Error  2.5% 97.5% Rhat Bulk_ESS Tail_ESS
#> x     0.07      0.97 -1.85  1.98 1.00     3924     2823
#> y     0.05      1.01 -1.95  1.99 1.00     3875     2727
#> z    -0.10      0.83 -1.84  1.38 1.00     3002     2731
#> 
#> Cross selection effects:
#>       Estimate Est.Error  2.5% 97.5% Rhat Bulk_ESS Tail_ESS
#> x ⟶ y     0.01      1.03 -2.04  1.96 1.00     4285     2897
#> x ⟶ z    -0.18      1.00 -2.14  1.74 1.00     1847     1295
#> y ⟶ x     0.05      0.99 -1.88  2.04 1.00     3295     3130
#> y ⟶ z    -0.12      0.98 -2.00  1.87 1.00     2516     2003
#> z ⟶ x    -0.17      0.95 -2.00  1.69 1.00     3166     1419
#> z ⟶ y    -0.08      0.97 -1.91  1.88 1.00     3223     2365
#> 
#> Drift scale parameters:
#>   Estimate Est.Error 2.5% 97.5% Rhat Bulk_ESS Tail_ESS
#> x     0.89      0.64 0.05  2.38 1.00     1877     2009
#> y     0.78      0.60 0.03  2.20 1.00     1909     1757
#> z     0.58      0.44 0.02  1.68 1.00     1759     1492
#> 
#> Continuous time intercept parameters:
#>   Estimate Est.Error  2.5% 97.5% Rhat Bulk_ESS Tail_ESS
#> x    -0.33      0.91 -2.10  1.47 1.00     4323     2973
#> y    -0.06      0.94 -1.88  1.77 1.00     4097     3033
#> z     0.62      0.82 -1.04  2.20 1.00     3166     2554
#> 
#> Ordinal cutpoint parameters:
#>      Estimate Est.Error  2.5% 97.5% Rhat Bulk_ESS Tail_ESS
#> y[1]    -1.66      1.04 -3.79  0.25 1.00     2245     2687
#> y[2]     1.95      1.11 -0.14  4.27 1.00     4332     3213
#> Warning: There were 15 divergent transitions after warmup.
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
