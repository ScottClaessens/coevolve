
<!-- README.md is generated from README.Rmd. Please edit that file -->

<img src="man/figures/logo.png" align="right" height="139" alt="" />

# coevolve

<!-- badges: start -->
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
`bernoulli_logit`, `ordered_logistic`, and `poisson_softmax`
respectively.

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
      z = "poisson_softmax"
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
#> Chain 3 finished in 410.7 seconds.
#> Chain 4 finished in 413.0 seconds.
#> Chain 1 finished in 413.9 seconds.
#> Chain 2 finished in 416.1 seconds.
#> 
#> All 4 chains finished successfully.
#> Mean chain execution time: 413.4 seconds.
#> Total execution time: 416.4 seconds.
```

The results can be investigated using:

``` r
summary(m)
#> Variables: x = bernoulli_logit 
#>            y = ordered_logistic 
#>            z = poisson_softmax 
#>      Data: d (Number of observations: 10)
#>     Draws: 4 chains, each with iter = 1000; warmup = 1000; thin = 1
#>            total post-warmup draws = 4000
#> 
#> Autoregressive selection effects:
#>   Estimate Est.Error  2.5% 97.5% Rhat Bulk_ESS Tail_ESS
#> x     0.05      1.01 -2.00  2.00 1.00     6749     3096
#> y     0.00      1.00 -1.94  1.96 1.00     6782     3359
#> z    -0.72      0.75 -2.32  0.63 1.00     7058     2913
#> 
#> Cross selection effects:
#>       Estimate Est.Error  2.5% 97.5% Rhat Bulk_ESS Tail_ESS
#> x ⟶ y     0.02      0.97 -1.92  1.90 1.00     7368     3142
#> x ⟶ z    -0.36      1.07 -2.36  1.84 1.00     4109     3381
#> y ⟶ x     0.01      1.00 -1.94  1.94 1.00     5984     2771
#> y ⟶ z    -0.20      1.10 -2.30  1.99 1.00     3267     3170
#> z ⟶ x    -0.28      0.84 -1.90  1.39 1.00     5231     3161
#> z ⟶ y    -0.19      0.91 -1.92  1.61 1.00     4562     3073
#> 
#> Drift scale parameters:
#>   Estimate Est.Error 2.5% 97.5% Rhat Bulk_ESS Tail_ESS
#> x     0.92      0.66 0.04  2.42 1.00     3005     2235
#> y     0.84      0.62 0.03  2.32 1.00     3447     2448
#> z     0.87      0.60 0.05  2.17 1.00     3078     2282
#> 
#> Continuous time intercept parameters:
#>   Estimate Est.Error  2.5% 97.5% Rhat Bulk_ESS Tail_ESS
#> x    -0.28      0.98 -2.21  1.64 1.00     7253     2935
#> y    -0.04      0.98 -1.89  1.92 1.00     6932     2545
#> z     0.92      0.90 -0.90  2.65 1.00     5979     3059
#> 
#> Ordinal cutpoint parameters:
#>      Estimate Est.Error  2.5% 97.5% Rhat Bulk_ESS Tail_ESS
#> y[1]    -1.81      1.17 -4.13  0.43 1.00     2894     2681
#> y[2]     1.83      1.17 -0.43  4.16 1.00     4462     3616
#> Warning: There were 3 divergent transitions after warmup.
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
