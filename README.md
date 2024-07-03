
<!-- README.md is generated from README.Rmd. Please edit that file -->

<img src="man/figures/logo.png" height="139" alt="coevolve Logo"/>[<img src="https://raw.githubusercontent.com/stan-dev/logos/master/logo_tm.png" align="right" height="139" alt="Stan Logo"/>](https://mc-stan.org/)

# coevolve

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

Then, simulate data for an ordinal trait and a count trait.

``` r
# simulate data
d <- 
  data.frame(
    # id to match dataset to tree tips
    id = tree$tip.label,
    # simulate variables
    x = ordered(sample(1:3, size = n, replace = TRUE)),
    y = rpois(n, 2)
  )

head(d)
#>   id x y
#> 1 t6 1 3
#> 2 t4 3 0
#> 3 t9 2 4
#> 4 t8 2 1
#> 5 t2 3 3
#> 6 t7 3 1
```

We can then fit our Bayesian dynamic coevolutionary model in `cmdstanr`
with the `coev_fit()` function. We declare both variables and set the
response distributions for ordinal and count variables as
`ordered_logistic` and `poisson_softplus` respectively.

``` r
# load the coevolve package
library(coevolve)

# fit model
m <-
  coev_fit(
    data = d,
    variables = list(
      x = "ordered_logistic",
      y = "poisson_softplus"
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
#> Chain 1 finished in 33.2 seconds.
#> Chain 4 finished in 33.4 seconds.
#> Chain 2 finished in 33.9 seconds.
#> Chain 3 finished in 39.9 seconds.
#> 
#> All 4 chains finished successfully.
#> Mean chain execution time: 35.1 seconds.
#> Total execution time: 40.2 seconds.
```

The results can be investigated using:

``` r
summary(m)
#> Variables: x = ordered_logistic 
#>            y = poisson_softplus 
#>      Data: d (Number of observations: 10)
#>     Draws: 4 chains, each with iter = 1000; warmup = 1000; thin = 1
#>            total post-warmup draws = 4000
#> 
#> Autoregressive selection effects:
#>   Estimate Est.Error  2.5% 97.5% Rhat Bulk_ESS Tail_ESS
#> x    -0.84      0.63 -2.39 -0.03 1.00     2500     1302
#> y    -0.53      0.43 -1.58 -0.02 1.00     2653     1555
#> 
#> Cross selection effects:
#>       Estimate Est.Error  2.5% 97.5% Rhat Bulk_ESS Tail_ESS
#> x ⟶ y     0.12      1.10 -2.08  2.21 1.00     2173     3174
#> y ⟶ x     0.08      0.87 -1.69  1.75 1.00     3149     2635
#> 
#> Drift scale parameters:
#>   Estimate Est.Error 2.5% 97.5% Rhat Bulk_ESS Tail_ESS
#> x     0.74      0.56 0.03  2.09 1.00     3592     1951
#> y     0.74      0.58 0.03  2.13 1.00     2612     1631
#> 
#> Continuous time intercept parameters:
#>   Estimate Est.Error  2.5% 97.5% Rhat Bulk_ESS Tail_ESS
#> x     0.12      0.96 -1.82  2.01 1.00     3300     3164
#> y     1.05      0.86 -0.64  2.69 1.00     3973     2885
#> 
#> Ordinal cutpoint parameters:
#>      Estimate Est.Error  2.5% 97.5% Rhat Bulk_ESS Tail_ESS
#> x[1]    -1.87      1.20 -4.34  0.40 1.00     2076     2780
#> x[2]     1.22      1.15 -0.96  3.57 1.00     2599     3161
#> Warning: There were 20 divergent transitions after warmup.
#> http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
```

We can infer whether these traits influence one another in their
evolution by plotting $\Delta\theta_{z}$, the expected change in the
equilibrium trait value of one variable from an absolute deviation
increase in another variable.

``` r
coev_plot_delta_theta(m)
```

<img src="man/figures/README-plot_delta_theta-1.png" width="60%" style="display: block; margin: auto;" />

This plot suggests that these variables do not coevolve over
evolutionary time.

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
