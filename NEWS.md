# coevolve - development version

### Bug Fixes

* Fixed `coev_plot_predictive_check()` when `nuts_sampler = "nutpie"`. The
  JAX backend now returns `yrep` posterior predictions alongside the other
  model parameters, matching the Stan backend (#118)
* `log_lik = TRUE` is now supported with `nuts_sampler = "nutpie"`. The JAX
  backend previously ignored the argument and returned no `log_lik`, so
  `loo` and `waic` could not be computed from nutpie fits (#118)
* Fixed `prior_only = TRUE` for models with gaussian variables and no
  repeated observations. `terminal_drift` only received a prior inside the
  likelihood block, leaving it improper when the likelihood was skipped.
  Sampling failed to initialise with `nuts_sampler = "nutpie"` and was
  unreliable with `nuts_sampler = "stan"` (#118)
* Fixed the pointwise `log_lik` and `yrep` for gaussian variables with
  missing data. The generated quantities block used the `-9999` missing
  value placeholder when constructing residuals, which corrupted the
  conditional densities of the other variables of the same observation

# coevolve 1.2.0

### New Features

* Added `coev_ancestral_states()` for extracting posterior estimates of
  ancestral trait values at internal phylogenetic nodes, on either the
  latent or response scale (#86)
* Added vignette "Ancestral state reconstruction"

### Bug Fixes

* Fixed issue with plotting functions when `nuts_sampler = "nutpie"` (#114)

### Other Changes

* Updated documentation (#111)
* Fixed cmdstanr location in documentation (#112)
* Responded to [rOpenSci](https://ropensci.org/) package review (#113)
* Added package citation (#116)

# coevolve 1.1.0

### New Features

* Allowed for single traits (#107)
* Added pure JAX/NumPyro backend via `nuts_sampler = "nutpie"`.
  Uses nutpie's Rust NUTS sampler with JAX gradients for ~5x
  faster sampling than Stan on typical models. Requires
  `pip install jax numpyro nutpie` (#109)

# coevolve 1.0.0

### New Features

* Added `nutpie` as an alternative sampler for the Stan models (#92)
* Cached matrix computations when branch lengths are identical (#93)
* Implemented Hilbert-space approximate Gaussian processes for spatial control,
  adding the `lon_lat` argument and deprecating the `dist_mat` argument in
  `coev_fit()` (#103)

### Bug Fixes

* Fixed issue with `summary()` when `estimate_residual = FALSE` (#95)

### Other Changes

* Updated license
* Submitted to [rOpenSci](https://ropensci.org/) for package review
* Linted package using [lintr](https://lintr.r-lib.org/)
* Implemented [srr](https://docs.ropensci.org/srr/) compliance checks
* Added continuous integration checks using GitHub Actions
* Reduced cyclomatic complexity for some functions
* Re-factored Stan code generation to use 
  [whisker](https://github.com/edwindj/whisker) templates (#99)
* Implemented automatic test fixture re-generation in GitHub Actions (#100)

# coevolve 0.1.0

* Initial release version.
