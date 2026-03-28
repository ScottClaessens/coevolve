# coevolve (development version)

### New Features

* Added pure JAX/NumPyro backend via `nuts_sampler = "nutpie"`.
  Uses nutpie's Rust NUTS sampler with JAX gradients for ~5x
  faster sampling than Stan on typical models. Requires
  `pip install jax numpyro nutpie`.

### Bug Fixes

* Fixed issue with `summary()` when `estimate_residual = FALSE` (#95)

### New Features

* Re-use matrix computations when branch lengths are identical (#93)

### Other Changes

* Updated license
* Submitted to [rOpenSci](https://ropensci.org/) for package review
* Linted package using [lintr](https://lintr.r-lib.org/)
* Implemented [srr](https://docs.ropensci.org/srr/) compliance checks
* Added continuous integration checks using GitHub Actions
* Reduced cyclomatic complexity for some functions

# coevolve 0.1.0

* Initial release version.
