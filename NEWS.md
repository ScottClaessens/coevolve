# coevolve (development version)

### New Features

* Added `coev_ancestral_states()` for extracting posterior estimates of
  ancestral trait values at internal phylogenetic nodes, on either the
  latent or response scale (#86)
* Added `coev_plot_node_labels()` for plotting the tree with internal
  node IDs displayed
* Added `coev_identify_nodes()` for interactive multi-click node selection,
  wrapping `ape::identify.phylo()` in a loop
* Added vignette "Ancestral State Reconstruction with coevolve"
* Allow for single traits (#107)

# coevolve 1.0.0

### Bug Fixes

* Fixed issue with `summary()` when `estimate_residual = FALSE` (#95)

### New Features

* Added `nutpie` as an alternative sampler for the Stan models (#92)
* Cached matrix computations when branch lengths are identical (#93)
* Implemented Hilbert-space approximate Gaussian processes for spatial control,
  adding the `lon_lat` argument and deprecating the `dist_mat` argument in 
  `coev_fit()` (#103)

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
