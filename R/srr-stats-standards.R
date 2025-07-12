#' srr_stats
#'
#' @srrstatsVerbose TRUE
#'
#' @srrstatsTODO {BS7.0} *Software should demonstrate and confirm recovery of parametric estimates of a prior distribution*
#' @srrstatsTODO {BS7.1} *Software should demonstrate and confirm recovery of a prior distribution in the absence of any additional data or information*
#' @srrstatsTODO {BS7.2} *Software should demonstrate and confirm recovery of a expected posterior distribution given a specified prior and some input data*
#' @srrstatsTODO {BS7.3} *Bayesian software should include tests which demonstrate and confirm the scaling of algorithmic efficiency with sizes of input data.*
#' @srrstatsTODO {BS7.4} *Bayesian software should implement tests which confirm that predicted or fitted values are on (approximately) the same scale as input values.*
#' @srrstatsTODO {BS7.4a} *The implications of any assumptions on scales on input objects should be explicitly tested in this context; for example that the scales of inputs which do not have means of zero will not be able to be recovered.*
#' @noRd
NULL

#' NA_standards
#'
#' @srrstatsNA {G1.5} We have not yet made performance claims in associated
#'   publications
#' @srrstatsNA {G1.6} There are no alternative implementations of this algorithm
#'   in other R packages
#' @srrstatsNA {G2.4d} We do not explicitly convert input data to factors in the
#'   code
#' @srrstatsNA {G2.9} The package does not do any type conversions that would
#'   result in lost data
#' @srrstatsNA {G2.14a} NAs are important information, and so the package does
#'   not error when they are present
#' @srrstatsNA {G3.0} The package does not compare floating points for equality
#' @srrstatsNA {G5.0} This package is designed to be used with specific data
#'   (i.e., phylogenetic data) so NIST datasets would not be applicable
#' @srrstatsNA {G5.4b} This package is not a new implementation of an existing
#'   method
#' @srrstatsNA {G5.4c} It is not necessary to draw stored values from published
#'   outputs
#' @srrstatsNA {G5.11, G5.11a} Extended tests do not require large data sets or
#'   additional asset downloads
#' @srrstatsNA {BS1.0} The term "hyperparameter" is not used
#' @srrstatsNA {BS1.3a, BS2.8} It is not possible to use the output of previous
#'   simulations as starting points of subsequent simulations with this package
#' @srrstatsNA {BS1.5, BS4.6, BS4.7, BS5.4} This package enables convergence
#'   checks through cmdstanr, but does not enable multiple types of convergence
#'   checks itself
#' @srrstatsNA {BS2.10, BS2.11} Seeds and starting values are handled by
#'   additional arguments to coev_fit() which are passed to
#'   cmdstanr::sample().
#' @srrstatsNA {BS3.1, BS3.2} This package does not check for (or error with)
#'   perfectly collinear input data, as these data could still be valid outcomes
#'   of a coevolutionary process (e.g., strong reciprocal coevolution)
#' @srrstatsNA {BS4.1} We do not compare to alternative samplers
#' @srrstatsNA {BS4.4} cmdstanr::sample() does not have an option to stop chains
#'   upon convergence
#'
#' @noRd
NULL
