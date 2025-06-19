#' Extract samples (draws) from a fitted \code{coevfit} object
#'
#' @param object An object of class \code{coevfit}
#'
#' @returns Samples in 'rethinking' style list format
#'
#' @export
extract_samples <- function(object) {
  UseMethod("extract_samples")
}

#' @export
extract_samples.coevfit <- function(object) {
  # get variables and rvars draws
  vars <- object$fit$metadata()$stan_variables
  draws <- posterior::as_draws_rvars(object$fit$draws())
  # reshape to 'rethinking' style list format
  lapply(
    vars,
    \(var_name){
      posterior::draws_of(
        draws[[var_name]],
        with_chains = FALSE
        )
      }
    ) |>
  stats::setNames(vars)
}
