#' Extract samples (draws) from a fitted \code{coevfit} object
#'
#' @param object An object of class \code{coevfit}
#'
#' @returns Samples in 'rethinking' style list format
#'
#' @examples
#' \dontrun{
#' # fit dynamic coevolutionary model
#' fit <- coev_fit(
#'   data = authority$data,
#'   variables = list(
#'     political_authority = "ordered_logistic",
#'     religious_authority = "ordered_logistic"
#'   ),
#'   id = "language",
#'   tree = authority$phylogeny,
#'   # additional arguments for cmdstanr::sample()
#'   chains = 4,
#'   parallel_chains = 4,
#'   seed = 1
#'   )
#'
#' # extract samples from fitted model
#' extract_samples(fit)
#' }
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
    \(var_name) {
      posterior::draws_of(
        draws[[var_name]],
        with_chains = FALSE
      )
    }
  ) |>
    stats::setNames(vars)
}
