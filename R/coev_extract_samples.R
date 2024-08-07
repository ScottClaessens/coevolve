#' Extract samples (draws) from a fitted \code{coevfit} object
#'
#' @param object An object of class \code{coevfit}
#'
#' @return Samples in 'rethinking' style list format
#' @export
#' 
coev_extract_samples <- function(object) {
  # stop if object is not of class coevfit
  if (!methods::is(object, "coevfit")) {
    stop2(
      paste0(
        "Argument 'object' must be a fitted coevolutionary model ",
        "of class coevfit."
        )
      )
  }
  # stop if rethinking package not installed
  if (!requireNamespace("rethinking", quietly = TRUE)) {
    stop2(
      paste0(
        "'rethinking' package must be installed to use this function. See https://github.com/rmcelreath/rethinking"
        )
      )
  }
    vars <- object$fit$metadata()$stan_variables
    draws <- posterior::as_draws_rvars(object$fit$draws())
    
    return(lapply(vars, \(var_name){  
      posterior::draws_of(draws[[var_name]], with_chains = FALSE)
    }) |> setNames(vars))
}
