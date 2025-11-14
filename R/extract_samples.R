#' Extract samples (draws) from a fitted \code{coevfit} object
#'
#' @srrstats {G1.4} Function is documented
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
  # get variables and draws
  vars <- object$fit$metadata()$stan_variables
  draws_array <- object$fit$draws()
  
  # Convert to rvars format (groups indexed variables into arrays)
  # e.g., A_diag[1], A_diag[2] -> A_diag array
  draws_rvars <- posterior::as_draws_rvars(draws_array)
  
  # Extract each base variable (not individual indexed variables)
  # rvars groups indexed variables, so we extract by base name
  base_vars <- names(draws_rvars)
  
  samples_list <- lapply(base_vars, function(var_name) {
    # Get the rvar object
    var_rvar <- draws_rvars[[var_name]]
    
    # Extract draws as array (collapses chains)
    # draws_of() returns arrays where first dimension is draws
    var_draws <- posterior::draws_of(var_rvar, with_chains = FALSE)
    
    # Return as vector or array depending on dimensions
    # For scalar variables: return as vector
    # For array variables: return as array/matrix
    if (is.array(var_draws)) {
      n_dims <- length(dim(var_draws))
      if (n_dims == 1) {
        # 1D array -> vector
        return(as.vector(var_draws))
      } else {
        # Multi-dimensional array -> keep as array
        # First dimension is draws, remaining dimensions are parameter dimensions
        return(var_draws)
      }
    } else {
      # Not an array -> return as vector
      return(as.vector(var_draws))
    }
  })
  
  # Set names
  names(samples_list) <- base_vars
  
  return(samples_list)
}
