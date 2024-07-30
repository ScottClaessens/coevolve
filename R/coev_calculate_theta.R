#' Calculate optimal trait values (theta) for a fitted \code{coevfit} object
#'
#' @param object An object of class \code{coevfit}
#' @param intervention_values A named list of variables and associated
#'   intervention values for calculating optimal trait values. All coevolving
#'   variables must be declared separately in the named list without repetition.
#'   If the intervention value for a particular variable is set to NA, this
#'   variable is treated as a free variable. Otherwise, if the intervention
#'   value for a particular variable is specified, the variable is held
#'   constant at this trait value in the calculation. At least one variable must
#'   be declared as a free variable and at least one variable must be held
#'   constant (e.g., \code{list(var1 = NA, var2 = 0)}).
#'
#' @return Posterior samples in matrix format
#' @export
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
#' # calculate theta given intervention values
#' coev_calculate_theta(
#'   object = fit,
#'   intervention_values = list(
#'     political_authority = NA,
#'     religious_authority = 0
#'     )
#'   )
#' }
coev_calculate_theta <- function(object, intervention_values) {
  # stop if object is not of class coevfit
  if (!methods::is(object, "coevfit")) {
    stop2(
      paste0(
        "Argument 'object' must be a fitted coevolutionary model ",
        "of class coevfit."
        )
      )
  }
  # stop if intervention_values argument is not a named list
  if (!is.list(intervention_values) | is.null(names(intervention_values))) {
    stop2("Argument 'intervention_values' is not a named list.")
  }
  # stop if intervention_list contains variables not included in the model
  if (!all(names(intervention_values) %in% names(object$variables))) {
    stop2(
      paste0(
        "At least one variable in 'intervention_values' is not included in ",
        "the fitted model."
      )
    )
  }
  # stop if any coevolving variables are not listed in intervention_values
  if (any(!(names(object$variables) %in% names(intervention_values)))) {
    stop2(
      paste0(
        "All coevolving variables must be included in ",
        "argument 'intervention_values'."
        )
      )
  }
  # stop if repetition in intervention_values
  if (any(duplicated(names(intervention_values)))) {
    stop2("Argument 'intervention_values' contains duplicated variable names.")
  }
  # stop if any values in intervention_list are not of length one
  if (any(unlist(lapply(intervention_values, length)) != 1)) {
    stop2("Values in 'intervention_values' must each be of length one.")
  }
  # stop if any values in intervention_list are not NA or numeric
  if (any(unlist(lapply(intervention_values,
                        function(x) !is.na(x) & !is.numeric(x))))) {
    stop2("Values in 'intervention_values' must each be NA or numeric.")
  }
  # stop if all variables are held constant in intervention_list
  if (mean(is.na(intervention_values)) == 0) {
    stop2(
      paste0(
        "Argument 'intervention_values' must have at least one NA value ",
        "declaring a free variable. If all variables are held constant, the ",
        "system is already at equilibrium and there is nothing to compute."
        )
      )
  }
  # stop if all variables in intervention_values are free
  if (mean(is.na(intervention_values)) == 1) {
    stop2(
      paste0(
        "Argument 'intervention_values' must have at least one variable ",
        "held constant (i.e., all values are NA)."
        )
      )
  }
  # construct intervention values x_hat
  x_hat <- unlist(intervention_values[names(object$variables)])
  # extract posterior draws
  post <- posterior::as_draws_rvars(object$fit$draws(variables = c("A", "b")))
  A <- posterior::draws_of(post$A)
  b <- posterior::draws_of(post$b)
  # partition the matrices
  held_indices <- which(!is.na(x_hat))
  free_indices <- which(is.na(x_hat))
  A_free_free  <- A[,free_indices, free_indices, drop = FALSE]
  A_free_held  <- A[,free_indices, held_indices, drop = FALSE]
  b_free       <- b[,free_indices, drop = FALSE]
  # construct theta matrix
  theta <- matrix(NA, nrow = posterior::ndraws(post), ncol = length(x_hat))
  for (i in 1:posterior::ndraws(post)) {
    # compute the equilibrium for the free variables
    if (length(held_indices) == 1) {
      equilibrium_free <-
        -solve( A_free_free[i,,] ) %*%
        (b_free[i,] + as.matrix(A_free_held[i,,]) %*% x_hat[!is.na(x_hat)])
    } else {
      equilibrium_free <-
        -solve( A_free_free[i,,] ) %*%
        (b_free[i,] + A_free_held[i,,] %*% x_hat[!is.na(x_hat)])
    }
    # initialize the equilibrium vector
    equilibrium <- rep(NA, length(x_hat))
    # fill in the computed equilibrium values for free variables
    equilibrium[free_indices] <- equilibrium_free
    # fill in the constant values for held variables
    equilibrium[held_indices] <- x_hat[!is.na(x_hat)]
    # add to theta matrix
    theta[i,] = equilibrium
  }
  # add column names to theta matrix
  colnames(theta) <- names(object$variables)
  return(theta)
}
