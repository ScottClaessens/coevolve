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
#'   be declared as a free variable (e.g., \code{list(var1 = NA, var2 = 0)}).
#'
#' @return Posterior samples in matrix format
#' @export
#'
#' @examples
#' \dontrun{
#' # simulate data
#' n <- 20
#' tree <- ape::rcoal(n)
#' d <- data.frame(
#'   id = tree$tip.label,
#'   x = rbinom(n, size = 1, prob = 0.5),
#'   y = ordered(sample(1:4, size = n, replace = TRUE))
#' )
#' # fit dynamic coevolutionary model
#' m <- coev_fit(
#'   data = d,
#'   variables = list(
#'     x = "bernoulli_logit",
#'     y = "ordered_logistic"
#'   ),
#'   id = "id",
#'   tree = tree,
#'   # additional arguments for cmdstanr::sample()
#'   chains = 4,
#'   parallel_chains = 4,
#'   iter_warmup = 500,
#'   seed = 1
#' )
#' # calculate theta given intervention values
#' coev_calculate_theta(m, list(x = NA, y = 0))
#' }
coev_calculate_theta <- function(object, x_hat) {
  # stop if object is not of class coevfit
  if (!methods::is(object, "coevfit")) {
    stop2("Argument 'object' must be a fitted coevolutionary model of class coevfit.")
  }

  # stop if all variables are held constant
  if (mean(is.na(x_hat)) == 0) {
      stop2("Argument 'x_hat' must have at least one NA value, signifying a free variable. If all variables held constant, the system is already at equilibrium--nothing to compute!.")
  }

   # warn all variables are free
   if (mean(is.na(x_hat)) == 1) {
    warnings("Argument 'x_hat' has no variables held constant (i.e., all values are NA), is this what you intended?")
   }

  # extract posterior draws
  post <- posterior::as_draws_rvars(object$fit$draws(variables = c("A", "b")))

  A <- posterior::draws_of(post$A)
  b <- posterior::draws_of(post$b)

  # stop if object if not enough x_hat values specified
  if (dim(A)[2] != length(x_hat)) {
    stop2("Argument 'x_hat' must be the same length as the number of variables in the fit model.")
  }

# Partition the matrices
held_indices <- which(!is.na(x_hat))
free_indices <- which(is.na(x_hat))
A_ff <- A[,free_indices, free_indices, drop = F]
b_f <- b[,free_indices, drop = F]
A_fc <- A[,free_indices, held_indices, drop = F]

theta <- matrix(NA, nrow = ndraws(post), ncol = length(x_hat))

for (i in 1:ndraws(post)) {
  # Compute the equilibrium for the free variables
  if (length(held_indices) == 1) equilibrium_f <- -solve( A_ff[i,,] ) %*% (b_f[i,] + as.matrix(A_fc[i,,]) %*% x_hat[!is.na(x_hat)])
  else equilibrium_f <- -solve( A_ff[i,,] ) %*% (b_f[i,] + A_fc[i,,] %*% x_hat[!is.na(x_hat)])

  # Initialize the equilibrium vector
  equilibrium <- rep(NA, length(x_hat))

  # Fill in the computed equilibrium values for free variables
  equilibrium[free_indices] <- equilibrium_f

  # Fill in the constant values for held variables
  equilibrium[held_indices] <- x_hat[!is.na(x_hat)]

  theta[i,] = equilibrium
}

return( theta )
}


