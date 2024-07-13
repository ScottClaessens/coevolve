#' Calculate optimal trait values (theta) given a fitted \code{coevfit} object and a vector of "intervention" trait values x-hat, which are assumed to be held constant.
#'
#' @param object An object of class \code{coevfit}
#' @param x_hat A numeric vector of intervention values where the elements correspond to the varibables of the coevfit. If an element is NA, then
#'
#' @return Posterior samples in matrix format
#' @export
#'

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


