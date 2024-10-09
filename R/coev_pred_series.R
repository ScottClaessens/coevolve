#' Predict a co-evolutionary time series from a fitted \code{coevfit}
#' object
#'
#' @param object An object of class \code{coevfit}
#' @param eta_anc If \code{NULL} default, the starting values for the latent states \eqn{\eta} will be set to the estimated ancestral states at the root from the the \code{coevfit}. Otherwise, a numeric vector with length equal to the number of variable specifying the initial \eqn{\eta}  values.
#' @param tmax A positive number indicating the total duration of time for which to predict. Set to 1 by default, corresponding to the entire time depth of the original phylogenetic tree(s).
#' @param tmax A positive number indicating the total duration of time for which to predict. Set to 1 by default, corresponding to the entire time depth of the original phylogenetic tree(s).
#' @param ntimes A positive integer indicating the total number of discrete time steps over which to make predictions. Each step corresponds to a time difference of dt = tmax/ntimes. Set to 30 by default.
#' @param ndraws An integer indicating the number of draws to return. The
#'   default and maximum number of draws is the size of the posterior sample.
#' @param stochastic Logical (defaults to \code{FALSE}); indicator of whether predictions should include only the expected co-evolutionary change (deterministic selection), or also stochastic drift.
#' @return An [ndraw, ntimes, nvariables] array of predicted \eqn{\eta}  values.
#'
#' @author Scott Claessens \email{scott.claessens@@gmail.com}, Erik Ringen
#'   \email{erikjacob.ringen@@uzh.ch}
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
#' # simulate trait co-evolution
#' sims <- coev_pred_series(
#'   object = fit,
#'   stochastic = T
#'   )
#' 
#' # expected trait co-evolution, no drift
#' epreds <- coev_pred_series(
#'   object = fit,
#'   stochastic = F
#'   )
#' }
#'
#' @export
coev_pred_series <- function(object, eta_anc = NULL, tmax = 1, ntimes = 30, ndraws = NULL, stochastic = FALSE){
  if (!methods::is(object, "coevfit")) {
    stop2(
      paste0(
        "Argument 'object' must be a fitted coevolutionary model of class coevfit."
        )
      )
  }
  # stop if eta_anc is neither NULL nor a vector of length equal to the number of variables
  if (!is.null(eta_anc) & (!is.numeric(eta_anc) | length(eta_anc) != length(object$variables))) {
    stop2("Argument 'eta_anc' must be numeric and equal in length to the number of variables.")
  }
  # stop if tmax is not numeric or not positive
  if (!is.numeric(tmax)|length(tmax) != 1) {
    stop2("Argument 'tmax' must be a single numeric value.")
  }
  else if (tmax <= 0) {
    stop2("Argument 'tmax' must be positive.")
  }
  # stop if ndraws is not a single integer between 1 and the total num draws
  if (!is.null(ndraws)) {
    if (!is.integer(ndraws) | length(ndraws) != 1) {
      stop2("Argument 'ndraws' must be a single integer.")
    } else if (ndraws < 1 | ndraws > nrow(object$fit$draws())) {
      stop2(
        "Argument 'ndraws' must be between 1 and the total number of draws."
        )
    }
  }
# stop if stochastic not logical
  if (!is.logical(stochastic)) {
    stop2("Argument 'stochastic' must be logical.")
  }

  post <- extract_samples(object)
  J <- length(object$variables)

  nsamps <- ifelse(is.null(ndraws), length(post$lp__), ndraws)

  preds <- array(NA, dim = c(nsamps, ntimes+1, J), dimnames = list(samps = 1:nsamps, time = 1:(ntimes+1), response = names(object$variables)))

  # User defined ancestral states
  if (!is.null(eta_anc)) {
    for (i in 1:nsamps) preds[i,1,] = eta_anc
  }

  # Model inferred ancestral states
  else {
    eta_anc_long <- post$eta_anc
    ntrees <- dim(eta_anc_long)[2]
    # shuffle dimenisions of tree dimension to get random draws
    eta_anc_long <- eta_anc_long[,sample(1:ntrees, size = ntrees, replace = F),]

    if (ntrees > 1) {
    # make into 2d matrix by stacking trees
    eta_anc_long2 <- eta_anc_long[,1,]
    for (t in 2:ntrees) eta_anc_long2 <- rbind(eta_anc_long2, eta_anc_long[,t,]) -> eta_anc_long
    }
    for (i in 1:nsamps) preds[i,1,] = eta_anc_long[i,]
    }

  for (i in 1:nsamps) {
    A <- post$A[i,,]
    A_delta <- as.matrix(Matrix::expm(A * tmax/ntimes))
    inv_A <- solve(A)
    I <- diag(rep(1, J))
    b <- post$b[i,]

    Q_inf <- post$Q_inf[i,,]
    VCV <- Q_inf - ((A_delta) %*% Q_inf %*% t(A_delta))
    chol_VCV <- t(chol(VCV))

    for (t in 1:ntimes) {
      preds[i, t+1, ] = (A_delta %*% preds[i, t, ] + (inv_A %*% (A_delta - I) %*% b))[,1]
      # add drift
      if (stochastic == T) preds[i, t+1, ] = preds[i, t+1, ] + chol_VCV %*% rnorm(J, 0, 1)
    }
  }
  # shuffle draws to de-correlate when doing simulations
  return(preds[sample(1:nsamps, size = nsamps, replace = F),,])
}
