#' Predict a co-evolutionary time series from a fitted \code{coevfit} object
#'
#' This function produces predicted values of traits over a time series, given
#' the estimated parameters for a \code{coevfit} model. This function is used
#' under the hood by the plotting function \code{\link{coev_pred_series}}.
#'
#' @param object An object of class \code{coevfit}
#' @param eta_anc If \code{NULL} (default), the starting values for the latent
#'   states \eqn{\eta} will be set to the estimated ancestral states at the
#'   root from the \code{coevfit} model. Otherwise, a named numeric vector with
#'   length equal to the number of variables specifying the initial \eqn{\eta}
#'   values. All variable names must be included in this vector.
#' @param tmax A positive number indicating the total duration of time for which
#'   to predict. Set to 1 by default, corresponding to the entire time depth of
#'   the original phylogenetic tree(s).
#' @param ntimes A positive integer indicating the total number of discrete time
#'   steps over which to make predictions. Each step corresponds to a time
#'   difference of dt = tmax/ntimes. Set to 30 by default.
#' @param ndraws An integer indicating the number of draws to return. The
#'   default and maximum number of draws is the size of the posterior sample.
#' @param stochastic Logical (defaults to \code{FALSE}); indicator of whether
#'   predictions should include only the expected co-evolutionary change
#'   due to deterministic selection (FALSE) or also stochastic drift (TRUE).
#'
#' @return An \[ndraw, ntimes, nvariables\] array of predicted \eqn{\eta}
#'   values.
#'
#' @author Scott Claessens \email{scott.claessens@@gmail.com}, Erik Ringen
#'   \email{erikjacob.ringen@@uzh.ch}
#'
#' @seealso \code{\link{coev_plot_pred_series}}
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
#'   stochastic = TRUE
#'   )
#'
#' # expected trait co-evolution, no drift, specify ancestral states (initial values)
#' epreds <- coev_pred_series(
#'   object = fit,
#'   stochastic = FALSE,
#'   eta_anc = c("political_authority" = -2, "religious_authority" = 1.5)
#'   )
#' }
#'
#' @export
coev_pred_series <- function(object, eta_anc = NULL, tmax = 1, ntimes = 30,
                             ndraws = NULL, stochastic = FALSE) {
  # stop if object is not of class coevfit
  if (!methods::is(object, "coevfit")) {
    stop2(
      paste0(
        "Argument 'object' must be a fitted coevolutionary model of class ",
        "coevfit."
        )
      )
  }
  if (!is.null(eta_anc)) {
    # stop if eta_anc is neither NULL nor
    # a named vector of length equal to the number of variables
    if (!is.numeric(eta_anc) | length(eta_anc) != length(object$variables)) {
        stop2(
          paste0(
            "Argument 'eta_anc' must be numeric and equal in length to the ",
            "number of variables."
          )
        )
    }
    # stop if eta_anc is not a named vector
    if (is.null(names(eta_anc))) {
      stop2("Argument 'eta_anc' is not a named vector.")
    }
    # stop if eta_anc does not have names equal to variables in the model
    if (!identical(sort(names(eta_anc)), sort(names(object$variables)))) {
      stop2(
        paste0(
          "Argument 'eta_anc' has names different to the variables included ",
          "in the model."
        )
      )
    }
  }
  # stop if tmax is not numeric or not positive
  if (!is.numeric(tmax) | length(tmax) != 1) {
    stop2("Argument 'tmax' must be a single numeric value.")
  } else if (tmax <= 0) {
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
  # get posterior samples and number of variables
  post <- extract_samples(object)
  J <- length(object$variables)
  # get eta_anc in the correct order
  eta_anc <- eta_anc[names(object$variables)]
  # get number of samples
  nsamps <- ifelse(is.null(ndraws), length(post$lp__), ndraws)
  # initialise empty array for predictions
  preds <-
    array(
      NA,
      dim = c(nsamps, ntimes + 1, J),
      dimnames = list(
        samps = 1:nsamps,
        time = 1:(ntimes + 1),
        response = names(object$variables)
        )
      )
  # if user defined ancestral states
  # otherwise, use model inferred ancestral states
  if (!is.null(eta_anc)) {
    for (i in 1:nsamps) preds[i, 1, ] = eta_anc
  } else {
    # set posterior ancestral states
    eta_anc_long <- post$eta_anc
    ntrees <- dim(eta_anc_long)[2]
    # shuffle dimensions of tree dimension to get random draws
    eta_anc_long <-
      eta_anc_long[, sample(1:ntrees, size = ntrees, replace = FALSE), ]
    if (ntrees > 1) {
      # make into 2d matrix by stacking trees
      eta_anc_long2 <- eta_anc_long[, 1, ]
      for (t in 2:ntrees) {
        eta_anc_long2 <- rbind(eta_anc_long2, eta_anc_long[,t,])
        eta_anc_long <- eta_anc_long2
      }
    }
    for (i in 1:nsamps) {
      preds[i, 1, ] = eta_anc_long[i, ]
    }
  }
  # calculate predictions
  for (i in 1:nsamps) {
    # selection parameters
    A <- post$A[i,,]
    A_delta <- as.matrix(Matrix::expm(A * tmax / ntimes))
    inv_A <- solve(A)
    I <- diag(rep(1, J))
    b <- post$b[i,]
    # drift parameters
    Q_inf <- post$Q_inf[i,,]
    VCV <- Q_inf - ((A_delta) %*% Q_inf %*% t(A_delta))
    chol_VCV <- t(chol(Matrix::nearPD(VCV)$mat))
    # calculate predictions over time
    for (t in 1:ntimes) {
      preds[i, t + 1, ] <- (A_delta %*% preds[i, t, ] +
                              (inv_A %*% (A_delta - I) %*% b))[,1]
      # add drift
      if (stochastic == TRUE) {
        preds[i, t + 1, ] <-
          preds[i, t + 1, ] + (chol_VCV %*% stats::rnorm(J, 0, 1))
      }
    }
  }
  # shuffle draws to de-correlate when doing simulations
  preds[sample(1:nsamps, size = nsamps, replace = FALSE), , ]
}
