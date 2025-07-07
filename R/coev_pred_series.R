#' Predict a co-evolutionary time series from a fitted \code{coevfit} object
#'
#' This function produces predicted values of traits over a time series, given
#' the estimated parameters for a \code{coevfit} model. This function is used
#' under the hood by the plotting function \code{\link{coev_pred_series}}.
#'
#' @srrstats {G1.3, G1.4, G2.1a} Function documentation begins here, with
#'   expected data types and definitions of statistical terminology and inputs
#'
#' @param object An object of class \code{coevfit}
#' @param eta_anc If \code{NULL} (default), the starting values for the latent
#'   states \eqn{\eta} will be set to the estimated ancestral states at the
#'   root from the \code{coevfit} model. Otherwise, a named list with
#'   length equal to the number of variables specifying the initial \eqn{\eta}
#'   values. All variable names must be included in this list.
#' @param intervention_values Either \code{NULL} (the default) or a named list
#'   of variables and associated intervention values. If \code{NULL}, all traits
#'   are free to vary. Otherwise, all coevolving variables must be
#'   declared separately in the named list without repetition. If the
#'   intervention value for a particular variable is set to NA, this variable is
#'   treated as a free variable. If the intervention value for a particular
#'   variable is specified, the variable is held constant at this trait value in
#'   the calculation. At least one variable must be declared as a free variable
#'   and at least one variable must be held constant (e.g.,
#'   \code{list(var1 = NA, var2 = 0)}).
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
#'   Incompatible with the 'intervention_values' argument.
#'
#' @returns An \[ndraw, ntimes, nvariables\] array of predicted \eqn{\eta}
#'   values
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
#' # expected trait co-evolution, no drift
#' epreds <- coev_pred_series(
#'   object = fit,
#'   stochastic = FALSE,
#'   eta_anc = list(political_authority = -2, religious_authority = 1.5)
#'   )
#'
#' # expected trait co-evolution under intervention, no drift
#' epreds_intervention <- coev_pred_series(
#'   object = fit,
#'   stochastic = FALSE,
#'   intervention_values = list(political_authority = NA,
#'                              religious_authority = 1.5)
#'   )
#' }
#'
#' @export
coev_pred_series <- function(object, eta_anc = NULL, intervention_values = NULL,
                             tmax = 1, ntimes = 30, ndraws = NULL,
                             stochastic = FALSE) {
  # run checks
  run_checks_pred_series(object, eta_anc, intervention_values,
                         tmax, ntimes, ndraws, stochastic)
  # get posterior samples and number of variables
  post <- extract_samples(object)
  j <- length(object$variables)
  # get number of samples
  nsamps <- ifelse(is.null(ndraws), length(post$lp__), ndraws)
  # initialize empty array for predictions
  preds <-
    array(
      NA,
      dim = c(nsamps, ntimes + 1, j),
      dimnames = list(
        samps = 1:nsamps,
        time = 1:(ntimes + 1),
        response = names(object$variables)
      )
    )
  # handle intervention_values
  if (!is.null(intervention_values)) {
    # construct intervention values vector x_hat
    x_hat <- unlist(intervention_values[names(object$variables)])
    # partition indices into held and free variables
    held_indices <- which(!is.na(x_hat))
    free_indices <- which(is.na(x_hat))
    # check if there is at least one free variable
    if (length(free_indices) == 0) {
      stop2(
        paste0(
          "At least one variable must be declared as a free variable ",
          "(with NA in 'intervention_values')."
        )
      )
    }
  } else {
    # all variables are free
    x_hat <- rep(NA, j)
    held_indices <- integer(0)
    free_indices <- 1:j
  }
  if (is.null(intervention_values)) {
    initial_values <- rep(NA, j)
  } else {
    initial_values <- unlist(intervention_values[names(object$variables)])
  }
  # handle eta_anc
  if (!is.null(eta_anc)) {
    # user provided initial values
    eta_anc <- unlist(eta_anc[names(object$variables)])
    # check for conflicts between eta_anc and intervention_values
    conflicting_vars <- intersect(
      names(eta_anc)[!is.na(eta_anc)],
      names(x_hat)[!is.na(x_hat)]
    )
    if (length(conflicting_vars) > 0) {
      message(
        paste0(
          "Note: For variable(s) ", paste(conflicting_vars, collapse = ", "),
          ", both 'eta_anc' and 'intervention_values' specify non-NA values. ",
          "The 'intervention_values' will take precedence for these ",
          "variable(s)."
        )
      )
      # override eta_anc values for conflicting variables
      # with intervention_values
      eta_anc[conflicting_vars] <- x_hat[conflicting_vars]
    }
    # set initial values for variables
    initial_values <- eta_anc
  }
  # get model inferred ancestral states, if necessary
  if (any(is.na(initial_values))) {
    # Use estimated ancestral states from the model
    eta_anc_long <- post$eta_anc
    ntrees <- dim(eta_anc_long)[2]
    # shuffle dimensions of tree dimension to get random draws
    eta_anc_long <-
      eta_anc_long[, sample(1:ntrees, size = ntrees, replace = FALSE), ]
    if (ntrees > 1) {
      # make into 2D matrix by stacking trees
      eta_anc_long2 <- eta_anc_long[, 1, ]
      for (t in 2:ntrees) {
        eta_anc_long2 <- rbind(eta_anc_long2, eta_anc_long[, t, ])
      }
      eta_anc_long <- eta_anc_long2
    }
  }
  for (j_index in 1:j) {
    for (i in 1:nsamps) {
      if (is.na(initial_values[j_index])) {
        preds[i, 1, j_index] <- eta_anc_long[i, j_index]
      } else {
        preds[i, 1, j_index] <- initial_values[j_index]
      }
    }
  }
  # iterate over each sample
  for (i in 1:nsamps) {
    # selection parameters for this sample
    a <- post$A[i, , ]
    b <- post$b[i, ]
    # partition A and b into free and held
    a_free_free <- a[free_indices, free_indices, drop = FALSE]
    if (length(held_indices) > 0) {
      a_free_held <- a[free_indices, held_indices, drop = FALSE]
    } else {
      a_free_held <- matrix(0, nrow = length(free_indices), ncol = 0)
    }
    b_free <- b[free_indices]
    # compute c = A_free_held * x_hat_held + b_free
    if (length(held_indices) > 0) {
      c <- a_free_held %*% x_hat[held_indices] + b_free
    } else {
      c <- b_free
    }
    # compute A_delta_free_free and inv_A_free_free
    a_delta_free_free <- as.matrix(Matrix::expm(a_free_free * tmax / ntimes))
    # ensure A_free_free is square
    if (nrow(a_free_free) != ncol(a_free_free)) {
      stop2("Matrix A_free_free must be square.")
    }
    # Invert A_free_free
    inv_a_free_free <- tryCatch(
      solve(a_free_free),
      error = function(e) {
        stop2("Matrix A_free_free is singular and cannot be inverted.")
      }
    )
    # identity matrix
    i_free_free <- diag(rep(1, length(free_indices)))
    # drift parameters, cannot be used in conjunction
    # with intervention values (currently)
    if (stochastic == TRUE) {
      q_inf <- post$Q_inf[i, , ]
      vcv <- q_inf - ((a_delta_free_free) %*% q_inf %*% t(a_delta_free_free))
      chol_vcv <- t(chol(Matrix::nearPD(vcv)$mat))
    }
    # initialize preds_free with current state
    preds_free <- preds[i, 1, free_indices]
    # iterate over each time step
    for (t in 1:ntimes) {
      # compute expected change for free variables
      preds_free <-
        (a_delta_free_free %*% preds_free +
         (inv_a_free_free %*% (a_delta_free_free - i_free_free) %*% c))[, 1]
      # add drift if stochastic
      if (stochastic == TRUE) {
        preds_free <-
          preds_free + (chol_vcv %*% stats::rnorm(length(free_indices), 0, 1))
      }
      # update preds for the next time point
      preds[i, t + 1, free_indices] <- preds_free
      # set held variables to intervention values
      if (length(held_indices) > 0) {
        preds[i, t + 1, held_indices] <- x_hat[held_indices]
      }
    }
  }
  return(preds)
}

#' Internal helper function for checking coev_pred_series() arguments
#'
#' @srrstats {G1.4a} Non-exported function documented here
#'
#' @description Checks arguments for coev_pred_series()
#'
#' @returns Error message if any of the checks fail
#'
#' @noRd
run_checks_pred_series <- function(object, eta_anc, intervention_values,
                                   tmax, ntimes, ndraws, stochastic) {
  #' @srrstats {G5.2, G5.2a} Unique error messages for each input
  # stop if object is not of class coevfit
  #' @srrstats {G2.1} Assertion on type of input
  if (!methods::is(object, "coevfit")) {
    stop2(
      paste0(
        "Argument 'object' must be a fitted coevolutionary model of class ",
        "coevfit."
      )
    )
  }
  if (!is.null(eta_anc)) {
    # stop if eta_anc argument is not a named list
    #' @srrstats {G2.1} Assertion on type of input
    if (!is.list(eta_anc) || is.null(names(eta_anc))) {
      stop2("Argument 'eta_anc' is not a named list.")
    }
    # stop if eta_anc contains variables not included in the model
    if (!all(names(eta_anc) %in% names(object$variables))) {
      stop2(
        paste0(
          "At least one variable in 'eta_anc' is not included in ",
          "the fitted model."
        )
      )
    }
    # stop if any coevolving variables are not listed in eta_anc
    if (any(!(names(object$variables) %in% names(eta_anc)))) {
      stop2(
        paste0(
          "All coevolving variables must be included in ",
          "argument 'eta_anc'."
        )
      )
    }
    # stop if repetition in eta_anc
    if (any(duplicated(names(eta_anc)))) {
      stop2(
        "Argument 'eta_anc' contains duplicated variable names."
      )
    }
    # stop if any values in eta_anc are not of length one
    if (any(unlist(lapply(eta_anc, length)) != 1)) {
      stop2("Values in 'eta_anc' must each be of length one.")
    }
    # stop if any values in eta_anc are not numeric
    if (any(unlist(lapply(eta_anc,
                          function(x) !is.numeric(x))))) {
      stop2("Values in 'eta_anc' must each be numeric.")
    }
  }
  # stop if tmax is not numeric or not positive
  #' @srrstats {G2.1} Assertion on type of input
  if (!is.numeric(tmax) || length(tmax) != 1) {
    stop2("Argument 'tmax' must be a single numeric value.")
  } else if (tmax <= 0) {
    stop2("Argument 'tmax' must be positive.")
  }
  # stop if ndraws is not a single integer between 1 and the total num draws
  if (!is.null(ndraws)) {
    if (!is.numeric(ndraws)) {
      #' @srrstats {G2.1} Assertion on type of input
      stop2("Argument 'ndraws' must be numeric.")
    } else if (!all(as.integer(ndraws) == ndraws) || length(ndraws) != 1) {
      #' @srrstats {G2.0, G2.1, G2.2, G2.4, G2.4a} Assertion on length and type
      #' of input, convert to integer
      stop2("Argument 'ndraws' must be a single integer.")
    } else if (ndraws < 1 || ndraws > nrow(object$fit$draws())) {
      stop2(
        "Argument 'ndraws' must be between 1 and the total number of draws."
      )
    }
  }
  # stop if stochastic not logical
  #' @srrstats {G2.1} Assertion on type of input
  if (!is.logical(stochastic)) {
    stop2("Argument 'stochastic' must be logical.")
  }
  if (!is.null(intervention_values)) {
    # stop if stochastic and intervention
    if (stochastic == TRUE) {
      stop2(
        paste0(
          "Argument 'stochastic' cannot be `TRUE` when intervention_values",
          " are set."
        )
      )
    }
    # stop if intervention_values argument is not a named list
    #' @srrstats {G2.1} Assertion on type of input
    if (!is.list(intervention_values) || is.null(names(intervention_values))) {
      stop2("Argument 'intervention_values' is not a named list.")
    }
    # stop if intervention_values contains variables not included in the model
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
      stop2(
        "Argument 'intervention_values' contains duplicated variable names."
      )
    }
    # stop if any values in intervention_values are not of length one
    if (any(unlist(lapply(intervention_values, length)) != 1)) {
      stop2("Values in 'intervention_values' must each be of length one.")
    }
    # stop if any values in intervention_values are not NA or numeric
    if (any(unlist(lapply(intervention_values,
                          function(x) !is.na(x) & !is.numeric(x))))) {
      stop2("Values in 'intervention_values' must each be NA or numeric.")
    }
    # stop if all variables are held constant in intervention_values
    if (mean(is.na(intervention_values)) == 0) {
      stop2(
        paste0(
          "Argument 'intervention_values' must have at least one NA value ",
          "declaring a free variable. If all variables are held constant, the ",
          "system is already at equilibrium and there is nothing to compute."
        )
      )
    }
    # stop if all variables in intervention_values are free (all NA)
    if (mean(is.na(intervention_values)) == 1) {
      stop2(
        paste0(
          "Argument 'intervention_values' must have at least one variable ",
          "held constant (i.e., not all values are NA)."
        )
      )
    }
  }
}
