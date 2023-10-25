#' Create a summary of a fitted model represented by a \code{coevfit} object
#'
#' @param object An object of class \code{coevfit}.
#' @param prob A value between 0 and 1 indicating the desired probability
#' to be covered by the uncertainty intervals. The default is 0.95.
#' @param ... Other potential arguments
#'
#' @details The convergence diagnostics \code{rhat}, \code{ess_bulk}, and
#' \code{ess_tail} are described in detail in Vehtari et al. (2020).
#'
#' @references
#' Aki Vehtari, Andrew Gelman, Daniel Simpson, Bob Carpenter, and
#' Paul-Christian Bürkner (2020). Rank-normalization, folding, and
#' localization: An improved R-hat for assessing convergence of
#' MCMC. *Bayesian Analysis*. 1–28. dpi:10.1214/20-BA1221
#'
#' @method summary coevfit
#' @export
summary.coevfit <- function(object, prob = 0.95, ...) {
  # stop if prob is outside of range 0 - 1
  if (prob <= 0 | prob >= 1) {
    stop2("Argument 'prob' is not between 0 and 1.")
  }
  # get summary of selection and drift parameters from cmdstanr
  s <- as.data.frame(object$fit$summary(variables = c("alpha","sigma")))
  # summarise autoregressive selection effects
  auto <- s[grepl("alpha", s$variable) & (substring(s$variable, 7, 7) == substring(s$variable, 9, 9)),]
  rownames(auto) <- names(object$variables)[as.integer(substring(auto$variable, 7, 7))]
  auto <- auto[,2:ncol(auto)]
  # summarise cross selection effects
  cross <- s[grepl("alpha", s$variable) & (substring(s$variable, 7, 7) != substring(s$variable, 9, 9)),]
  rownames(cross) <-
    paste0(
      names(object$variables)[as.integer(substring(cross$variable, 9, 9))],
      " \U27F6 ",
      names(object$variables)[as.integer(substring(cross$variable, 7, 7))]
    )
  cross <- cross[,2:ncol(cross)]
  # summarise drift parameters
  drift <- s[grepl("sigma", s$variable),]
  rownames(drift) <- names(object$variables)[as.integer(substring(drift$variable, 7, 7))]
  drift <- drift[,2:ncol(drift)]
  # create summary list
  out <-
    list(
      variables     = object$variables,
      data_name     = object$data_name,
      nobs          = nrow(object$data),
      chains        = object$fit$num_chains(),
      iter          = object$fit$metadata()$iter_sampling,
      warmup        = object$fit$metadata()$iter_sampling,
      thin          = object$fit$metadata()$thin,
      auto          = auto,
      cross         = cross,
      drift         = drift,
      num_divergent = sum(object$fit$diagnostic_summary("divergences", quiet = TRUE)$num_divergent),
      rhats         = object$fit$summary(NULL, "rhat")$rhat
    )
  class(out) <- "coevsummary"
  return(out)
}

#' Print a summary of a fitted model represented by a \code{coevfit} object
#'
#' @aliases print.coevsummary
#'
#' @param x An object of class \code{coevfit}
#' @param digits The number of significant digits for printing out the summary; defaults to 2
#' @param ... Additional arguments that would be passed
#'  to method \code{summary} of \code{coevfit}.
#'
#' @seealso \code{\link{summary.coevfit}}
#'
#' @export
print.coevfit <- function(x, digits = 2, ...) {
  print(summary(x, ...), digits = digits, ...)
}

#' @export
print.coevsummary <- function(x, digits = 2, ...) {
  # print variables and response distributions
  cat("Variables: ")
  for (j in 1:length(x$variables)) {
    if (j == 1) {
      cat(names(x$variables)[j], "=", as.character(x$variables)[j], "\n")
    } else {
      cat("          ",
          names(x$variables)[j], "=", as.character(x$variables)[j], "\n")
    }
  }
  # print data name and number of observations
  cat(paste0("     Data: ", x$data_name, " (Number of observations: ", x$nobs, ")\n"))
  # print mcmc settings
  cat(paste0("    Draws: ", x$chains, " chains, each with iter = ",
             x$iter, "; warmup = ",
             x$warmup, "; thin = ",
             x$thin, "\n           total post-warmup draws = ",
             x$iter * x$chains, "\n\n"))
  # print autoregressive effects
  cat("Autoregressive selection effects:\n")
  print_format(x$auto, digits = digits)
  cat("\n")
  # print cross effects
  cat("Cross selection effects:\n")
  print_format(x$cross, digits = digits)
  cat("\n")
  # print drift
  cat("Drift scale parameters:\n")
  print_format(x$drift, digits = digits)
  cat("\n")
  # note
  cat("Note: Not all model parameters are displayed in this summary.")
  # warnings for high rhats or divergences
  if (max(x$rhats, na.rm = TRUE) > 1.05) {
    cat("\n")
    warning2(
      paste0(
        "Parts of the model have not converged (some Rhats are > 1.05). ",
        "Be careful when analysing the results! We recommend running ",
        "more iterations and/or setting stronger priors."
      )
    )
  }
  if (x$num_divergent > 0) {
    cat("\n")
    warning2(
      paste0(
        "There were ", x$num_divergent, " divergent transitions after warmup. ",
        "http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup"
      )
    )
  }
  invisible(x)
}

# helper function to print summary data frames in nice format
# also displays -0.00 as a result of round negative values to zero
# @param x object to be printed
# @param digits number of digits to show
# @param no_digits names of columns for which no digits should be shown
print_format <- function(x, digits = 2, no_digits = c("ess_bulk", "ess_tail")) {
  digits <- as.numeric(digits)
  if (length(digits) != 1L) {
    stop2("'digits' should be a single numeric value.")
  }
  out <- x
  fmt <- paste0("%.", digits, "f")
  for (i in 1:ncol(x)) {
    if (colnames(x)[i] %in% no_digits) {
      out[,i] <- sprintf("%.0f", x[,i])
    } else {
      out[,i] <- sprintf(fmt, x[,i])
    }
  }
  print(out, quote = FALSE, right = TRUE)
  invisible(x)
}
