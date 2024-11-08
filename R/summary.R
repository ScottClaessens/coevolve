#' Create a summary of a fitted model represented by a \code{coevfit} object
#'
#' @param object An object of class \code{coevfit}.
#' @param prob A value between 0 and 1 indicating the desired probability
#'   to be covered by the uncertainty intervals. The default is 0.95.
#' @param robust If \code{FALSE} (the default) the mean is used as the measure
#'   of central tendency and the standard deviation as the measure of
#'   variability. If \code{TRUE}, the median and the median absolute deviation
#'   are applied instead.
#' @param ... Other potential arguments
#'
#' @details The convergence diagnostics \code{rhat}, \code{ess_bulk}, and
#'   \code{ess_tail} are described in detail in Vehtari et al. (2020).
#'
#' @references
#' Aki Vehtari, Andrew Gelman, Daniel Simpson, Bob Carpenter, and
#' Paul-Christian Bürkner (2020). Rank-normalization, folding, and
#' localization: An improved R-hat for assessing convergence of
#' MCMC. *Bayesian Analysis*. 1–28. dpi:10.1214/20-BA1221
#'
#' @method summary coevfit
#' @export
summary.coevfit <- function(object, prob = 0.95, robust = FALSE, ...) {
  # stop if prob is outside of range 0 - 1
  if (prob <= 0 | prob >= 1) {
    stop2("Argument 'prob' is not between 0 and 1.")
  }
  # percentiles to print
  probs <- c(((1 - prob) / 2), 1 - ((1 - prob) / 2))
  # get summary of selection and drift parameters from cmdstanr
  s <-
    as.data.frame(
      object$fit$summary(
        NULL,
        Estimate = ifelse(robust, "median", "mean"),
        `Est.Error` = ifelse(robust, "mad", "sd"),
        ~quantile(.x, probs = probs),
        Rhat = "rhat",
        Bulk_ESS = "ess_bulk",
        Tail_ESS = "ess_tail"
        )
    )
  # summarise autoregressive selection effects
  A <- s[stringr::str_starts(s$variable, pattern = "A\\["),]
  equal <-
    readr::parse_number(
      stringr::str_extract(A$variable, pattern = "A\\[(\\d+)\\,")
      ) ==
    readr::parse_number(
      stringr::str_extract(A$variable, pattern = "\\,\\d+\\]")
      )
  auto <- A[equal,]
  rownames(auto) <- names(object$variables)[
    readr::parse_number(
      stringr::str_extract(auto$variable, pattern = "A\\[(\\d+)\\,")
      )
    ]
  auto <- auto[, 2:ncol(auto)]
  # summarise cross selection effects
  cross <- A[!equal,]
  rownames(cross) <-
    paste0(
      names(object$variables)[
        readr::parse_number(
          stringr::str_extract(cross$variable, pattern = "\\,\\d+\\]")
          )
        ],
      " \U27F6 ",
      names(object$variables)[
        readr::parse_number(
          stringr::str_extract(cross$variable, pattern = "A\\[(\\d+)\\,")
          )
        ]
    )
  cross <- cross[, 2:ncol(cross)]
  # only include cross selection effects that have been estimated in summary
  cross <- cross[!is.na(cross$Rhat),]
  # summarise drift sd parameters
  sd_drift <- s[stringr::str_starts(s$variable, pattern = "Q_sigma\\["),]
  rownames(sd_drift) <- paste0(
    "sd(",
    names(object$variables)[
      readr::parse_number(
        stringr::str_extract(sd_drift$variable, pattern = "Q_sigma\\[(\\d+)\\]")
      )
    ],
    ")"
  )
  sd_drift <- sd_drift[, 2:ncol(sd_drift)]
  # summarise drift cor parameters
  cor_drift <- NULL
  if (object$estimate_Q_offdiag) {
    cor_drift <- s[stringr::str_starts(s$variable, "cor_R"),]
    for (i in 1:length(object$variables)) {
      for (j in 1:length(object$variables)) {
        if (i >= j) {
          var <- paste0("cor_R[", i, ",", j, "]")
          cor_drift <- cor_drift[cor_drift$variable != var,]
        }
      }
    }
    rownames(cor_drift) <- paste0(
      "cor(",
      names(object$variables)[
        readr::parse_number(
          stringr::str_extract(
            cor_drift$variable,
            pattern = "cor\\_R\\[(\\d+)\\,"
          )
        )
      ],
      ",",
      names(object$variables)[
        readr::parse_number(
          stringr::str_extract(
            cor_drift$variable,
            pattern = "\\,(\\d+)\\]"
          )
        )
      ],
      ")"
    )
    cor_drift <- cor_drift[, 2:ncol(cor_drift)]
  }
  # summarise SDE intercepts
  sde_intercepts <- s[stringr::str_starts(s$variable, "b"),]
  rownames(sde_intercepts) <- names(object$variables)[
    readr::parse_number(sde_intercepts$variable)
    ]
  sde_intercepts <- sde_intercepts[, 2:ncol(sde_intercepts)]
  # summarise ordinal cutpoints
  cutpoints <- NULL
  if ("ordered_logistic" %in% object$variables) {
    cutpoints <- s[stringr::str_starts(s$variable, "c") &
                     !stringr::str_starts(s$variable, "cor_group") &
                     !stringr::str_starts(s$variable, "cor_R"),]
    rownames(cutpoints) <- paste0(
      names(object$variables)[readr::parse_number(cutpoints$variable)],
      stringr::str_extract(cutpoints$variable, pattern = "\\[\\d+\\]")
      )
    cutpoints <- cutpoints[, 2:ncol(cutpoints)]
  }
  # summarise overdispersion parameters
  phi <- NULL
  if ("negative_binomial_softplus" %in% object$variables) {
    phi <- s[stringr::str_starts(s$variable, "phi"),]
    rownames(phi) <- names(object$variables)[readr::parse_number(phi$variable)]
    phi <- phi[, 2:ncol(phi)]
  }
  # summarise shape parameters
  shape <- NULL
  if ("gamma_log" %in% object$variables) {
    shape <- s[stringr::str_starts(s$variable, "shape"),]
    rownames(shape) <- names(object$variables)[readr::parse_number(shape$variable)]
    shape <- shape[, 2:ncol(shape)]
  }
  # summarise gaussian process parameters
  gpterms <- NULL
  if (!is.null(object$dist_mat)) {
    gpterms <- s[stringr::str_starts(s$variable, "rho_dist") |
                   stringr::str_starts(s$variable, "sigma_dist"),]
    gpvars <- stringr::str_extract(gpterms$variable, pattern = "[^_]+")
    rownames(gpterms) <- paste0(
      ifelse(gpvars == "sigma", "sdgp", gpvars),
      "(", names(object$variables)[readr::parse_number(gpterms$variable)], ")"
    )
    gpterms <- gpterms[, 2:ncol(gpterms)]
  }
  # summarise group-level hyperparameters
  sd_group <- NULL
  cor_group <- NULL
  if (any(duplicated(object$data[,object$id]))) {
    # sd parameters
    sd_group <- s[stringr::str_starts(s$variable, "sigma_group"),]
    rownames(sd_group) <- paste0(
      "sd(",
      names(object$variables)[readr::parse_number(sd_group$variable)],
      ")"
      )
    sd_group <- sd_group[, 2:ncol(sd_group)]
    # correlation parameters
    cor_group <- s[stringr::str_starts(s$variable, "cor_group"),]
    for (i in 1:length(object$variables)) {
      for (j in 1:length(object$variables)) {
        if (i >= j) {
          var <- paste0("cor_group[", i, ",", j, "]")
          cor_group <- cor_group[cor_group$variable != var,]
        }
      }
    }
    rownames(cor_group) <- paste0(
      "cor(",
      names(object$variables)[
        readr::parse_number(
          stringr::str_extract(
            cor_group$variable,
            pattern = "cor\\_group\\[(\\d+)\\,"
          )
        )
      ],
      ",",
      names(object$variables)[
        readr::parse_number(
          stringr::str_extract(
            cor_group$variable,
            pattern = "\\,(\\d+)\\]"
          )
        )
      ],
      ")"
    )
    cor_group <- cor_group[, 2:ncol(cor_group)]
  }
  # create summary list
  out <-
    list(
      data           = object$data,
      variables      = object$variables,
      data_name      = object$data_name,
      nobs           = object$stan_data$N_obs,
      ntrees         = object$stan_data$N_tree,
      tree_name      = object$tree_name,
      chains         = object$fit$num_chains(),
      iter           = object$fit$metadata()$iter_sampling,
      warmup         = object$fit$metadata()$iter_warmup,
      thin           = object$fit$metadata()$thin,
      auto           = auto,
      cross          = cross,
      sd_drift       = sd_drift,
      cor_drift      = cor_drift,
      sde_intercepts = sde_intercepts,
      cutpoints      = cutpoints,
      phi            = phi,
      shape          = shape,
      gpterms        = gpterms,
      sd_group       = sd_group,
      cor_group      = cor_group,
      num_divergent  = sum(
        object$fit$diagnostic_summary("divergences", quiet = TRUE)$num_divergent
        ),
      rhats          = object$fit$summary(NULL, "rhat")$rhat
    )
  class(out) <- "coevsummary"
  return(out)
}

#' Print a summary of a fitted model represented by a \code{coevfit} object
#'
#' @aliases print.coevsummary
#'
#' @param x An object of class \code{coevfit}
#' @param digits The number of significant digits for printing out the summary;
#'   defaults to 2
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
  cat(paste0("     Data: ", x$data_name,
             " (Number of observations: ", x$nobs, ")\n"))
  cat(paste0("Phylogeny: ", x$tree_name,
             " (Number of trees: ", x$ntrees, ")\n"))
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
  # print cross effects (if any have been estimated)
  if (nrow(x$cross) > 0) {
    cat("Cross selection effects:\n")
    print_format(x$cross, digits = digits)
    cat("\n")
  }
  # print drift
  cat("Drift parameters:\n")
  print_format(rbind(x$sd_drift, x$cor_drift), digits = digits)
  cat("\n")
  # print SDE intercepts
  cat("Continuous time intercept parameters:\n")
  print_format(x$sde_intercepts, digits = digits)
  # print ordinal cutpoints
  if (!is.null(x$cutpoints)) {
    cat("\n")
    cat("Ordinal cutpoint parameters:\n")
    print_format(x$cutpoints, digits = digits)
  }
  # print overdispersion parameters
  if (!is.null(x$phi)) {
    cat("\n")
    cat("Overdispersion parameters:\n")
    print_format(x$phi, digits = digits)
  }
  # print shape parameters
  if (!is.null(x$shape)) {
    cat("\n")
    cat("Shape parameters:\n")
    print_format(x$shape, digits = digits)
  }
  # print gaussian process parameters
  if (!is.null(x$gpterms)) {
    cat("\n")
    cat("Gaussian Process parameters for distances:\n")
    print_format(x$gpterms, digits = digits)
  }
  # print group-level hyperparameters
  if (!is.null(x$sd_group) & !is.null(x$cor_group)) {
    cat("\n")
    cat("Group-level hyperparameters:\n")
    print_format(rbind(x$sd_group, x$cor_group), digits = digits)
  }
  # warnings for high rhats or divergences
  if (max(x$rhats, na.rm = TRUE) > 1.05) {
    warning2(
      paste0(
        "Parts of the model have not converged (some Rhats are > 1.05). ",
        "Be careful when analysing the results! We recommend running ",
        "more iterations and/or setting stronger priors."
      )
    )
    cat("\n")
  }
  if (x$num_divergent > 0) {
    warning2(
      paste0(
        "There were ", x$num_divergent, " divergent transitions after warmup. ",
        "http://mc-stan.org/misc/warnings.html",
        "#divergent-transitions-after-warmup"
      )
    )
    cat("\n")
  }
  # warnings if data excluded
  if (nrow(x$data) != x$nobs) {
    warning2(
      paste0(
        "Rows with NAs for all coevolving variables were ",
        "excluded from the model."
      )
    )
    cat("\n")
  }
  # return
  invisible(x)
}

# helper function to print summary data frames in nice format
# also displays -0.00 as a result of round negative values to zero
# @param x object to be printed
# @param digits number of digits to show
# @param no_digits names of columns for which no digits should be shown
print_format <- function(x, digits = 2, no_digits = c("Bulk_ESS", "Tail_ESS")) {
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
