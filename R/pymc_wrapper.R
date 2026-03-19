#' Create a wrapper object for PyMC fits that mimics cmdstanr interface
#'
#' @param trace_result PyMC trace object (Python object via reticulate).
#' @param draws_array A draws_array object from the posterior package.
#' @param stan_variables Character vector of Stan-compatible variable names.
#' @param iter_sampling Integer. Number of sampling iterations per chain.
#' @param iter_warmup Integer. Number of warmup iterations per chain.
#' @param chains Integer. Number of chains.
#' @param seed Integer. Random seed used.
#'
#' @returns An object of class 'pymc_fit' that mimics cmdstanr's interface.
#'
#' @noRd
create_pymc_wrapper <- function(trace_result, draws_array, stan_variables,
                                iter_sampling, iter_warmup, chains,
                                seed = NULL) {
  wrapper <- list(
    trace_result = trace_result,
    draws_array = draws_array,
    stan_variables = stan_variables,
    iter_sampling = iter_sampling,
    iter_warmup = iter_warmup,
    chains = chains,
    seed = seed
  )
  class(wrapper) <- "pymc_fit"

  wrapper$draws <- function(variables = NULL, ...) {
    draws.pymc_fit(wrapper, variables = variables, ...)
  }
  wrapper$summary <- function(variables = NULL, ...) {
    summary.pymc_fit(wrapper, variables = variables, ...)
  }
  wrapper$metadata <- function() {
    metadata.pymc_fit(wrapper)
  }
  wrapper$num_chains <- function() {
    wrapper$chains
  }
  wrapper$diagnostic_summary <- function(diagnostic = "divergences",
                                         quiet = TRUE) {
    if (diagnostic == "divergences") {
      tryCatch({
        if (reticulate::py_has_attr(wrapper$trace_result, "sample_stats")) {
          sample_stats <- wrapper$trace_result$sample_stats
          if (!is.null(sample_stats) &&
                reticulate::py_has_attr(sample_stats, "diverging")) {
            div_arr <- reticulate::py_to_r(sample_stats$diverging$values)
            return(list(num_divergent = as.integer(sum(div_arr, na.rm = TRUE))))
          }
        } else if (reticulate::py_has_attr(wrapper$trace_result, "get_sampler_stats")) {
          div_arr <- reticulate::py_to_r(
            wrapper$trace_result$get_sampler_stats(
              "diverging",
              combine = TRUE,
              squeeze = TRUE
            )
          )
          return(list(num_divergent = as.integer(sum(div_arr, na.rm = TRUE))))
        }
        list(num_divergent = 0L)
      }, error = function(e) {
        if (!quiet) {
          warning("Could not extract divergence info from PyMC trace: ",
                  conditionMessage(e))
        }
        list(num_divergent = 0L)
      })
    } else {
      list()
    }
  }
  wrapper
}

#' Extract draws from pymc_fit object
#'
#' @param x A pymc_fit object.
#' @param variables Character vector of variable names to extract.
#' @param ... Additional arguments (ignored).
#'
#' @returns A draws_array object.
#'
#' @method draws pymc_fit
#' @export
draws.pymc_fit <- function(x, variables = NULL, ...) {
  if (is.null(variables)) {
    x$draws_array
  } else {
    posterior::subset_draws(x$draws_array, variable = variables)
  }
}

#' Summary statistics for pymc_fit object
#'
#' @param object A pymc_fit object.
#' @param variables Character vector of variable names to summarize.
#' @param ... Named summary spec arguments (cmdstanr style).
#'
#' @returns A data.frame with summary statistics.
#'
#' @method summary pymc_fit
#' @export
summary.pymc_fit <- function(object, variables = NULL, ...) {
  if (missing(variables) || is.null(variables)) {
    draws <- object$draws_array
  } else {
    draws <- posterior::subset_draws(object$draws_array, variable = variables)
  }

  summary_specs <- list(...)
  function_map <- list(
    "mean" = "mean", "median" = "median", "sd" = "sd",
    "mad" = "mad", "rhat" = "rhat",
    "ess_bulk" = "ess_bulk", "ess_tail" = "ess_tail"
  )

  posterior_args <- list()
  for (col_name in names(summary_specs)) {
    spec <- summary_specs[[col_name]]
    if (is.character(spec) && length(spec) == 1 &&
          spec %in% names(function_map)) {
      posterior_args[[col_name]] <- function_map[[spec]]
    } else {
      posterior_args[[col_name]] <- spec
    }
  }

  if (length(posterior_args) == 0) {
    posterior_args <- list(
      mean = "mean", sd = "sd", rhat = "rhat",
      ess_bulk = "ess_bulk", ess_tail = "ess_tail"
    )
  }

  summary_df <- do.call(
    posterior::summarise_draws,
    c(list(.x = draws), posterior_args)
  )
  if (!"variable" %in% names(summary_df)) {
    summary_df$variable <- rownames(summary_df)
  }
  summary_df
}

#' Extract metadata from pymc_fit object
#'
#' @param x A pymc_fit object.
#'
#' @returns A list containing metadata.
#'
#' @method metadata pymc_fit
#' @export
metadata.pymc_fit <- function(x) {
  list(
    stan_variables = x$stan_variables,
    model_params = x$stan_variables,
    iter_sampling = x$iter_sampling,
    iter_warmup = x$iter_warmup,
    num_chains = x$chains,
    chains = x$chains,
    thin = 1L,
    seed = x$seed
  )
}
