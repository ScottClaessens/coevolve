#' Create a wrapper object for nutpie fits that mimics cmdstanr interface
#'
#' @srrstats {G1.4a} Non-exported function documented here
#'
#' @description Creates an S3 object that wraps a nutpie trace and provides
#'   methods compatible with cmdstanr's fit object interface. This allows
#'   nutpie-fitted models to work with existing coevolve functions.
#'
#' @param trace Nutpie trace object (Python object via reticulate).
#' @param draws_array A draws_array object from the posterior package.
#' @param stan_variables Character vector of Stan variable names.
#' @param iter_sampling Integer. Number of sampling iterations per chain.
#' @param iter_warmup Integer. Number of warmup iterations per chain.
#' @param chains Integer. Number of chains.
#' @param seed Integer. Random seed used.
#'
#' @returns An object of class 'nutpie_fit' that mimics cmdstanr's interface.
#'
#' @noRd
create_nutpie_wrapper <- function(trace, draws_array, stan_variables,
                                  iter_sampling, iter_warmup, chains,
                                  seed = NULL) {
  # create wrapper object
  wrapper <- list(
    trace = trace,
    draws_array = draws_array,
    stan_variables = stan_variables,
    iter_sampling = iter_sampling,
    iter_warmup = iter_warmup,
    chains = chains,
    seed = seed
  )
  class(wrapper) <- "nutpie_fit"
  # add $draws() method for compatibility with cmdstanr interface
  # this allows code like object$fit$draws() to work for both samplers
  wrapper$draws <- function(variables = NULL, ...) {
    draws.nutpie_fit(wrapper, variables = variables, ...)
  }
  # add $summary() method for compatibility with cmdstanr interface
  # this allows code like object$fit$summary() to work for both samplers
  wrapper$summary <- function(variables = NULL, ...) {
    summary.nutpie_fit(wrapper, variables = variables, ...)
  }
  # add $metadata() method for compatibility with cmdstanr interface
  # this allows code like object$fit$metadata() to work for both samplers
  wrapper$metadata <- function() {
    metadata.nutpie_fit(wrapper)
  }
  # add $num_chains() method for compatibility with cmdstanr interface
  wrapper$num_chains <- function() {
    wrapper$chains
  }
  # add $diagnostic_summary() method for compatibility with cmdstanr interface
  # nutpie stores divergence information in trace$sample_stats$diverging
  wrapper$diagnostic_summary <- function(diagnostic = "divergences",
                                         quiet = TRUE) {
    if (diagnostic == "divergences") {
      # extract divergence information from nutpie InferenceData
      # divergences are stored in trace$sample_stats$diverging
      tryCatch({
        if (reticulate::py_has_attr(wrapper$trace, "sample_stats")) {
          sample_stats <- wrapper$trace$sample_stats
          if (reticulate::py_has_attr(sample_stats, "data_vars")) {
            # check if diverging variable exists
            py_builtins <- reticulate::import_builtins()
            keys_obj <- sample_stats$data_vars$keys()
            var_names <- reticulate::py_to_r(py_builtins$list(keys_obj))
            if ("diverging" %in% var_names) {
              # extract diverging array
              diverging_xarray <- sample_stats[["diverging"]]
              diverging_array <- diverging_xarray$values
              diverging_r <- reticulate::py_to_r(diverging_array)
              # count total number of divergences across all chains and draws
              # diverging_r should be [chains, draws] or similar
              num_divergent <- sum(diverging_r, na.rm = TRUE)
              return(list(num_divergent = as.integer(num_divergent)))
            }
          }
        }
        # if we can't find divergence info, return 0
        return(list(num_divergent = 0L))
      }, error = function(e) {
        # if there's an error extracting divergence info, return 0
        if (!quiet) {
          warning(
            "Could not extract divergence information from nutpie trace: ",
            conditionMessage(e)
          )
        }
        list(num_divergent = 0L)
      })
    } else {
      # for other diagnostics, return empty structure
      return(list())
    }
  }
  wrapper
}

#' Extract draws from nutpie_fit object
#'
#' @srrstats {G1.4a} Non-exported function documented here
#'
#' @description Extracts draws from a nutpie_fit object, compatible with
#'   cmdstanr's draws() method.
#'
#' @param x A nutpie_fit object.
#' @param variables Character vector of variable names to extract. If NULL,
#'   extracts all variables.
#' @param ... Additional arguments (currently ignored).
#'
#' @returns A draws_array object.
#'
#' @method draws nutpie_fit
#' @export
draws.nutpie_fit <- function(x, variables = NULL, ...) {
  if (is.null(variables)) {
    return(x$draws_array)
  } else {
    # Subset to requested variables
    return(posterior::subset_draws(x$draws_array, variable = variables))
  }
}

#' Compute summary statistics for nutpie_fit object
#'
#' @srrstats {G1.4a} Non-exported function documented here
#'
#' @description Computes summary statistics from a nutpie_fit object,
#'   compatible with cmdstanr's summary() method. This method mimics
#'   cmdstanr's interface where summary() is called with named arguments
#'   where names are column names and values are function names or formulas.
#'
#' @param object A nutpie_fit object.
#' @param variables Character vector of variable names to summarize. If NULL,
#'   summarizes all variables.
#' @param ... Named arguments where names are output column names and values
#'   are function names (strings) or formulas for computing statistics.
#'
#' @returns A data.frame with summary statistics, compatible with cmdstanr
#'   format.
#'
#' @method summary nutpie_fit
#' @export
summary.nutpie_fit <- function(object, variables = NULL, ...) {
  # handle cmdstanr-style call where first argument might be NULL
  # if variables is explicitly NULL or missing, use all variables
  if (missing(variables) || is.null(variables)) {
    draws <- object$draws_array
  } else {
    draws <- posterior::subset_draws(object$draws_array, variable = variables)
  }
  # extract summary arguments from ...
  # cmdstanr format: column_name = "function_name" or column_name = ~formula
  summary_specs <- list(...)
  # map cmdstanr function names to posterior function names (as strings)
  # posterior::default_summary_measures() returns character vectors, not
  # functions
  function_map <- list(
    "mean" = "mean",
    "median" = "median",
    "sd" = "sd",
    "mad" = "mad",
    "rhat" = "rhat",
    "ess_bulk" = "ess_bulk",
    "ess_tail" = "ess_tail"
  )
  # convert summary specs to posterior format
  posterior_args <- list()
  for (col_name in names(summary_specs)) {
    spec <- summary_specs[[col_name]]
    if (is.character(spec) && length(spec) == 1) {
      # string function name - map to posterior function name
      if (spec %in% names(function_map)) {
        posterior_args[[col_name]] <- function_map[[spec]]
      } else {
        # unknown function name - try to use as-is
        posterior_args[[col_name]] <- spec
      }
    } else if (inherits(spec, "formula")) {
      # formula - use as-is
      posterior_args[[col_name]] <- spec
    } else {
      # other - use as-is
      posterior_args[[col_name]] <- spec
    }
  }
  # if no arguments provided, use default summary statistics
  if (length(posterior_args) == 0) {
    posterior_args <- list(
      mean = "mean",
      sd = "sd",
      rhat = "rhat",
      ess_bulk = "ess_bulk",
      ess_tail = "ess_tail"
    )
  }
  # compute summary statistics using posterior package
  summary_df <- do.call(
    posterior::summarise_draws,
    c(list(.x = draws), posterior_args)
  )
  # ensure 'variable' column exists (posterior uses this by default)
  if (!"variable" %in% names(summary_df)) {
    # if variable column doesn't exist, it means we have a single variable
    # or the format is different - this shouldn't happen but handle it
    summary_df$variable <- rownames(summary_df)
  }
  summary_df
}

#' Extract metadata from nutpie_fit object
#'
#' @srrstats {G1.4a} Non-exported function documented here
#'
#' @description Extracts metadata from a nutpie_fit object, compatible with
#'   cmdstanr's metadata() method.
#'
#' @param x A nutpie_fit object.
#'
#' @returns A list containing metadata.
#'
#' @method metadata nutpie_fit
#' @export
metadata.nutpie_fit <- function(x) {
  # extract model parameters (all variables except generated quantities)
  # for now, we'll include all variables
  # in the future, we might want to distinguish parameters from generated
  # quantities
  model_params <- x$stan_variables
  # return metadata list compatible with cmdstanr
  # cmdstanr uses num_chains, not chains
  list(
    stan_variables = x$stan_variables,
    model_params = model_params,
    iter_sampling = x$iter_sampling,
    iter_warmup = x$iter_warmup,
    num_chains = x$chains,  # cmdstanr uses num_chains
    chains = x$chains,  # also include chains for compatibility
    thin = 1L,  # nutpie doesn't use thinning by default
    seed = x$seed
  )
}

#' Generic method for extracting draws
#'
#' @srrstats {G1.4} Function is documented
#'
#' @description Generic function for extracting draws from fit objects.
#'   Currently supports nutpie_fit objects.
#'
#' @param x A fit object.
#' @param ... Additional arguments passed to methods.
#'
#' @returns A draws_array object.
#'
#' @noRd
draws <- function(x, ...) {
  UseMethod("draws")
}

#' Generic method for extracting metadata
#'
#' @srrstats {G1.4} Function is documented
#'
#' @description Generic function for extracting metadata from fit objects.
#'   Currently supports nutpie_fit objects.
#'
#' @param x A fit object.
#'
#' @returns A list containing metadata.
#'
#' @noRd
metadata <- function(x) {
  UseMethod("metadata")
}
