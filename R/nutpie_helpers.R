#' Check if nutpie is available via reticulate
#'
#' @srrstats {G1.4} Function is documented
#'
#' @description Checks if nutpie can be imported via reticulate. This function
#'   attempts to import nutpie using the currently configured Python environment
#'   (via reticulate). Users must configure reticulate to use a Python
#'   environment with nutpie installed before using this function.
#'
#' @details To configure nutpie, you must:
#' \enumerate{
#'   \item Install nutpie in a Python environment: \code{pip install
#'     "nutpie[stan]"}
#'   \item Configure reticulate to use that Python environment:
#'     \itemize{
#'       \item Set \code{RETICULATE_PYTHON} environment variable to the Python
#'         executable path
#'       \item Or use \code{reticulate::use_python()} or
#'         \code{reticulate::use_virtualenv()}
#'     }
#' }
#'
#' @returns Logical. TRUE if nutpie is available, FALSE otherwise.
#'
#' @examples
#' \dontrun{
#' # Configure reticulate first (choose one):
#' # Option 1: Set environment variable
#' Sys.setenv(RETICULATE_PYTHON = "/path/to/python")
#'
#' # Option 2: Use Python directly
#' reticulate::use_python("/path/to/python")
#'
#' # Option 3: Use virtual environment
#' reticulate::use_virtualenv("~/.venvs/nutpie-env")
#'
#' # Then check availability
#' if (check_nutpie_available()) {
#'   # Use nutpie
#' } else {
#'   # Fall back to cmdstanr
#' }
#' }
#'
#' @export
check_nutpie_available <- function() {
  # check if reticulate is available
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    return(FALSE)
  }
  # try to import nutpie using the currently configured Python
  # no automatic searching - user must configure reticulate explicitly
  tryCatch({
    reticulate::import("nutpie", convert = FALSE)
    # If we get here, nutpie is available
    TRUE
  }, error = function(e) {
    # nutpie not available
    FALSE
  })
}

#' Compile a Stan model using nutpie
#'
#' @srrstats {G1.4} Function is documented
#'
#' @description Compiles a Stan model string using nutpie's compile_stan_model
#'   function. This is a wrapper around nutpie's Python API via reticulate.
#'
#' @param stan_code Character string containing Stan model code.
#' @param ... Additional arguments passed to \code{nutpie.compile_stan_model()}.
#'   Valid parameters include:
#'   \itemize{
#'     \item \code{extra_stanc_args}: List of strings. Arguments passed to
#'       stanc3 (Stan compiler). Common options include \code{"--O0"} (no
#'       optimization), \code{"--O1"} (optimization level 1). Example:
#'       \code{extra_stanc_args = list("--O1")}.
#'     \item \code{extra_compile_args}: List of strings. Arguments passed to
#'       Make (C++ compiler). Common options include \code{"STAN_THREADS=true"}
#'       to enable threading. Example: \code{extra_compile_args =
#'       list("STAN_THREADS=true")}.
#'     \item \code{filename}: Character string. Path to Stan file (alternative
#'       to \code{stan_code}).
#'     \item \code{dims}: Named list. Dimension information.
#'     \item \code{coords}: Named list. Coordinate information.
#'     \item \code{model_name}: Character string. Name for the compiled model.
#'     \item \code{cleanup}: Logical. Whether to clean up temporary files
#'       (default: TRUE).
#'   }
#'
#' @returns Compiled Stan model object (Python object via reticulate).
#'
#' @examples
#' \dontrun{
#' stan_code <- "
#' data {
#'   int<lower=0> N;
#'   vector[N] y;
#' }
#' parameters {
#'   real mu;
#' }
#' model {
#'   mu ~ normal(0, 1);
#'   y ~ normal(mu, 1);
#' }
#' "
#' # compile Stan model with default settings
#' compiled <- nutpie_compile_stan_model(stan_code)
#'
#' # compile with optimization level 1
#' compiled_o1 <- nutpie_compile_stan_model(
#'   stan_code,
#'   extra_stanc_args = list("--O1")
#' )
#'
#' # compile with threading enabled
#' compiled_threads <- nutpie_compile_stan_model(
#'   stan_code,
#'   extra_compile_args = list("STAN_THREADS=true")
#' )
#' }
#'
#' @export
nutpie_compile_stan_model <- function(stan_code, ...) {
  # check if nutpie is available
  if (!check_nutpie_available()) {
    stop2(
      "nutpie is not available. Please install nutpie with: ",
      "pip install 'nutpie[stan]'"
    )
  }
  # import nutpie
  nutpie <- reticulate::import("nutpie", convert = FALSE)
  # prepare compilation arguments
  compile_args <- list(code = stan_code)
  # add any additional arguments from ...
  additional_args <- list(...)
  compile_args <- c(compile_args, additional_args)
  # compile Stan model
  tryCatch({
    do.call(nutpie$compile_stan_model, compile_args)
  }, error = function(e) {
    # convert Python error to R error
    error_msg <- conditionMessage(e)
    stop2(
      "Failed to compile Stan model with nutpie: ",
      error_msg
    )
  })
}

#' Convert R data to Python format for nutpie
#'
#' @srrstats {G1.4a} Non-exported function documented here
#'
#' @description Converts an R list (Stan data) to Python format suitable for
#'   nutpie. Handles matrices, arrays, integers, and numeric vectors.
#'
#' @param data_list Named list containing Stan data.
#'
#' @returns Python dict (via reticulate) with converted data.
#'
#' @noRd
convert_r_to_python_data <- function(data_list) {
  # import numpy for array conversion
  np <- reticulate::import("numpy", convert = FALSE)
  # convert each element of the data list
  python_data <- list()
  for (name in names(data_list)) {
    value <- data_list[[name]]
    # handle different R types
    if (is.matrix(value) || is.array(value)) {
      # check if matrix/array should be integer type
      # Stan integer variables: tip, node_seq, parent, effects_mat, tip_id, etc.
      integer_vars <- c("tip", "node_seq", "parent", "effects_mat", "tip_id",
                        "N_tips", "N_tree", "N_obs", "J", "N_seg",
                        "num_effects", "prior_only", "miss")
      # check if this is a known integer variable or if values are integers
      is_known_integer <- name %in% integer_vars
      is_integer_type <- is.integer(value)
      # check if all non-NA values are integers (within tolerance)
      value_no_na <- value[!is.na(value)]
      all_integers <- if (length(value_no_na) > 0) {
        all(abs(value_no_na - round(value_no_na)) < 1e-10)
      } else {
        TRUE  # empty or all NA - treat as integer if known integer var
      }
      # convert to integer if it's a known integer variable
      # or if all values are integers
      if (is_known_integer || (is_integer_type && all_integers)) {
        # convert to integer array, preserving matrix/array structure
        # don't use as.integer() directly as it may flatten matrices
        # instead, convert element-wise while preserving dimensions
        value_int <- array(as.integer(value), dim = dim(value))
        # NAs become NA_integer_ which numpy will handle
        # for Stan arrays like array[N_tree, N_seg], ensure 2D structure is
        # preserved even when N_tree=1 (numpy might squeeze single dimensions)
        if (name %in% c("node_seq", "parent", "ts", "tip") &&
              length(dim(value_int)) == 2) {
          # use ndmin=2 to prevent numpy from squeezing dimensions
          python_data[[name]] <- np$array(value_int, dtype = "int32",
                                          ndmin = 2L)
        } else {
          python_data[[name]] <- np$array(value_int, dtype = "int32")
        }
      } else {
        # convert to float array, preserving structure
        # for multi-dim arrays like ts, ensure 2D structure is preserved
        if (name %in% c("ts") && length(dim(value)) == 2) {
          python_data[[name]] <- np$array(array(as.double(value),
                                                dim = dim(value)),
                                          dtype = "float64", ndmin = 2L)
        } else {
          python_data[[name]] <- np$array(array(as.double(value),
                                                dim = dim(value)),
                                          dtype = "float64")
        }
      }
    } else if (is.integer(value) || (is.numeric(value) &&
                                       all(abs(value - round(value)) < 1e-10,
                                           na.rm = TRUE))) {
      # convert to Python integer or integer array
      if (length(value) == 1) {
        python_data[[name]] <- as.integer(value)
      } else {
        python_data[[name]] <- np$array(as.integer(value), dtype = "int32")
      }
    } else if (is.numeric(value)) {
      # convert to Python float or float array
      if (length(value) == 1) {
        python_data[[name]] <- as.double(value)
      } else {
        python_data[[name]] <- np$array(as.double(value), dtype = "float64")
      }
    } else {
      # try to convert as-is (might be list/other)
      python_data[[name]] <- value
    }
  }
  # convert to Python dict
  reticulate::r_to_py(python_data)
}

#' Sample from a Stan model using nutpie
#'
#' @srrstats {G1.4} Function is documented
#'
#' @description Samples from a compiled Stan model using nutpie. This function
#'   handles compilation, data conversion, and sampling.
#'
#' @param stan_code Character string containing Stan model code, OR a compiled
#'   nutpie model object.
#' @param data_list Named list containing Stan data.
#' @param num_chains Integer. Number of chains to run.
#' @param num_samples Integer. Number of post-warmup samples per chain.
#' @param num_warmup Integer. Number of warmup samples per chain.
#' @param seed Integer. Random seed for reproducibility.
#' @param target_accept Numeric between 0 and 1. Target acceptance probability
#'   for the NUTS sampler (default: NULL, uses nutpie default of 0.8). This
#'   corresponds to Stan's \code{adapt_delta} parameter.
#' @param low_rank_modified_mass_matrix Logical. If \code{TRUE}, enables
#'   low-rank modified mass matrix adaptation for models with strong parameter
#'   correlations (default: \code{FALSE}). This is an experimental feature in
#'   nutpie. When enabled, consider increasing \code{target_accept} (e.g., 0.95
#'   or 0.99) and adjusting \code{mass_matrix_gamma} and
#'   \code{mass_matrix_eigval_cutoff} via \code{...} to reduce divergences.
#' @param ... Additional arguments passed to \code{nutpie.sample()} or
#'   \code{nutpie.compile_stan_model()} (when \code{stan_code} is a character
#'   string). Useful parameters when using \code{low_rank_modified_mass_matrix =
#'   TRUE} include: \code{mass_matrix_gamma} (regularization parameter, default
#'   ~0.05) and \code{mass_matrix_eigval_cutoff} (eigenvalue cutoff, default
#'   ~0.01). For compilation, valid parameters include:
#'   \itemize{
#'     \item \code{extra_stanc_args}: List of strings. Arguments for stanc3
#'       compiler. Example: \code{extra_stanc_args = list("--O1")} for
#'       optimization level 1.
#'     \item \code{extra_compile_args}: List of strings. Arguments for Make/C++
#'       compiler. Example: \code{extra_compile_args =
#'       list("STAN_THREADS=true")} to enable threading.
#'     \item Other compilation parameters: \code{filename}, \code{dims},
#'       \code{coords}, \code{model_name}, \code{cleanup}.
#'   }
#'
#' @returns Nutpie trace object (Python object via reticulate).
#'
#' @examples
#' \dontrun{
#' stan_code <- "
#' data {
#'   int<lower=0> N;
#'   vector[N] y;
#' }
#' parameters {
#'   real mu;
#' }
#' model {
#'   mu ~ normal(0, 1);
#'   y ~ normal(mu, 1);
#' }
#' "
#' data_list <- list(N = 10L, y = rnorm(10))
#' trace <- nutpie_sample(stan_code, data_list,
#'                       num_chains = 4L,
#'                       num_samples = 1000L,
#'                       num_warmup = 500L,
#'                       seed = 12345L)
#' }
#'
#' @export
nutpie_sample <- function(stan_code, data_list,
                          num_chains = 4L,
                          num_samples = 1000L,
                          num_warmup = 500L,
                          seed = NULL,
                          target_accept = NULL,
                          low_rank_modified_mass_matrix = FALSE,
                          ...) {
  # check if nutpie is available
  if (!check_nutpie_available()) {
    stop2(
      "nutpie is not available. Please install nutpie with: ",
      "pip install 'nutpie[stan]'"
    )
  }
  # import nutpie
  nutpie <- reticulate::import("nutpie", convert = FALSE)
  # separate compilation arguments from sampling arguments
  all_additional_args <- list(...)
  # known compilation arguments for compile_stan_model
  # these are the valid parameters for nutpie.compile_stan_model
  compile_param_names <- c("filename", "extra_compile_args", "extra_stanc_args",
                           "dims", "coords", "model_name", "cleanup")
  # split args into compilation and sampling args
  compile_args <- list()
  sampling_args <- list()
  for (arg_name in names(all_additional_args)) {
    if (arg_name %in% compile_param_names) {
      compile_args[[arg_name]] <- all_additional_args[[arg_name]]
    } else {
      # assume it's a sampling argument
      sampling_args[[arg_name]] <- all_additional_args[[arg_name]]
    }
  }
  # compile model if stan_code is a character string
  if (is.character(stan_code)) {
    # pass compilation args to compile_stan_model
    compiled <- do.call(
      nutpie_compile_stan_model,
      c(list(stan_code = stan_code), compile_args)
    )
  } else {
    # assume it's already compiled
    compiled <- stan_code
  }
  # convert R data to Python format
  python_data <- convert_r_to_python_data(data_list)
  # prepare sampling arguments with correct nutpie parameter names
  # nutpie uses: draws, tune, chains (not num_samples, num_warmup, num_chains)
  # based on docstring: "draws: The number of draws after tuning in each chain"
  # and "tune: The number of tuning (warmup) draws in each chain"
  # so both should be per-chain values, matching cmdstanr's behavior.
  sample_args <- list(
    chains = as.integer(num_chains),
    draws = as.integer(num_samples),
    tune = as.integer(num_warmup)
  )
  # add seed if provided
  if (!is.null(seed)) {
    sample_args$seed <- as.integer(seed)
  }
  # add target_accept if provided (maps from adapt_delta in Stan)
  if (!is.null(target_accept)) {
    sample_args$target_accept <- as.double(target_accept)
  }
  # add low_rank_modified_mass_matrix if provided
  if (!is.null(low_rank_modified_mass_matrix)) {
    sample_args$low_rank_modified_mass_matrix <-
      as.logical(low_rank_modified_mass_matrix)
  }
  # nutpie returns InferenceData by default (not raw trace)
  # we'll extract draws from trace$posterior in convert_nutpie_draws()
  # only set return_raw_trace if user explicitly requests it
  # for sampling, use the sampling_args we separated out
  additional_args <- sampling_args
  # ensure tune, draws, and chains are not overridden by additional_args
  # these are core parameters that should come from the function arguments
  if ("tune" %in% names(additional_args)) {
    additional_args$tune <- NULL
  }
  if ("draws" %in% names(additional_args)) {
    additional_args$draws <- NULL
  }
  if ("chains" %in% names(additional_args)) {
    additional_args$chains <- NULL
  }
  if (!"return_raw_trace" %in% names(additional_args)) {
    # default: use InferenceData (return_raw_trace = FALSE)
    sample_args <- c(sample_args, additional_args)
  } else {
    # user explicitly set return_raw_trace, use their value
    sample_args$return_raw_trace <- additional_args$return_raw_trace
    additional_args$return_raw_trace <- NULL
    sample_args <- c(sample_args, additional_args)
  }
  # prepare data using with_data method
  tryCatch({
    # nutpie uses with_data() method to attach data
    # in Python: compiled.with_data(**data_dict)
    # in R via reticulate, we need to unpack the dict as keyword arguments
    # we'll define a Python helper function to handle **kwargs
    # define Python helper function to handle **kwargs unpacking
    # this is needed because reticulate doesn't directly support **kwargs syntax
    # we define it each time (it's lightweight and ensures it's available)
    reticulate::py_run_string("
def nutpie_call_with_data(compiled_model, data_dict):
    return compiled_model.with_data(**data_dict)
")
    # call the helper function
    call_func <- reticulate::py$nutpie_call_with_data
    model_with_data <- call_func(compiled, python_data)
    # sample
    trace <- do.call(nutpie$sample, c(list(model_with_data), sample_args))
    return(trace)
  }, error = function(e) {
    # Convert Python error to R error
    error_msg <- conditionMessage(e)
    stop2(
      "Failed to sample from Stan model with nutpie: ",
      error_msg
    )
  })
}

#' Convert nutpie draws to posterior draws_array format
#'
#' @srrstats {G1.4} Function is documented
#'
#' @description Converts draws from a nutpie trace object to a posterior
#'   draws_array format compatible with cmdstanr output. This function handles
#'   dimension ordering and variable extraction.
#'
#' @param trace Nutpie trace object (Python object via reticulate).
#'
#' @returns A draws_array object from the posterior package.
#'
#' @examples
#' \dontrun{
#' # after sampling with nutpie
#' trace <- nutpie_sample(stan_code, data_list)
#' draws <- convert_nutpie_draws(trace)
#' }
#'
#' @export
convert_nutpie_draws <- function(trace) {
  # check if posterior package is available
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop2("The 'posterior' package is required for draw conversion.")
  }
  # extract draws from nutpie InferenceData object
  # nutpie returns InferenceData (ArviZ) with posterior as xarray Dataset
  # format: trace$posterior$data_vars contains variables with dims (chain,
  # draw, ...)
  tryCatch({
    # check if this is InferenceData (has posterior attribute)
    if (reticulate::py_has_attr(trace, "posterior")) {
      posterior_ds <- trace$posterior
      # get variable names from data_vars
      # keys() returns a Python dict_keys object,
      # convert to list using Python's list()
      keys_obj <- posterior_ds$data_vars$keys()
      py_builtins <- reticulate::import_builtins()
      var_names <- reticulate::py_to_r(py_builtins$list(keys_obj))
    } else {
      # try raw trace with draws() method (if return_raw_trace=TRUE was used)
      draws_dict <- trace$draws()
      var_names <- names(reticulate::py_to_r(draws_dict))
    }
  }, error = function(e) {
    stop2("Failed to extract draws from nutpie trace: ", conditionMessage(e))
  })
  if (length(var_names) == 0) {
    stop2("No variables found in nutpie trace.")
  }
  # process each variable
  # cmdstanr flattens all parameter dimensions into separate scalar variables
  # e.g., z[1,2,3] becomes separate variables z[1,2,3]
  # so we need to flatten arrays and create separate variables for each element
  draws_arrays <- list()
  for (var_name in var_names) {
    tryCatch({
      if (reticulate::py_has_attr(trace, "posterior")) {
        # extract from InferenceData posterior Dataset
        var_xarray <- trace$posterior[[var_name]]
        # convert xarray DataArray to numpy array
        # xarray format: (chain, draw, ...) -> numpy array
        var_array <- var_xarray$values
        var_r <- reticulate::py_to_r(var_array)
        # ensure it's fully converted (not a Python reference)
        if (inherits(var_r, "python.builtin.object")) {
          var_r <- reticulate::py_to_r(var_r)
        }
      } else {
        # extract from raw trace draws dict
        draws_dict <- trace$draws()
        draws_list <- reticulate::py_to_r(draws_dict)
        var_array <- draws_list[[var_name]]
        var_r <- reticulate::py_to_r(var_array)
        # ensure it's fully converted (not a Python reference)
        if (inherits(var_r, "python.builtin.object")) {
          var_r <- reticulate::py_to_r(var_r)
        }
      }
    }, error = function(e) {
      stop2("Failed to convert variable '", var_name, "' to R format: ",
            conditionMessage(e))
    })
    # check for NULL or empty variables
    if (is.null(var_r)) {
      stop2("Variable '", var_name, "' is NULL after conversion from nutpie.")
    }
    if (length(var_r) == 0) {
      stop2("Variable '", var_name, "' is empty (length 0) after conversion from nutpie.")
    }
    # nutpie/ArviZ format: [chains, draws, dims...]
    # posterior draws_array format: [iterations, chains, dims...]
    # so we need to swap first two dimensions:
    # [chains, draws, ...] -> [draws, chains, ...]
    if (is.array(var_r)) {
      dims_var <- dim(var_r)
      n_dims <- length(dims_var)
      if (n_dims == 1) {
        # 1D array: unexpected format from nutpie
        # nutpie should return [chains, draws, ...] format
        # If we get 1D, assume it's malformed and reshape to [length, 1]
        # representing [draws=length, chains=1]
        # This will be caught by dimension validation if inconsistent
        var_length <- length(var_r)
        var_r <- array(var_r, dim = c(var_length, 1))
        draws_arrays[[var_name]] <- var_r
      } else if (n_dims == 2) {
        # scalar variable: [chains, draws] -> [draws, chains]
        var_r <- aperm(var_r, c(2, 1))
        # ensure it's still a 2D array (R might drop dimensions)
        if (is.null(dim(var_r)) || length(dim(var_r)) != 2) {
          var_r <- array(var_r, dim = dims_var[c(2, 1)])
        }
        # add as single variable
        draws_arrays[[var_name]] <- var_r
      } else if (n_dims > 2) {
        # array variable: [chains, draws, dims...] -> [draws, chains, dims...]
        perm <- c(2, 1, 3:n_dims)
        var_r <- aperm(var_r, perm)
        # get dimensions after permutation
        dims_permuted <- dim(var_r)
        n_draws_perm <- dims_permuted[1]
        n_chains_perm <- dims_permuted[2]
        # flatten parameter dimensions into separate variables
        # reshape to [iterations, chains, param_elements]
        param_dims <- dims_permuted[3:n_dims]
        n_param_elements <- prod(param_dims)
        var_reshaped <- array(var_r, dim = c(n_draws_perm, n_chains_perm,
                                             n_param_elements))
        # create separate variable for each parameter element
        for (i in seq_len(n_param_elements)) {
          # get indices for this element
          indices <- arrayInd(i, param_dims)
          # create variable name with indices
          if (length(indices) == 1) {
            var_name_indexed <- paste0(var_name, "[", indices[1], "]")
          } else {
            var_name_indexed <- paste0(var_name, "[",
                                       paste(indices, collapse = ","), "]")
          }
          # extract this element:
          # [iterations, chains, element] -> [iterations, chains]
          # R will drop dimensions when extracting, so we need to ensure
          # the result is always a 2D array [draws, chains]
          extracted <- var_reshaped[, , i]
          # check dimensions and restore if dropped
          extracted_dims <- dim(extracted)
          if (is.null(extracted_dims)) {
            # dimension was dropped - it's now a vector
            # restore as [draws, chains=1] or [draws=1, chains]
            # we know original shape was [n_draws_perm, n_chains_perm]
            extracted <- array(extracted, dim = c(n_draws_perm, n_chains_perm))
          } else if (length(extracted_dims) == 1) {
            # one dimension was dropped
            # restore as [draws, chains=1] or [draws=1, chains]
            extracted <- array(extracted, dim = c(n_draws_perm, n_chains_perm))
          } else if (length(extracted_dims) != 2) {
            # unexpected number of dimensions
            stop2(
              "Unexpected dimensions when extracting variable '", var_name_indexed,
              "'. Expected 2D array [draws, chains], got ", length(extracted_dims),
              " dimensions: [", paste(extracted_dims, collapse = ", "), "]"
            )
          }
          draws_arrays[[var_name_indexed]] <- extracted
        }
      }
    } else {
      # not an array - convert to array format
      # This shouldn't happen with nutpie (should always return arrays)
      # but handle gracefully for debugging
      tryCatch({
        # try to convert to numeric array
        var_numeric <- as.numeric(var_r)
        if (any(is.na(var_numeric)) && !is.null(var_r)) {
          stop2(
            "Variable '", var_name, "' could not be converted to numeric. ",
            "Type: ", typeof(var_r), ", Class: ", paste(class(var_r), collapse = ", ")
          )
        }
        # create 1x1 array [draws=1, chains=1]
        # dimension validation will catch if this is inconsistent
        var_r <- array(var_numeric, dim = c(1, 1))
        draws_arrays[[var_name]] <- var_r
      }, error = function(e) {
        stop2(
          "Failed to convert variable '", var_name, "' to array format. ",
          "Original error: ", conditionMessage(e), ". ",
          "Variable type: ", typeof(var_r), ", Class: ", paste(class(var_r), collapse = ", ")
        )
      })
    }
  }
  # validate that all variables have consistent dimensions
  # posterior requires all variables to have the same number of iterations and
  # chains
  if (length(draws_arrays) == 0) {
    stop2("No variables found in draws_arrays.")
  }
  # get dimensions from first variable
  first_var <- draws_arrays[[1]]
  if (!is.array(first_var)) {
    # provide detailed error message for debugging
    var_info <- paste(
      sapply(names(draws_arrays), function(nm) {
        var_obj <- draws_arrays[[nm]]
        paste0(
          nm, ": type=", typeof(var_obj), ", class=",
          paste(class(var_obj), collapse = ","), ", is.array=", is.array(var_obj)
        )
      }),
      collapse = "; "
    )
    stop2(
      "Unexpected format for nutpie draws. ",
      "First variable is not an array. ",
      "Variable details: ", var_info
    )
  }
  first_dims <- dim(first_var)
  n_draws <- first_dims[1]  # after permutation: [draws, chains, ...]
  n_chains <- first_dims[2]
  # validate all variables have the same first two dimensions
  for (var_name in names(draws_arrays)) {
    var_array <- draws_arrays[[var_name]]
    if (!is.array(var_array)) {
      stop2("Variable '", var_name, "' is not an array.")
    }
    var_dims <- dim(var_array)
    if (length(var_dims) < 2) {
      stop2("Variable '", var_name, "' has fewer than 2 dimensions.")
    }
    if (var_dims[1] != n_draws || var_dims[2] != n_chains) {
      stop2(
        "Variable '", var_name, "' has inconsistent dimensions. ",
        "Expected [", n_draws, ", ", n_chains, ", ...], ",
        "but got [", paste(var_dims, collapse = ", "), "]."
      )
    }
  }
  # debug: Print dimensions before creating draws_array
  # this helps identify which variable is causing issues
  if (getOption("coevolve.debug", FALSE)) {
    cat("Variable dimensions before creating draws_array:\n")
    for (var_name in names(draws_arrays)) {
      cat("  ", var_name, ": [", paste(dim(draws_arrays[[var_name]]),
                                       collapse = ", "), "]\n")
    }
  }
  # create draws_array using posterior package
  # posterior::as_draws_array() expects a 3D array
  # [iterations, chains, variables]
  # so we need to combine all 2D arrays [iterations, chains] into
  # a single 3D array
  tryCatch({
    # Get variable names in order
    var_names_ordered <- names(draws_arrays)
    n_vars <- length(var_names_ordered)
    # combine all arrays into a single 3D array [iterations, chains, variables]
    # all arrays should have the same first two dimensions [iterations, chains]
    # set proper iteration and chain indices for ESS calculation
    combined_array <- array(
      dim = c(n_draws, n_chains, n_vars),
      dimnames = list(
        iteration = seq_len(n_draws),
        chain = seq_len(n_chains),
        variable = var_names_ordered
      )
    )
    for (i in seq_along(var_names_ordered)) {
      var_name <- var_names_ordered[i]
      combined_array[, , i] <- draws_arrays[[var_name]]
    }
    # create draws_array from the combined 3D array
    draws_array <- posterior::as_draws_array(combined_array)
    # set variable names properly
    dimnames(draws_array)$variable <- var_names_ordered
  }, error = function(e) {
    # provide more detailed error message
    error_msg <- conditionMessage(e)
    dims_info <- paste(
      sapply(names(draws_arrays), function(nm) {
        paste0(nm, ": [", paste(dim(draws_arrays[[nm]]), collapse = ", "), "]")
      }),
      collapse = "; "
    )
    stop2(
      "Failed to create draws_array: ", error_msg, "\n",
      "Variable dimensions: ", dims_info
    )
  })
  draws_array
}
