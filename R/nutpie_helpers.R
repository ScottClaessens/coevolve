#' Check if nutpie is available via reticulate
#'
#' @srrstats {G1.4} Function is documented
#'
#' @description Checks if nutpie can be imported via reticulate. This function
#'   attempts to import nutpie using the currently configured Python environment
#'   (via reticulate). Users must configure reticulate to use a Python environment
#'   with nutpie installed before using this function.
#'
#' @details To configure nutpie, you must:
#' \enumerate{
#'   \item Install nutpie in a Python environment: \code{pip install "nutpie[stan]"}
#'   \item Configure reticulate to use that Python environment:
#'     \itemize{
#'       \item Set \code{RETICULATE_PYTHON} environment variable to the Python executable path
#'       \item Or use \code{reticulate::use_python()} or \code{reticulate::use_virtualenv()}
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
  # Check if reticulate is available
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    return(FALSE)
  }
  
  # Try to import nutpie using the currently configured Python
  # No automatic searching - user must configure reticulate explicitly
  tryCatch({
    nutpie <- reticulate::import("nutpie", convert = FALSE)
    # If we get here, nutpie is available
    return(TRUE)
  }, error = function(e) {
    # nutpie not available
    return(FALSE)
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
#' compiled <- nutpie_compile_stan_model(stan_code)
#' }
#'
#' @export
nutpie_compile_stan_model <- function(stan_code) {
  # Check if nutpie is available
  if (!check_nutpie_available()) {
    stop2(
      "nutpie is not available. Please install nutpie with: ",
      "pip install 'nutpie[stan]'"
    )
  }
  
  # Import nutpie
  nutpie <- reticulate::import("nutpie", convert = FALSE)
  
  # Compile Stan model
  tryCatch({
    compiled <- nutpie$compile_stan_model(code = stan_code)
    return(compiled)
  }, error = function(e) {
    # Convert Python error to R error
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
  # Import numpy for array conversion
  np <- reticulate::import("numpy", convert = FALSE)
  
  # Convert each element of the data list
  python_data <- list()
  
  for (name in names(data_list)) {
    value <- data_list[[name]]
    
    # Handle different R types
    if (is.matrix(value) || is.array(value)) {
      # Check if matrix/array should be integer type
      # Stan integer variables: tip, node_seq, parent, effects_mat, tip_id, etc.
      integer_vars <- c("tip", "node_seq", "parent", "effects_mat", "tip_id", 
                        "N_tips", "N_tree", "N_obs", "J", "N_seg", "num_effects",
                        "prior_only", "miss")
      
      # Check if this is a known integer variable or if values are integers
      is_known_integer <- name %in% integer_vars
      is_integer_type <- is.integer(value)
      
      # Check if all non-NA values are integers (within tolerance)
      value_no_na <- value[!is.na(value)]
      all_integers <- if (length(value_no_na) > 0) {
        all(abs(value_no_na - round(value_no_na)) < 1e-10)
      } else {
        TRUE  # Empty or all NA - treat as integer if known integer var
      }
      
      # Convert to integer if it's a known integer variable or all values are integers
      if (is_known_integer || (is_integer_type && all_integers)) {
        # Convert to integer array, preserving matrix/array structure
        # Don't use as.integer() directly as it may flatten matrices
        # Instead, convert element-wise while preserving dimensions
        value_int <- array(as.integer(value), dim = dim(value))
        # NAs become NA_integer_ which numpy will handle
        # For Stan arrays like array[N_tree, N_seg], ensure 2D structure is preserved
        # even when N_tree=1 (numpy might squeeze single dimensions)
        if (name %in% c("node_seq", "parent", "ts", "tip") && length(dim(value_int)) == 2) {
          # Use ndmin=2 to prevent numpy from squeezing dimensions
          python_data[[name]] <- np$array(value_int, dtype = "int32", ndmin = 2L)
        } else {
          python_data[[name]] <- np$array(value_int, dtype = "int32")
        }
      } else {
        # Convert to float array, preserving structure
        # For multi-dim arrays like ts, ensure 2D structure is preserved
        if (name %in% c("ts") && length(dim(value)) == 2) {
          python_data[[name]] <- np$array(array(as.double(value), dim = dim(value)), 
                                          dtype = "float64", ndmin = 2L)
        } else {
          python_data[[name]] <- np$array(array(as.double(value), dim = dim(value)), 
                                          dtype = "float64")
        }
      }
    } else if (is.integer(value) || (is.numeric(value) && 
                                      all(abs(value - round(value)) < 1e-10, na.rm = TRUE))) {
      # Convert to Python integer or integer array
      if (length(value) == 1) {
        python_data[[name]] <- as.integer(value)
      } else {
        python_data[[name]] <- np$array(as.integer(value), dtype = "int32")
      }
    } else if (is.numeric(value)) {
      # Convert to Python float or float array
      if (length(value) == 1) {
        python_data[[name]] <- as.double(value)
      } else {
        python_data[[name]] <- np$array(as.double(value), dtype = "float64")
      }
    } else {
      # Try to convert as-is (might be list/other)
      python_data[[name]] <- value
    }
  }
  
  # Convert to Python dict
  return(reticulate::r_to_py(python_data))
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
#'   nutpie.
#' @param ... Additional arguments passed to nutpie.sample().
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
  # Check if nutpie is available
  if (!check_nutpie_available()) {
    stop2(
      "nutpie is not available. Please install nutpie with: ",
      "pip install 'nutpie[stan]'"
    )
  }
  
  # Import nutpie
  nutpie <- reticulate::import("nutpie", convert = FALSE)
  
  # Compile model if stan_code is a character string
  if (is.character(stan_code)) {
    compiled <- nutpie_compile_stan_model(stan_code)
  } else {
    # Assume it's already compiled
    compiled <- stan_code
  }
  
  # Convert R data to Python format
  python_data <- convert_r_to_python_data(data_list)
  
  # Prepare sampling arguments with correct nutpie parameter names
  # nutpie uses: draws, tune, chains (not num_samples, num_warmup, num_chains)
  # Based on docstring: "draws: The number of draws after tuning in each chain"
  # and "tune: The number of tuning (warmup) draws in each chain"
  # So both should be per-chain values, matching cmdstanr's behavior.
  sample_args <- list(
    chains = as.integer(num_chains),
    draws = as.integer(num_samples),
    tune = as.integer(num_warmup)
  )
  
  # Add seed if provided
  if (!is.null(seed)) {
    sample_args$seed <- as.integer(seed)
  }
  
  # Add target_accept if provided (maps from adapt_delta in Stan)
  if (!is.null(target_accept)) {
    sample_args$target_accept <- as.double(target_accept)
  }
  
  # Add low_rank_modified_mass_matrix if provided
  if (!is.null(low_rank_modified_mass_matrix)) {
    sample_args$low_rank_modified_mass_matrix <- as.logical(low_rank_modified_mass_matrix)
  }
  
  # nutpie returns InferenceData by default (not raw trace)
  # We'll extract draws from trace$posterior in convert_nutpie_draws()
  # Only set return_raw_trace if user explicitly requests it
  additional_args <- list(...)
  
  # Ensure tune, draws, and chains are not overridden by additional_args
  # These are core parameters that should come from the function arguments
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
    # Default: use InferenceData (return_raw_trace = FALSE)
    sample_args <- c(sample_args, additional_args)
  } else {
    # User explicitly set return_raw_trace, use their value
    sample_args$return_raw_trace <- additional_args$return_raw_trace
    additional_args$return_raw_trace <- NULL
    sample_args <- c(sample_args, additional_args)
  }
  
  # Prepare data using with_data method
  tryCatch({
    # nutpie uses with_data() method to attach data
    # In Python: compiled.with_data(**data_dict)
    # In R via reticulate, we need to unpack the dict as keyword arguments
    # We'll define a Python helper function to handle **kwargs
    
    # Define Python helper function to handle **kwargs unpacking
    # This is needed because reticulate doesn't directly support **kwargs syntax
    # We define it each time (it's lightweight and ensures it's available)
    reticulate::py_run_string("
def nutpie_call_with_data(compiled_model, data_dict):
    return compiled_model.with_data(**data_dict)
")
    
    # Call the helper function
    call_func <- reticulate::py$nutpie_call_with_data
    model_with_data <- call_func(compiled, python_data)
    
    # Sample
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
#' # After sampling with nutpie
#' trace <- nutpie_sample(stan_code, data_list)
#' draws <- convert_nutpie_draws(trace)
#' }
#'
#' @export
convert_nutpie_draws <- function(trace) {
  # Check if posterior package is available
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop2("The 'posterior' package is required for draw conversion.")
  }
  
  # Extract draws from nutpie InferenceData object
  # nutpie returns InferenceData (ArviZ) with posterior as xarray Dataset
  # Format: trace$posterior$data_vars contains variables with dims (chain, draw, ...)
  tryCatch({
    # Check if this is InferenceData (has posterior attribute)
    if (reticulate::py_has_attr(trace, "posterior")) {
      posterior_ds <- trace$posterior
      # Get variable names from data_vars
      # keys() returns a Python dict_keys object, convert to list using Python's list()
      keys_obj <- posterior_ds$data_vars$keys()
      py_builtins <- reticulate::import_builtins()
      var_names <- reticulate::py_to_r(py_builtins$list(keys_obj))
    } else {
      # Try raw trace with draws() method (if return_raw_trace=TRUE was used)
      draws_dict <- trace$draws()
      var_names <- names(reticulate::py_to_r(draws_dict))
    }
  }, error = function(e) {
    stop2("Failed to extract draws from nutpie trace: ", conditionMessage(e))
  })
  
  if (length(var_names) == 0) {
    stop2("No variables found in nutpie trace.")
  }
  
  # Process each variable
  # cmdstanr flattens all parameter dimensions into separate scalar variables
  # e.g., z[1,2,3] becomes separate variables z[1,2,3]
  # So we need to flatten arrays and create separate variables for each element
  draws_arrays <- list()
  
  for (var_name in var_names) {
    tryCatch({
      if (reticulate::py_has_attr(trace, "posterior")) {
        # Extract from InferenceData posterior Dataset
        var_xarray <- trace$posterior[[var_name]]
        # Convert xarray DataArray to numpy array
        # xarray format: (chain, draw, ...) -> numpy array
        var_array <- var_xarray$values
        var_r <- reticulate::py_to_r(var_array)
      } else {
        # Extract from raw trace draws dict
        draws_dict <- trace$draws()
        draws_list <- reticulate::py_to_r(draws_dict)
        var_array <- draws_list[[var_name]]
        var_r <- reticulate::py_to_r(var_array)
      }
    }, error = function(e) {
      stop2("Failed to convert variable '", var_name, "' to R format: ", conditionMessage(e))
    })
    
    # nutpie/ArviZ format: [chains, draws, dims...]
    # posterior draws_array format: [iterations, chains, dims...]
    # So we need to swap first two dimensions: [chains, draws, ...] -> [draws, chains, ...]
    if (is.array(var_r)) {
      dims_var <- dim(var_r)
      n_dims <- length(dims_var)
      
      if (n_dims == 2) {
        # Scalar variable: [chains, draws] -> [draws, chains]
        var_r <- aperm(var_r, c(2, 1))
        # Add as single variable
        draws_arrays[[var_name]] <- var_r
      } else if (n_dims > 2) {
        # Array variable: [chains, draws, dims...] -> [draws, chains, dims...]
        perm <- c(2, 1, 3:n_dims)
        var_r <- aperm(var_r, perm)
        
        # Get dimensions after permutation
        dims_permuted <- dim(var_r)
        n_draws_perm <- dims_permuted[1]
        n_chains_perm <- dims_permuted[2]
        
        # Flatten parameter dimensions into separate variables
        # Reshape to [iterations, chains, param_elements]
        param_dims <- dims_permuted[3:n_dims]
        n_param_elements <- prod(param_dims)
        var_reshaped <- array(var_r, dim = c(n_draws_perm, n_chains_perm, n_param_elements))
        
        # Create separate variable for each parameter element
        for (i in seq_len(n_param_elements)) {
          # Get indices for this element
          indices <- arrayInd(i, param_dims)
          # Create variable name with indices
          if (length(indices) == 1) {
            var_name_indexed <- paste0(var_name, "[", indices[1], "]")
          } else {
            var_name_indexed <- paste0(var_name, "[", paste(indices, collapse = ","), "]")
          }
          # Extract this element: [iterations, chains, element] -> [iterations, chains]
          draws_arrays[[var_name_indexed]] <- var_reshaped[, , i]
        }
      }
    } else {
      # Not an array - add as-is
      draws_arrays[[var_name]] <- var_r
    }
  }
  
  # Validate that all variables have consistent dimensions
  # posterior requires all variables to have the same number of iterations and chains
  if (length(draws_arrays) == 0) {
    stop2("No variables found in draws_arrays.")
  }
  
  # Get dimensions from first variable
  first_var <- draws_arrays[[1]]
  if (!is.array(first_var)) {
    stop2("Unexpected format for nutpie draws.")
  }
  
  first_dims <- dim(first_var)
  n_draws <- first_dims[1]  # After permutation: [draws, chains, ...]
  n_chains <- first_dims[2]
  
  # Validate all variables have the same first two dimensions
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
  
  # Debug: Print dimensions before creating draws_array
  # This helps identify which variable is causing issues
  if (getOption("coevolve.debug", FALSE)) {
    cat("Variable dimensions before creating draws_array:\n")
    for (var_name in names(draws_arrays)) {
      cat("  ", var_name, ": [", paste(dim(draws_arrays[[var_name]]), collapse = ", "), "]\n")
    }
  }
  
  # Create draws_array using posterior package
  # posterior::as_draws_array() expects a 3D array [iterations, chains, variables]
  # So we need to combine all 2D arrays [iterations, chains] into a single 3D array
  tryCatch({
    # Get variable names in order
    var_names_ordered <- names(draws_arrays)
    n_vars <- length(var_names_ordered)
    
    # Combine all arrays into a single 3D array [iterations, chains, variables]
    # All arrays should have the same first two dimensions [iterations, chains]
    combined_array <- array(
      dim = c(n_draws, n_chains, n_vars),
      dimnames = list(
        iteration = NULL,
        chain = NULL,
        variable = var_names_ordered
      )
    )
    
    for (i in seq_along(var_names_ordered)) {
      var_name <- var_names_ordered[i]
      combined_array[, , i] <- draws_arrays[[var_name]]
    }
    
    # Create draws_array from the combined 3D array
    draws_array <- posterior::as_draws_array(combined_array)
    
    # Set variable names properly
    dimnames(draws_array)$variable <- var_names_ordered
    
  }, error = function(e) {
    # Provide more detailed error message
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
  
  return(draws_array)
}

