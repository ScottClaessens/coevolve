#' Run PyMC MCMC and return draws in posterior format
#'
#' @description Executes a PyMC model via reticulate. The generated Python
#'   model code is sourced into the Python environment, data is converted, and
#'   NUTS sampling is run. Samples are returned as a posterior draws_array.
#'
#' @param pymc_code Character string of generated Python model code.
#' @param data_list Named R list (PyMC-ready, from standata_to_pymc).
#' @param num_chains Integer. Number of MCMC chains.
#' @param num_samples Integer. Post-warmup draws per chain.
#' @param num_warmup Integer. Warmup draws per chain.
#' @param seed Integer. Random seed.
#' @param target_accept Numeric. Target acceptance probability for NUTS.
#' @param compile_mode Character. PyTensor compilation mode: \code{"cpu"}
#'   (default) or \code{"mlx"} for Apple Silicon GPU via MLX.
#'
#' @returns A list with elements \code{trace} (ArviZ InferenceData) and
#'   \code{draws_array} (posterior::draws_array).
#'
#' @noRd
pymc_run_mcmc <- function(pymc_code, data_list,
                          num_chains = 4L, num_samples = 1000L,
                          num_warmup = 1000L, seed = 0L,
                          target_accept = 0.95,
                          compile_mode = "cpu") {
  compile_mode <- tolower(compile_mode)
  use_mlx <- identical(compile_mode, "mlx")

  if (use_mlx) {
    mlx_ok <- mlx_backend_viable()
    if (!mlx_ok) {
      message(
        "MLX backend requested but not available. ",
        "Falling back to default (C) backend. ",
        "Install MLX with: pip install mlx"
      )
      use_mlx <- FALSE
    }
  }

  Sys.setenv(PYTHONIOENCODING = "utf-8")

  pm <- reticulate::import("pymc", convert = FALSE)
  np <- reticulate::import("numpy", convert = FALSE)

  reticulate::py_run_string(pymc_code)
  build_fn <- reticulate::py$build_model

  py_data <- convert_r_to_python_data_pymc(data_list)
  model <- build_fn(py_data, compile_mode)

  sample_kwargs <- list(
    model         = model,
    draws         = as.integer(num_samples),
    tune          = as.integer(num_warmup),
    chains        = as.integer(num_chains),
    target_accept = target_accept,
    random_seed   = as.integer(seed),
    cores         = 1L,
    init          = "adapt_full",
    return_inferencedata = FALSE,
    progressbar   = FALSE,
    compute_convergence_checks = FALSE
  )

  if (use_mlx) {
    message("Running PyMC with MLX backend (Apple Silicon GPU).")
    sample_kwargs$compile_kwargs <- list(mode = "MLX")
  } else {
    message("Running PyMC with default (C) backend.")
  }

  trace <- do.call(pm$sample, sample_kwargs)

  draws_array <- convert_pymc_draws(trace, num_chains)

  list(trace = trace, draws_array = draws_array)
}

#' Run PyMC MCMC via nutpie (Rust NUTS sampler)
#'
#' @description Builds a PyMC model, compiles it with nutpie, and samples
#'   using nutpie's Rust-based NUTS implementation. Returns draws in the same
#'   format as \code{pymc_run_mcmc}.
#'
#' @inheritParams pymc_run_mcmc
#' @param low_rank_modified_mass_matrix Logical. If TRUE, use nutpie's
#'   low-rank modified mass matrix adaptation (default: FALSE).
#'
#' @returns A list with elements \code{trace} (ArviZ InferenceData) and
#'   \code{draws_array} (posterior::draws_array).
#'
#' @noRd
pymc_run_nutpie <- function(pymc_code, data_list,
                            num_chains = 4L, num_samples = 1000L,
                            num_warmup = 1000L, seed = 0L,
                            target_accept = 0.95,
                            nutpie_backend = "jax",
                            nutpie_args = list()) {

  Sys.setenv(PYTHONIOENCODING = "utf-8")

  # Force JAX to use CPU — GPU adds overhead for the small matrix ops typical
  # of phylogenetic coevolutionary models (J ~ 2-10, N_tips ~ 50-500).
  # Benchmarks show nutpie JAX on CPU is 2-4x faster than GPU for these sizes.
  if (identical(tolower(nutpie_backend), "jax")) {
    Sys.setenv(JAX_PLATFORMS = "cpu")
  }

  nutpie <- reticulate::import("nutpie", convert = FALSE)
  np     <- reticulate::import("numpy",  convert = FALSE)

  reticulate::py_run_string(pymc_code)
  build_fn <- reticulate::py$build_model

  py_data <- convert_r_to_python_data_pymc(data_list)
  model   <- build_fn(py_data, "cpu")

  # JAX backend compiles scan graphs much faster than numba for large models
  use_backend <- nutpie_backend
  if (!("gradient_backend" %in% names(nutpie_args))) {
    nutpie_args$gradient_backend <- use_backend
  }
  compile_arg_names <- c(
    "gradient_backend", "initial_points", "jitter_rvs",
    "default_initialization_strategy", "var_names", "freeze_model"
  )
  compile_args <- nutpie_args[intersect(names(nutpie_args), compile_arg_names)]
  sample_args_extra <- nutpie_args[setdiff(names(nutpie_args), compile_arg_names)]
  if (use_backend == "jax") {
    # nutpie's JAX compiler rejects PyMC models with custom initvals.
    # Clear them here so base PyMC can keep its calibrated defaults.
    reticulate::py_run_string("
def _clear_nutpie_initvals(model):
    for _rv in model.free_RVs:
        model.rvs_to_initial_values[_rv] = None
    return model
")
    model <- reticulate::py$`_clear_nutpie_initvals`(model)
  }
  message(sprintf("Compiling PyMC model with nutpie (%s backend) ...", use_backend))
  compiled <- tryCatch({
    do.call(
      nutpie$compile_pymc_model,
      c(list(model = model, backend = use_backend), compile_args)
    )
  }, error = function(e) {
    if (use_backend == "jax") {
      message("JAX backend failed, falling back to numba ...")
      do.call(
        nutpie$compile_pymc_model,
        c(list(model = model, backend = "numba"), compile_args)
      )
    } else {
      stop(e)
    }
  })

  # nutpie's Rust layer can't handle StringArray coordinates; replace with ints
  reticulate::py_run_string("
import numpy as _np
def _fix_nutpie_coords(compiled):
    for k in list(compiled._coords):
        v = compiled._coords[k]
        if hasattr(v, 'dtype') and v.dtype.kind in ('U', 'S', 'O'):
            compiled._coords[k] = _np.arange(len(v))
    return compiled
")
  compiled <- reticulate::py$`_fix_nutpie_coords`(compiled)

  message("Sampling with nutpie (Rust NUTS) ...")
  sample_kwargs <- list(
    compiled_model = compiled,
    draws          = as.integer(num_samples),
    tune           = as.integer(num_warmup),
    chains         = as.integer(num_chains),
    seed           = as.integer(seed)
  )
  if (!is.null(target_accept)) {
    sample_kwargs$target_accept <- as.double(target_accept)
  }
  if (length(sample_args_extra)) {
    sample_kwargs[names(sample_args_extra)] <- sample_args_extra
  }

  trace <- do.call(nutpie$sample, sample_kwargs)

  draws_array <- convert_pymc_draws(trace, num_chains)

  list(trace = trace, draws_array = draws_array)
}


#' Probe MLX backend viability
#'
#' Checks whether MLX is importable and PyTensor's MLX mode is available.
#'
#' @returns Logical.
#' @noRd
mlx_backend_viable <- function() {
  tryCatch({
    reticulate::py_run_string("
import mlx.core as mx
import pymc as pm
import pytensor.tensor as pt
from pytensor.compile.mode import get_mode
mlx_mode = get_mode('MLX')

# Test that pm.sample can compile with MLX by building a trivial model
with pm.Model() as _test_m:
    _x = pm.Normal('_x', 0, 1)
pm.sample(model=_test_m, draws=2, tune=2, chains=1,
          random_seed=0, progressbar=False,
          compile_kwargs=dict(mode='MLX'))
", convert = FALSE)
    TRUE
  }, error = function(e) {
    FALSE
  })
}

#' Convert PyMC trace to posterior draws_array
#'
#' @param trace PyMC trace object (MultiTrace or ArviZ InferenceData).
#' @param num_chains Integer. Number of chains used.
#'
#' @returns A draws_array from the posterior package.
#'
#' @noRd
convert_pymc_draws <- function(trace, num_chains) {
  draws_list <- list()
  if (reticulate::py_has_attr(trace, "posterior")) {
    posterior_ds <- trace$posterior
    py_builtins <- reticulate::import_builtins()
    var_names <- reticulate::py_to_r(
      py_builtins$list(posterior_ds$data_vars$keys())
    )
    var_names <- var_names[!grepl("^_", var_names)]

    for (vname in var_names) {
      var_xarray <- posterior_ds[[vname]]
      arr <- reticulate::py_to_r(var_xarray$values)
      if (!is.array(arr) && !is.matrix(arr)) {
        arr <- array(as.numeric(arr), dim = length(arr))
      }
      dims <- dim(arr)
      ndim <- length(dims)

      # ArviZ format: (chain, draw, ...) -> need (draw, chain, ...)
      if (ndim == 2) {
        arr <- aperm(arr, c(2, 1))
        draws_list[[vname]] <- arr
      } else if (ndim > 2) {
        arr <- aperm(arr, c(2, 1, 3:ndim))
        dims_perm <- dim(arr)
        n_draws  <- dims_perm[1]
        n_chains <- dims_perm[2]
        param_dims <- dims_perm[3:ndim]
        n_elem <- prod(param_dims)
        arr_flat <- array(arr, dim = c(n_draws, n_chains, n_elem))

        for (k in seq_len(n_elem)) {
          idx <- arrayInd(k, param_dims)
          idx_str <- paste(idx, collapse = ",")
          slice <- arr_flat[, , k, drop = FALSE]
          dim(slice) <- c(n_draws, n_chains)
          draws_list[[paste0(vname, "[", idx_str, "]")]] <- slice
        }
      }
    }
  } else {
    var_names <- reticulate::py_to_r(trace$varnames)
    var_names <- var_names[!grepl("^_", var_names)]

    for (vname in var_names) {
      arr_list <- reticulate::py_to_r(
        trace$get_values(vname, combine = FALSE, squeeze = FALSE)
      )
      arr <- simplify2array(arr_list)
      dims <- dim(arr)
      ndim <- length(dims)

      # MultiTrace list simplifies to (draw, chain, ...) for arrays
      if (ndim == 2) {
        draws_list[[vname]] <- arr
      } else if (ndim > 2) {
        n_draws  <- dims[1]
        n_chains <- dims[2]
        param_dims <- dims[3:ndim]
        n_elem <- prod(param_dims)
        arr_flat <- array(arr, dim = c(n_draws, n_chains, n_elem))

        for (k in seq_len(n_elem)) {
          idx <- arrayInd(k, param_dims)
          idx_str <- paste(idx, collapse = ",")
          slice <- arr_flat[, , k, drop = FALSE]
          dim(slice) <- c(n_draws, n_chains)
          draws_list[[paste0(vname, "[", idx_str, "]")]] <- slice
        }
      }
    }
  }

  if (length(draws_list) == 0) {
    stop2("No variables found in PyMC trace.")
  }

  first <- draws_list[[1]]
  n_draws  <- as.integer(dim(first)[1])
  n_chains <- as.integer(dim(first)[2])
  n_vars   <- length(draws_list)
  var_ordered <- names(draws_list)

  combined <- array(
    dim = c(n_draws, n_chains, n_vars),
    dimnames = list(
      iteration = seq_len(n_draws),
      chain = seq_len(n_chains),
      variable = var_ordered
    )
  )
  for (i in seq_along(var_ordered)) {
    combined[, , i] <- draws_list[[var_ordered[i]]]
  }

  posterior::as_draws_array(combined)
}
