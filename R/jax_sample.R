#' Run MCMC via pure JAX model + nutpie Rust sampler
#'
#' @description Builds a pure JAX log-density function, compiles it with
#'   nutpie via \code{from_pyfunc}, and samples using nutpie's Rust NUTS.
#'
#' @param data_list Named R list (JAX-ready, from standata_to_jax).
#' @param num_chains Integer. Number of MCMC chains.
#' @param num_samples Integer. Post-warmup draws per chain.
#' @param num_warmup Integer. Warmup draws per chain.
#' @param seed Integer. Random seed.
#' @param target_accept Numeric. Target acceptance probability for NUTS.
#' @param nutpie_args Named list of extra arguments passed to
#'   \code{nutpie.sample}.
#'
#' @returns A list with elements \code{trace} (ArviZ InferenceData) and
#'   \code{draws_array} (posterior::draws_array).
#'
#' @noRd
jax_run_nutpie <- function(data_list,
                           num_chains = 4L,
                           num_samples = 1000L,
                           num_warmup = 1000L,
                           seed = 0L,
                           target_accept = 0.95,
                           nutpie_args = list()) {
  Sys.setenv(PYTHONIOENCODING = "utf-8")
  Sys.setenv(PYTHONUNBUFFERED = "1")
  Sys.setenv(JAX_PLATFORMS = "cpu")

  py_data <- convert_r_to_python_data_jax(data_list) # nolint
  jax_mod <- load_jax_model_module(convert = FALSE) # nolint

  message("Compiling JAX model with nutpie ...")
  build_result <- jax_mod$build_nutpie_model(py_data)
  compiled <- build_result[[0]]

  message("Sampling with nutpie (Rust NUTS + JAX gradients) ...")
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
  if (length(nutpie_args)) {
    sample_kwargs[names(nutpie_args)] <- nutpie_args
  }

  # Start sampling in background thread so R can show progress
  do.call(jax_mod$start_sampling, sample_kwargs)

  while (TRUE) {
    status <- reticulate::py_to_r(jax_mod$check_sampling())
    if (isTRUE(status$done)) break
    prog <- status$progress_text
    if (nzchar(prog)) {
      cat(sprintf("\r%s", prog))
    } else {
      cat(sprintf(
        "\r  Sampling %d chains x %d draws | %.0fs elapsed",
        num_chains, num_warmup + num_samples, status$elapsed
      ))
    }
    utils::flush.console()
    Sys.sleep(0.2)
  }
  cat("\n")

  if (!is.null(status$error)) {
    stop2("nutpie sampling failed: ", status$error)
  }

  trace <- jax_mod$collect_result()

  draws_array <- convert_jax_draws(trace, num_chains)

  list(trace = trace, draws_array = draws_array)
}

#' Convert JAX/nutpie trace to posterior draws_array
#'
#' @param trace Nutpie trace object (ArviZ InferenceData).
#' @param num_chains Integer. Number of chains used.
#'
#' @returns A draws_array from the posterior package.
#'
#' @noRd
convert_jax_draws <- function(trace, num_chains) {
  draws_list <- list()

  if (!reticulate::py_has_attr(trace, "posterior")) {
    stop2("Expected ArviZ InferenceData with posterior group.")
  }

  posterior_ds <- trace$posterior
  py_builtins <- reticulate::import_builtins()
  var_names <- reticulate::py_to_r(
    py_builtins$list(posterior_ds$data_vars$keys())
  )

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
      n_draws <- dims_perm[1]
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

  if (length(draws_list) == 0) {
    stop2("No variables found in JAX/nutpie trace.")
  }

  first <- draws_list[[1]]
  n_draws <- as.integer(dim(first)[1])
  n_chains <- as.integer(dim(first)[2])
  n_vars <- length(draws_list)
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
