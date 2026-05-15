#' Compare Stan and JAX log density at the same unconstrained point
#'
#' Internal diagnostic used in development and the package test suite to
#' verify that the JAX log-density agrees with Stan's. Compiles both
#' models, evaluates log-density and gradient at the same unconstrained
#' parameter vector, and checks that the log-density offset is constant
#' across points and gradients agree to a tolerance.
#'
#' @importFrom stats rnorm
#' @noRd
compare_stan_jax_logprob <- function(
    data,
    variables,
    id,
    tree,
    effects_mat = NULL,
    complete_cases = FALSE,
    lon_lat = NULL,
    dist_k = NA,
    dist_cov = "exp_quad",
    measurement_error = NULL,
    prior = NULL,
    scale = TRUE,
    estimate_correlated_drift = FALSE,
    estimate_residual = TRUE,
    prior_only = TRUE,
    seed = 1L,
    n_points = 5L,
    grad_tol = 1e-6) {

  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop2("Package 'cmdstanr' is required.")
  }
  stop_if_jax_not_available() # nolint

  distributions <- as.character(variables)

  # --- Stan model ---
  sc <- coev_make_stancode(
    data, variables, id, tree, effects_mat,
    complete_cases, lon_lat, dist_k, dist_cov,
    measurement_error, prior, scale,
    estimate_correlated_drift, estimate_residual,
    log_lik = FALSE, prior_only = prior_only
  )
  sd <- coev_make_standata(
    data, variables, id, tree, effects_mat,
    complete_cases, lon_lat, dist_k, dist_cov,
    measurement_error, prior, scale,
    estimate_correlated_drift, estimate_residual,
    log_lik = FALSE, prior_only = prior_only
  )

  mod <- cmdstanr::cmdstan_model(
    cmdstanr::write_stan_file(sc), force_recompile = TRUE
  )
  fit <- suppressWarnings(mod$sample(
    data = sd, chains = 1L, iter_warmup = 1L,
    iter_sampling = 1L, seed = 1L, refresh = 0,
    show_messages = FALSE, init = 0
  ))
  n_upars <- ncol(posterior::as_draws_matrix(
    fit$unconstrain_draws(draws = fit$draws())
  ))

  # --- JAX model ---
  cfg <- coev_make_model_config( # nolint
    data, variables, id, tree, effects_mat,
    complete_cases, lon_lat, dist_k, dist_cov,
    measurement_error, prior, scale,
    estimate_correlated_drift, estimate_residual,
    prior_only = prior_only
  )
  sd_jax <- embed_model_config( # nolint
    standata_to_jax(sd, distributions), cfg # nolint
  )
  py_data <- convert_r_to_python_data_jax(sd_jax) # nolint

  jax_mod <- load_jax_model_module(convert = FALSE) # nolint
  jax <- reticulate::import("jax", convert = FALSE)
  jax$config$update("jax_platforms", "cpu")
  jax$config$update("jax_enable_x64", TRUE)

  model_obj <- jax_mod$CoevJaxModel()
  model_obj$build(py_data)
  jax_ndim <- reticulate::py_to_r(model_obj$ndim)

  if (n_upars != jax_ndim) {
    stop2(
      "Dimension mismatch: Stan has ", n_upars,
      " unconstrained params, JAX has ", jax_ndim
    )
  }

  np <- reticulate::import("numpy", convert = FALSE)
  jnp <- reticulate::import("jax.numpy", convert = FALSE)

  # Set up JIT-compiled value_and_grad
  py_main <- reticulate::import("__main__", convert = FALSE)
  py_main$model_obj <- model_obj
  reticulate::py_run_string(
    "import jax; vg = jax.jit(jax.value_and_grad(model_obj.log_density))"
  )

  # Evaluate at random points
  set.seed(as.integer(seed))
  offsets <- numeric(n_points)
  max_grad_diffs <- numeric(n_points)

  for (i in seq_len(n_points)) {
    u <- rnorm(n_upars, 0, 0.3)

    logp_stan <- as.numeric(
      fit$log_prob(u, jacobian = TRUE)
    )
    stan_grad <- fit$grad_log_prob(u, jacobian = TRUE)

    u_jax <- jnp$array(np$array(u, dtype = "float64"))
    py_main$u_jax <- u_jax
    reticulate::py_run_string("
import numpy as np
lp, g = vg(u_jax)
lp.block_until_ready()
_lp_val = float(lp.item())
_grad_np = np.array(g)
")
    logp_jax <- reticulate::py_to_r(py_main$`_lp_val`)
    jax_grad <- reticulate::py_to_r(py_main$`_grad_np`)

    offsets[i] <- logp_stan - logp_jax
    max_grad_diffs[i] <- max(abs(stan_grad - jax_grad))
  }

  offset_sd <- if (n_points > 1) sd(offsets) else 0
  max_grad_diff <- max(max_grad_diffs)

  if (offset_sd > 1e-6) {
    warning2(
      "logp offset is NOT constant (sd = ",
      format(offset_sd, digits = 4),
      "). Stan and JAX log-densities disagree."
    )
  }
  if (max_grad_diff > grad_tol) {
    warning2(
      "Max gradient discrepancy = ",
      format(max_grad_diff, digits = 4),
      " exceeds tolerance ", grad_tol
    )
  }

  message(
    "Stan vs JAX: ndim=", n_upars,
    ", constant_offset=", format(offsets[1], digits = 6),
    ", offset_sd=", format(offset_sd, digits = 4),
    ", max_grad_diff=", format(max_grad_diff, digits = 4)
  )

  list(
    n_params         = n_upars,
    constant_offset  = offsets[1],
    offset_sd        = offset_sd,
    max_grad_diff    = max_grad_diff,
    mean_grad_diff   = mean(max_grad_diffs)
  )
}
