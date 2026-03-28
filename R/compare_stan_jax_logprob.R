#' Compare Stan and JAX log density at the same parameters
#'
#' Compiles both Stan and JAX models, evaluates their log densities at
#' the same parameter point, and returns the difference. No MCMC sampling
#' is used: the evaluation point is a seeded random draw in Stan's
#' unconstrained space.
#'
#' @inheritParams coev_make_stancode
#' @param seed Integer seed for the evaluation point (default \code{1L}).
#' @param tol_warn Absolute difference above which a warning is issued
#'   (default \code{0.1} on log scale).
#'
#' @returns A list with \code{logprob_stan}, \code{logprob_jax},
#'   \code{abs_diff}, and \code{n_params}.
#'
#' @export
compare_stan_jax_logprob <- function(
    data,
    variables,
    id,
    tree,
    effects_mat = NULL,
    complete_cases = FALSE,
    dist_mat = NULL,
    dist_cov = "exp_quad",
    measurement_error = NULL,
    prior = NULL,
    scale = TRUE,
    estimate_correlated_drift = FALSE,
    estimate_residual = TRUE,
    prior_only = TRUE,
    seed = 1L,
    tol_warn = 0.1) {

  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop2(
      "Package 'cmdstanr' is required for ",
      "compare_stan_jax_logprob()."
    )
  }
  stop_if_jax_not_available()

  distributions <- as.character(variables)

  sc <- coev_make_stancode(
    data, variables, id, tree, effects_mat,
    complete_cases, dist_mat, dist_cov,
    measurement_error, prior, scale,
    estimate_correlated_drift, estimate_residual,
    log_lik = FALSE, prior_only = prior_only
  )
  cfg <- coev_make_model_config(
    data, variables, id, tree, effects_mat,
    complete_cases, dist_mat, dist_cov,
    measurement_error, prior, scale,
    estimate_correlated_drift, estimate_residual,
    prior_only = prior_only
  )
  sd <- coev_make_standata(
    data, variables, id, tree, effects_mat,
    complete_cases, dist_mat, dist_cov,
    measurement_error, prior, scale,
    estimate_correlated_drift, estimate_residual,
    log_lik = FALSE, prior_only = prior_only
  )
  sd_jax <- embed_model_config(
    standata_to_jax(sd, distributions), cfg
  )

  # Compile Stan and get unconstrained dimension
  mod <- cmdstanr::cmdstan_model(
    cmdstanr::write_stan_file(sc),
    force_recompile = TRUE
  )
  fit <- mod$sample(
    data = sd,
    chains = 1L,
    parallel_chains = 1L,
    iter_warmup = 1L,
    iter_sampling = 1L,
    seed = 1L,
    refresh = 0,
    show_messages = FALSE,
    init = 0
  )

  n_upars <- ncol(
    posterior::as_draws_matrix(
      fit$unconstrain_draws(draws = fit$draws())
    )
  )

  # Evaluate Stan at two points:
  # reference (u=0) and target (seeded random)
  u_0 <- numeric(n_upars)
  set.seed(as.integer(seed))
  u_1 <- rnorm(n_upars, 0, 0.5)

  logp_stan_0 <- as.numeric(
    fit$log_prob(u_0, jacobian = FALSE)
  )
  logp_stan_1 <- as.numeric(
    fit$log_prob(u_1, jacobian = FALSE)
  )

  # Build JAX model and evaluate at the same points
  py_data <- convert_r_to_python_data_jax(sd_jax)
  jax_mod <- load_jax_model_module(convert = FALSE)
  jax <- reticulate::import("jax", convert = FALSE)
  jax$config$update("jax_platforms", "cpu")
  jax$config$update("jax_enable_x64", TRUE)

  model_obj <- jax_mod$CoevJaxModel()
  model_obj$build(py_data)

  # JAX model uses same unconstrained param vector
  # but has different dimension (excludes Stan's
  # tip-edge z_drift, uses different LKJ param)
  # So we evaluate both at their own unconstrained=0
  # reference, then compare the DIFFERENCE
  jax_ndim <- reticulate::py_to_r(model_obj$ndim)

  np <- reticulate::import("numpy", convert = FALSE)
  jnp <- reticulate::import("jax.numpy", convert = FALSE)

  jax_u_0 <- jnp$zeros(as.integer(jax_ndim), dtype = "float64")
  logp_jax_0 <- reticulate::py_to_r(
    model_obj$log_density(jax_u_0)
  )

  set.seed(as.integer(seed))
  jax_u_1_r <- rnorm(jax_ndim, 0, 0.5)
  jax_u_1 <- jnp$array(
    np$array(jax_u_1_r, dtype = "float64")
  )
  logp_jax_1 <- reticulate::py_to_r(
    model_obj$log_density(jax_u_1)
  )

  # Compare log-density DIFFERENCES to cancel
  # normalization constants
  stan_diff <- logp_stan_1 - logp_stan_0
  jax_diff <- logp_jax_1 - logp_jax_0

  abs_diff <- abs(stan_diff - jax_diff)

  if (abs_diff > tol_warn) {
    warning(
      "Large |delta_logp Stan - delta_logp JAX| = ",
      format(abs_diff, digits = 5),
      ". Stan n_upars=", n_upars,
      ", JAX ndim=", jax_ndim,
      call. = FALSE
    )
  }

  list(
    logprob_stan_diff = stan_diff,
    logprob_jax_diff  = jax_diff,
    abs_diff          = abs_diff,
    n_params_stan     = n_upars,
    n_params_jax      = jax_ndim
  )
}
