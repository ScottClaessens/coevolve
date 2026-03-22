#' Compare Stan and PyMC log density at the same parameters
#'
#' Compiles both Stan and PyMC implementations and evaluates their log
#' densities at the same fixed parameter point.  No MCMC sampling is used
#' to choose the evaluation point: the target is a reproducible normal draw
#' in unconstrained space (seeded by \code{seed}); the reference is
#' unconstrained = 0.  A minimal one-step Stan run is used only to
#' initialise \code{log_prob} and \code{constrain_variables} methods.
#'
#' **Normalization / Jacobian:** Both Stan and PyMC are evaluated with
#' \code{jacobian = FALSE}.  The normalization-constant offset (PyMC
#' normalized vs Stan unnormalized) is cancelled via
#' \code{norm_const = logp_pymc_0 - logp_stan_0} computed at the same
#' constrained reference point (Stan's constrained values at unconstrained = 0).
#' By evaluating the reference in both backends at the same constrained
#' values, likelihood terms cancel exactly and only per-distribution
#' normalization constants remain.
#'
#' **Scope:** Stan stores extra \code{z_drift} components on tip edges; those
#' UNNORMALIZED kernels are subtracted into \code{stan_tip_z_log_prior}.  When
#' \code{estimate_correlated_drift = TRUE}, Stan's \code{L_R} (Cholesky
#' factor) is inverted to PyMC's stick-breaking parameters \code{_L_R_raw}.
#'
#' Intended use: one-off debugging and regression tests when Stan and PyMC
#' timings or posteriors disagree.
#'
#' @inheritParams coev_make_stancode
#' @param seed Integer seed for the target unconstrained evaluation point
#'   (default \code{1L}).  Different seeds exercise different parameter
#'   values without requiring MCMC.
#' @param compile_mode PyTensor mode passed to \code{build_model()}
#'   (\code{"cpu"} or \code{"mlx"}).
#' @param tol_warn Absolute difference above which a warning is issued
#'   (default \code{0.05} on log scale).
#'
#' @returns A list with \code{logprob_stan}, \code{logprob_stan_adj},
#'   \code{stan_tip_z_log_prior}, \code{stan_terminal_drift_log_prior},
#'   \code{logprob_pymc}, \code{abs_diff}, \code{n_params}, and
#'   \code{stan_params}.
#'
#' @export
compare_stan_pymc_logprob <- function(
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
    compile_mode = "cpu",
    tol_warn = 0.05) {

  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop2("Package 'cmdstanr' is required for compare_stan_pymc_logprob().")
  }
  stop_if_pymc_not_available()
  if (!is.null(dist_mat)) {
    stop2(
      "GP models (dist_mat) are not yet supported in the PyMC backend: ",
      "the spatial covariance translation is currently being refactored. ",
      "Use nuts_sampler = \"stan\" for GP models."
    )
  }

  distributions <- as.character(variables)

  sc <- coev_make_stancode(
    data, variables, id, tree, effects_mat, complete_cases, dist_mat,
    dist_cov, measurement_error, prior, scale,
    estimate_correlated_drift, estimate_residual,
    log_lik = FALSE, prior_only = prior_only
  )
  cfg <- coev_make_pymc(
    data, variables, id, tree, effects_mat, complete_cases, dist_mat,
    dist_cov, measurement_error, prior, scale,
    estimate_correlated_drift, estimate_residual, prior_only = prior_only
  )
  sd <- coev_make_standata(
    data, variables, id, tree, effects_mat, complete_cases, dist_mat,
    dist_cov, measurement_error, prior, scale,
    estimate_correlated_drift, estimate_residual,
    log_lik = FALSE, prior_only = prior_only
  )
  sd_pymc <- embed_pymc_config(standata_to_pymc(sd, distributions), cfg)

  mod <- cmdstanr::cmdstan_model(
    cmdstanr::write_stan_file(sc),
    force_recompile = TRUE
  )

  # Minimal fit (1 warmup + 1 sample at init = 0) solely to initialise
  # model methods (log_prob, constrain_variables).  The draw is not used.
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

  # Unconstrained parameter dimension.
  n_upars <- ncol(
    posterior::as_draws_matrix(fit$unconstrain_draws(draws = fit$draws()))
  )

  # Reference point: unconstrained = 0.
  u_0 <- numeric(n_upars)
  # Target point: seeded normal draw (deterministic, independent of HMC).
  set.seed(as.integer(seed))
  u_1 <- rnorm(n_upars, 0, 0.5)

  logp_stan_0 <- as.numeric(fit$log_prob(u_0, jacobian = FALSE))
  logp_stan   <- as.numeric(fit$log_prob(u_1, jacobian = FALSE))

  .clean_cv <- function(u) {
    cv <- fit$constrain_variables(
      unconstrained_variables = u,
      transformed_parameters  = FALSE,
      generated_quantities    = FALSE
    )
    cv <- cv[vapply(cv, function(x) !is.null(x), logical(1))]
    cv[!(names(cv) == "" | grepl("^\\.", names(cv)))]
  }
  cv   <- .clean_cv(u_1)
  cv_0 <- .clean_cv(u_0)

  py_file <- system.file("python", "coev_pymc_logprob.py", package = "coevolve")
  if (!nzchar(py_file)) {
    stop2("inst/python/coev_pymc_logprob.py not found (broken package install?).")
  }
  logp_mod <- reticulate::import_from_path(
    tools::file_path_sans_ext(basename(py_file)),
    path = normalizePath(dirname(py_file), winslash = "/", mustWork = TRUE),
    convert = TRUE
  )

  py_data     <- convert_r_to_python_data_pymc(sd_pymc)
  stan_py     <- reticulate::r_to_py(cv,   convert = TRUE)
  stan_ref_py <- reticulate::r_to_py(cv_0, convert = TRUE)

  res_py <- logp_mod$pymc_logprob_at_stan_primary_params(
    py_data, stan_py, compile_mode, stan_ref_py
  )
  logp_pymc   <- as.numeric(res_py[["logp_pymc"]])
  logp_pymc_0 <- as.numeric(res_py[["logp_pymc_0"]])
  stan_tip_z  <- as.numeric(res_py[["stan_tip_z_log_prior"]])
  stan_tdrift <- as.numeric(res_py[["stan_terminal_drift_log_prior"]])

  # norm_const cancels per-distribution normalization constants (PyMC
  # normalized vs Stan unnormalized).  Both reference evaluations use
  # Stan's constrained values at unconstrained = 0 so likelihood terms
  # cancel exactly.
  norm_const  <- logp_pymc_0 - logp_stan_0
  lp_stan_adj <- logp_stan - stan_tip_z - stan_tdrift
  lp_pymc_adj <- logp_pymc - norm_const

  diff <- abs(lp_stan_adj - lp_pymc_adj)
  if (diff > tol_warn) {
    warning(
      "Large |log p Stan (adj) - log p PyMC (adj)| = ", format(diff, digits = 5),
      " (adj subtracts Stan-only priors on tip-edge z_drift and observed-normal",
      " terminal_drift, plus normalization constant offset; see",
      " ?compare_stan_pymc_logprob).",
      call. = FALSE
    )
  }

  list(
    logprob_stan                  = logp_stan,
    logprob_stan_adj              = lp_stan_adj,
    stan_tip_z_log_prior          = stan_tip_z,
    stan_terminal_drift_log_prior = stan_tdrift,
    logprob_pymc                  = lp_pymc_adj,
    abs_diff                      = diff,
    n_params                      = n_upars,
    stan_params                   = names(cv)
  )
}
