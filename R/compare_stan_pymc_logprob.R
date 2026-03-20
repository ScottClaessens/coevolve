#' Compare Stan and PyMC log density at the same parameters
#'
#' Draws one post-warmup posterior draw with CmdStan, maps its **primary**
#' (constrained) parameters to the corresponding PyMC free random variables, and
#' evaluates \code{compile_logp(jacobian = TRUE)} on the PyMC model. This checks
#' that the two implementations target the same posterior \strong{at points
#' reachable by both parameterizations} (diagnostic; not a proof of global
#' equivalence).
#'
#' **Scope:** Use \code{estimate_correlated_drift = FALSE} (diagonal \eqn{Q}) so
#' there is no LKJ Cholesky reparameterization mismatch. Stan stores extra
#' \code{z_drift} components on tip edges; those priors are subtracted into
#' \code{stan_tip_z_log_prior}. If Stan declares parameters absent from the PyMC
#' model (e.g. \code{terminal_drift} when PyMC folds them into the Gaussian
#' likelihood for all-normal, fully-observed data), totals will not match---use
#' at least one missing Gaussian response so both include \code{terminal_drift},
#' or interpret differences carefully.
#'
#' Intended use: one-off debugging and regression tests when Stan and PyMC
#' timings or posteriors disagree.
#'
#' @inheritParams coev_make_stancode
#' @param iter_warmup,iter_sampling Passed to \code{cmdstanr::sample()}. A short
#'   run is enough to obtain one valid unconstrained draw in support.
#' @param chains Number of chains (default \code{1}).
#' @param seed Integer seed for CmdStan and for reproducibility.
#' @param compile_mode PyTensor mode passed to \code{build_model()} (\code{"cpu"}
#'   or \code{"mlx"}).
#' @param tol_warn Absolute difference \code{abs(logp_stan - logp_pymc)} above
#'   which a warning is issued (default \code{0.05} on log scale).
#'
#' @returns A list with \code{logprob_stan} (from \code{fit$log_prob}),
#'   \code{stan_tip_z_log_prior} (independent normal log density for Stan
#'   \code{z_drift} components that have no PyMC counterpart),
#'   \code{stan_terminal_drift_log_prior} (std_normal log density for Stan
#'   \code{terminal_drift} elements when PyMC omits that parameter),
#'   \code{logprob_stan_adj}
#'   (= \code{logprob_stan - stan_tip_z_log_prior - stan_terminal_drift_log_prior}),
#'   \code{logprob_pymc},
#'   \code{abs_diff} (= \code{|logprob_stan_adj - logprob_pymc|}), \code{draw_idx},
#'   \code{n_params}, and \code{stan_params}.
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
    iter_warmup = 150L,
    iter_sampling = 1L,
    chains = 1L,
    seed = 1L,
    compile_mode = "cpu",
    tol_warn = 0.05) {

  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop2("Package 'cmdstanr' is required for compare_stan_pymc_logprob().")
  }
  stop_if_pymc_not_available()

  distributions <- as.character(variables)

  sc <- coev_make_stancode(
    data, variables, id, tree, effects_mat, complete_cases, dist_mat,
    dist_cov, measurement_error, prior, scale,
    estimate_correlated_drift, estimate_residual,
    log_lik = FALSE, prior_only = prior_only
  )
  pc <- coev_make_pymc(
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
  sd_pymc <- standata_to_pymc(sd, distributions)

  mod <- cmdstanr::cmdstan_model(
    cmdstanr::write_stan_file(sc),
    force_recompile = FALSE
  )

  fit <- mod$sample(
    data = sd,
    chains = as.integer(chains),
    parallel_chains = 1L,
    iter_warmup = as.integer(iter_warmup),
    iter_sampling = as.integer(iter_sampling),
    seed = as.integer(seed),
    refresh = 0,
    show_messages = FALSE,
    init = 0,
    adapt_delta = 0.95
  )

  udraws <- fit$unconstrain_draws(draws = fit$draws())
  um <- posterior::as_draws_matrix(udraws)
  u <- as.numeric(um[1, ])

  logp_stan <- as.numeric(fit$log_prob(u, jacobian = TRUE))

  cv <- fit$constrain_variables(
    unconstrained_variables = u,
    transformed_parameters = FALSE,
    generated_quantities = FALSE
  )
  cv <- cv[vapply(cv, function(x) !is.null(x), logical(1))]
  bad <- names(cv) == "" | grepl("^\\.", names(cv))
  if (any(bad)) {
    cv <- cv[!bad]
  }

  py_file <- system.file("python", "coev_pymc_logprob.py", package = "coevolve")
  if (!nzchar(py_file)) {
    stop2("inst/python/coev_pymc_logprob.py not found (broken package install?).")
  }
  logp_mod <- reticulate::import_from_path(
    tools::file_path_sans_ext(basename(py_file)),
    path = normalizePath(dirname(py_file), winslash = "/", mustWork = TRUE),
    convert = TRUE
  )

  py_data <- convert_r_to_python_data_pymc(sd_pymc)
  stan_py <- reticulate::r_to_py(cv, convert = TRUE)

  res_py <- logp_mod$pymc_logprob_at_stan_primary_params(
    pc,
    py_data,
    stan_py,
    compile_mode
  )
  logp_pymc <- as.numeric(res_py[["logp_pymc"]])
  stan_tip_z <- as.numeric(res_py[["stan_tip_z_log_prior"]])
  stan_tdrift <- as.numeric(res_py[["stan_terminal_drift_log_prior"]])
  lp_stan_adj <- logp_stan - stan_tip_z - stan_tdrift

  diff <- abs(lp_stan_adj - logp_pymc)
  if (diff > tol_warn) {
    warning(
      "Large |log p Stan (adj) - log p PyMC| = ", format(diff, digits = 5),
      " (adj subtracts Stan-only priors on tip-edge z_drift and observed-normal",
      " terminal_drift; see ?compare_stan_pymc_logprob).",
      call. = FALSE
    )
  }

  list(
    logprob_stan = logp_stan,
    logprob_stan_adj = lp_stan_adj,
    stan_tip_z_log_prior = stan_tip_z,
    stan_terminal_drift_log_prior = stan_tdrift,
    logprob_pymc = logp_pymc,
    abs_diff = diff,
    draw_idx = 1L,
    n_params = length(u),
    stan_params = names(cv)
  )
}
