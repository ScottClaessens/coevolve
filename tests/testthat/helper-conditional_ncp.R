# Helpers for checking that the conditionally non-centred terminal drift
# parameterisation defines exactly the same model as the centred one.
#
# The centred Stan code is regenerated from the same generator by mocking
# use_conditional_ncp() to FALSE. Points are drawn in the non-centred space
# and mapped to the centred space with an independent R implementation of
# the transform, so both the log density (up to the change-of-variables
# Jacobian) and the pointwise log likelihood can be compared directly.

# compile a model exposing log_prob and constrain/unconstrain methods
compile_with_methods <- function(stan_code, stan_data) {
  mod <- cmdstanr::cmdstan_model(
    cmdstanr::write_stan_file(stan_code),
    compile_model_methods = TRUE,
    force_recompile = TRUE
  )
  fit <- suppressWarnings(mod$sample(
    data = stan_data, chains = 1, iter_warmup = 1, iter_sampling = 1,
    refresh = 0, show_messages = FALSE, seed = 1, init = 0,
    diagnostics = NULL
  ))
  fit$init_model_methods()
  fit
}

# map a constrained non-centred point to the centred parameterisation,
# returning the centred terminal drift and the log absolute Jacobian
ncp_to_centred <- function(cons, stan_data, distributions,
                           measurement_error = FALSE) {
  n_tree <- stan_data$N_tree
  n_seg <- stan_data$N_seg
  n_tips <- stan_data$N_tips
  n_var <- stan_data$J
  vcv <- array(cons$VCV_tips, c(n_tree, n_seg, n_var, n_var))
  eta <- array(cons$eta, c(n_tree, n_seg, n_var))
  dist_v <- if (!is.null(cons$dist_v)) {
    array(cons$dist_v, c(n_tips, n_var))
  } else {
    array(0, c(n_tips, n_var))
  }
  terminal_drift <- array(cons$terminal_drift, c(n_tree, n_tips, n_var))
  is_normal <- distributions == "normal"
  log_jacobian <- 0
  for (t in seq_len(n_tree)) {
    for (i in seq_len(stan_data$N_obs)) {
      tip <- stan_data$tip_id[i]
      observed <- is_normal & stan_data$miss[i, ] == 0
      perm <- c(which(observed), which(!observed))
      n_obs <- sum(observed)
      if (n_obs == n_var) next
      sigma <- vcv[t, tip, , ]
      if (measurement_error) sigma <- sigma + diag(stan_data$se[i, ], n_var)
      chol_perm <- t(chol(sigma[perm, perm]))
      latent <- (n_obs + 1):n_var
      z <- terminal_drift[t, tip, perm[latent]]
      drift_latent <- chol_perm[latent, latent, drop = FALSE] %*% z
      if (n_obs > 0) {
        drift_obs <- stan_data$y[i, perm[seq_len(n_obs)]] -
          eta[t, tip, perm[seq_len(n_obs)]] -
          dist_v[tip, perm[seq_len(n_obs)]]
        chol_obs <- chol_perm[seq_len(n_obs), seq_len(n_obs), drop = FALSE]
        drift_latent <- drift_latent +
          chol_perm[latent, seq_len(n_obs), drop = FALSE] %*%
          forwardsolve(chol_obs, drift_obs)
      }
      terminal_drift[t, tip, perm[latent]] <- drift_latent
      log_jacobian <- log_jacobian + sum(log(diag(chol_perm)[latent]))
    }
  }
  list(terminal_drift = terminal_drift, log_jacobian = log_jacobian)
}

# compare centred and conditionally non-centred models at random points
compare_ncp_centred <- function(data, variables, id, tree, ...,
                                n_points = 5L, seed = 1L) {
  args <- list(data = data, variables = variables, id = id, tree = tree, ...)
  args$log_lik <- TRUE
  code_ncp <- do.call(coev_make_stancode, args)
  code_centred <- testthat::with_mocked_bindings(
    do.call(coev_make_stancode, args),
    use_conditional_ncp = function(...) FALSE
  )
  stan_data <- do.call(coev_make_standata, args)
  fit_ncp <- compile_with_methods(code_ncp, stan_data)
  fit_centred <- compile_with_methods(code_centred, stan_data)
  n_upars <- ncol(posterior::as_draws_matrix(
    fit_ncp$unconstrain_draws(draws = fit_ncp$draws())
  ))
  withr::with_seed(seed, {
    res <- lapply(seq_len(n_points), function(k) {
      u <- stats::rnorm(n_upars, sd = 0.5)
      cons <- fit_ncp$constrain_variables(
        u, transformed_parameters = TRUE, generated_quantities = TRUE
      )
      mapped <- ncp_to_centred(
        cons, stan_data, as.character(variables),
        measurement_error = !is.null(args$measurement_error)
      )
      pars_centred <- fit_ncp$constrain_variables(u)
      pars_centred$terminal_drift <- mapped$terminal_drift
      u_centred <- fit_centred$unconstrain_variables(pars_centred)
      cons_centred <- fit_centred$constrain_variables(
        u_centred, generated_quantities = TRUE
      )
      list(
        lp_ncp = fit_ncp$log_prob(u, jacobian = TRUE),
        lp_centred_plus_jacobian =
          fit_centred$log_prob(u_centred, jacobian = TRUE) +
          mapped$log_jacobian,
        log_lik_ncp = cons$log_lik,
        log_lik_centred = cons_centred$log_lik
      )
    })
  })
  list(
    code_ncp = code_ncp,
    code_centred = code_centred,
    lp_ncp = vapply(res, `[[`, numeric(1), "lp_ncp"),
    lp_centred_plus_jacobian =
      vapply(res, `[[`, numeric(1), "lp_centred_plus_jacobian"),
    log_lik_ncp = do.call(rbind, lapply(res, `[[`, "log_lik_ncp")),
    log_lik_centred = do.call(rbind, lapply(res, `[[`, "log_lik_centred"))
  )
}
