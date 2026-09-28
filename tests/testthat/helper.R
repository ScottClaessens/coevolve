# reload fit from cmdstan csv in fixtures folder
reload_fit <- function(coevfit, filename) {
  coevfit$fit <-
    cmdstanr::as_cmdstan_fit(
      testthat::test_path("fixtures", filename)
    )
  coevfit
}

# Setup nutpie for testing
# This function configures reticulate to use a Python environment with nutpie.
# It requires explicit configuration via environment variables.
#
# Configuration options:
# 1. Set NUTPIE_PYTHON environment variable to point to Python executable:
#    export NUTPIE_PYTHON=/path/to/python
#    or in R: Sys.setenv(NUTPIE_PYTHON = "/path/to/python")
#
# 2. Set NUTPIE_VENV environment variable to point to virtual environment:
#    export NUTPIE_VENV=~/.venvs/nutpie-env
#    or in R: Sys.setenv(NUTPIE_VENV = "~/.venvs/nutpie-env")
#
# 3. Set RETICULATE_PYTHON environment variable (reticulate's standard
#    variable): export RETICULATE_PYTHON=/path/to/python
#
# Note: This function does NOT search for nutpie automatically. Explicit
# configuration is required to ensure reproducibility and avoid version
# conflicts.
setup_nutpie_for_tests <- function() {
  # only try if reticulate is available
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    return(FALSE)
  }
  # check NUTPIE_PYTHON environment variable (highest priority)
  nutpie_python <- Sys.getenv("NUTPIE_PYTHON", unset = "")
  if (nutpie_python != "" && file.exists(nutpie_python)) {
    tryCatch({
      # set RETICULATE_PYTHON before reticulate initializes
      Sys.setenv(RETICULATE_PYTHON = nutpie_python)
      # configure reticulate if Python hasn't been initialized yet
      if (!reticulate::py_available()) {
        reticulate::use_python(nutpie_python, required = FALSE)
      } else {
        # try to reconfigure (may not work if already initialized)
        reticulate::use_python(nutpie_python, required = FALSE)
      }
      if (coevolve:::check_nutpie_available()) {
        return(TRUE)
      }
    }, error = function(e) NULL)
  }
  # check NUTPIE_VENV environment variable
  nutpie_venv <- Sys.getenv("NUTPIE_VENV", unset = "")
  if (nutpie_venv != "") {
    expanded_venv <- path.expand(nutpie_venv)
    if (dir.exists(expanded_venv)) {
      tryCatch({
        reticulate::use_virtualenv(expanded_venv, required = FALSE)
        if (coevolve:::check_nutpie_available()) {
          return(TRUE)
        }
      }, error = function(e) NULL)
    }
  }
  # check RETICULATE_PYTHON (reticulate's standard variable)
  reticulate_python <- Sys.getenv("RETICULATE_PYTHON", unset = "")
  if (reticulate_python != "" && file.exists(reticulate_python)) {
    # if RETICULATE_PYTHON is set, reticulate should use it automatically
    # just check if nutpie is available
    if (coevolve:::check_nutpie_available()) {
      return(TRUE)
    }
  }
  # if nutpie is already available (from previous configuration), return TRUE
  if (coevolve:::check_nutpie_available()) {
    return(TRUE)
  }
  # if we get here, nutpie is not available
  FALSE
}

# manually fix parameters in stan code
manually_fix_parameters <- function(scode) {
  scode |>
    stringr::str_remove(
      stringr::fixed(
        paste0(
          "  vector<upper=0>[J] A_diag; // autoregressive terms of A\n",
          "  vector[num_effects - J] A_offdiag; // cross-lagged terms of A\n",
          "  vector<lower=0>[J] Q_sigma; // std deviation parameters of the ",
          "Q mat\n",
          "  vector[J] b; // SDE intercepts\n",
          "  array[N_tree] vector[J] eta_anc; // ancestral states\n"
        )
      )
    ) |>
    stringr::str_replace(
      pattern = stringr::fixed(
        paste0(
          "  matrix[J,J] A = diag_matrix(A_diag); // selection matrix\n",
          "  matrix[J,J] Q = diag_matrix(Q_sigma^2); // drift matrix\n"
        )
      ),
      replacement = paste0(
        "  array[N_tree] vector[J] eta_anc;\n",
        "  vector[J] b = rep_vector(0.0, J);\n",
        "  matrix[J,J] A = diag_matrix(rep_vector(-0.5, J));\n",
        "  matrix[J,J] Q = diag_matrix(rep_vector(1.5, J));\n"
      )
    ) |>
    stringr::str_replace(
      pattern = stringr::fixed(
        paste0(
          "  // fill off diagonal of A matrix\n",
          "  {\n",
          "    int ticker = 1;\n",
          "    for (i in 1:J) {\n",
          "      for (j in 1:J) {\n",
          "        if (i != j) {\n",
          "          if (effects_mat[i,j] == 1) {\n",
          "            A[i,j] = A_offdiag[ticker];\n",
          "            ticker += 1;\n",
          "          } else if (effects_mat[i,j] == 0) {\n",
          "            A[i,j] = 0;\n",
          "          }\n",
          "        }\n",
          "      }\n",
          "    }\n",
          "  }\n"
        )
      ),
      replacement = paste0(
        "  A[2,1] = 1;\n",
        "  A[1,2] = 0;\n",
        "  for (t in 1:N_tree) eta_anc[t] = rep_vector(0.0, J);\n"
      )
    ) |>
    stringr::str_remove(stringr::fixed("  b ~ std_normal();\n")) |>
    stringr::str_remove(stringr::fixed("    eta_anc[t] ~ std_normal();\n")) |>
    stringr::str_remove(stringr::fixed("  A_offdiag ~ std_normal();\n")) |>
    stringr::str_remove(stringr::fixed("  A_diag ~ std_normal();\n")) |>
    stringr::str_remove(stringr::fixed("  Q_sigma ~ std_normal();\n"))
}

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
