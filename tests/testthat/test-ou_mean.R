# The expected trait values at the tips must follow the exact solution of the
# Ornstein-Uhlenbeck mean equation dm/dt = A m + b. With stochastic drift set
# to zero and an ultrametric tree of height one, every tip should equal
# mu + exp(A) (eta_anc - mu), where mu = -A^-1 b. This holds for any A,
# including asymmetric A with cross-selection effects and non-zero b.
test_that("tip means follow the exact OU solution for asymmetric A", {
  skip_on_cran()
  withr::with_seed(1, {
    n <- 10
    tree <- ape::rcoal(n)
    d <- data.frame(id = tree$tip.label, x = rnorm(n), y = rnorm(n))
  })
  variables <- list(x = "normal", y = "normal")
  mod <- cmdstanr::cmdstan_model(
    cmdstanr::write_stan_file(coev_make_stancode(d, variables, "id", tree)),
    compile_model_methods = TRUE,
    force_recompile = TRUE
  )
  stan_data <- coev_make_standata(d, variables, "id", tree)
  fit <- suppressWarnings(mod$sample(
    data = stan_data, chains = 1, iter_warmup = 1, iter_sampling = 1,
    refresh = 0, show_messages = FALSE, seed = 1, init = 0,
    diagnostics = NULL
  ))
  fit$init_model_methods()
  n_upars <- ncol(posterior::as_draws_matrix(
    fit$unconstrain_draws(draws = fit$draws())
  ))
  withr::with_seed(2, {
    for (k in 1:3) {
      # random parameter values with stochastic drift switched off
      pars <- fit$constrain_variables(rnorm(n_upars))
      pars$z_drift <- pars$z_drift * 0
      cons <- fit$constrain_variables(
        fit$unconstrain_variables(pars),
        transformed_parameters = TRUE
      )
      a <- matrix(cons$A, 2, 2)
      expect_gt(abs(a[1, 2] - a[2, 1]), 0.01)
      mu <- -solve(a, cons$b)
      exp_a <- as.matrix(Matrix::expm(Matrix::Matrix(a)))
      expected <- mu + exp_a %*% (as.vector(cons$eta_anc) - mu)
      eta_tips <- matrix(cons$eta, stan_data$N_seg, 2)[seq_len(n), ]
      expect_equal(
        eta_tips,
        matrix(expected, n, 2, byrow = TRUE),
        tolerance = 1e-8
      )
    }
  })
})
