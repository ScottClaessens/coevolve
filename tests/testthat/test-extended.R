#' @srrstats {G5.10} Flag extended tests
run_extended_tests <- identical(Sys.getenv("COEVOLVE_EXTENDED_TESTS"), "true")

#' @srrstats {G5.4, G5.4a, G5.5} Extended correctness test with fixed seed
test_that("coev_fit() estimates correct direction of selection", {
  skip_if_not(run_extended_tests)
  # simulate data with fixed seed, where x -> y but not vice versa
  withr::with_seed(1, {
    sim <-
      coev_simulate_coevolution(
        n = 100,
        variables = c("x", "y"),
        selection_matrix = matrix(
          c(0.9, 0.0,
            0.9, 0.9),
          nrow = 2,
          byrow = TRUE,
          dimnames = list(c("x", "y"), c("x", "y"))
        ),
        drift = c("x" = 0.01, "y" = 0.01),
        prob_split = 0.05
      )
  })
  # fit model
  fit <-
    coev_fit(
      data = sim$data,
      variables = list(
        x = "normal",
        y = "normal"
      ),
      id = "species",
      tree = sim$tree,
      estimate_correlated_drift = FALSE,
      chains = 1,
      seed = 1,
      refresh = 0
    )
  # get delta theta values
  delta_theta_x_to_y <-
    coev_calculate_delta_theta(
      object = fit,
      predictor = "x",
      response = "y"
    )
  delta_theta_y_to_x <-
    coev_calculate_delta_theta(
      object = fit,
      predictor = "y",
      response = "x"
    )
  # for x -> y, 95% CI for delta theta should be greater than zero
  expect_true(
    quantile(delta_theta_x_to_y, 0.025) > 0 &&
      quantile(delta_theta_x_to_y, 0.975) > 0
  )
  # for y -> x, 95% CI for delta theta should include zero
  expect_true(
    quantile(delta_theta_y_to_x, 0.025) < 0 &&
      quantile(delta_theta_y_to_x, 0.975) > 0
  )
})

#' @srrstats {G5.6, G5.6a, G5.6b, G5.9, G5.9a, G5.9b} Extended parameter
#'   recovery tests with multiple fixed seeds for data simulation and cmdstanr
for (seed in 1:3) {
  test_that(paste0("coev_fit() recovers parameters (seed = ", seed, ")"), {
    skip_if_not(run_extended_tests)
    # get dummy data
    withr::with_seed(seed, {
      n <- 50
      tree <- ape::rcoal(n)
      d <- data.frame(
        id = tree$tip.label,
        x = rnorm(n),
        y = rnorm(n)
      )
    })
    # get stan data list with prior_only
    sdata <- coev_make_standata(
      data = d,
      variables = list(
        x = "normal",
        y = "normal"
      ),
      id = "id",
      tree = tree,
      estimate_correlated_drift = FALSE,
      prior_only = TRUE
    )
    # manually fix parameters by editing the stan code
    scode <- coev_make_stancode(
      data = d,
      variables = list(
        x = "normal",
        y = "normal"
      ),
      id = "id",
      tree = tree,
      estimate_correlated_drift = FALSE
    )
    scode_fixed <- manually_fix_parameters(scode)
    # simulate dataset from model with fixed parameters
    sim <- cmdstanr::cmdstan_model(
      stan_file = cmdstanr::write_stan_file(scode_fixed)
    )$sample(
      data = sdata,
      chains = 1,
      refresh = 0,
      seed = seed,
      iter_warmup = 50,
      iter_sampling = 1
    )
    draws <- posterior::as_draws_rvars(sim)
    d_sim <- data.frame(
      id = tree$tip.label,
      x = posterior::draws_of(draws$yrep)[1, 1, 1:n, 1],
      y = posterior::draws_of(draws$yrep)[1, 1, 1:n, 2]
    )
    # fit model to simulated data
    effects_mat <- matrix(TRUE, nrow = 2, ncol = 2,
                          dimnames = list(c("x", "y"), c("x", "y")))
    effects_mat[1, 2] <- FALSE
    fit <- coev_fit(
      data = d_sim,
      variables = list(
        x = "normal",
        y = "normal"
      ),
      id = "id",
      tree = tree,
      estimate_correlated_drift = FALSE,
      effects_mat = effects_mat,
      prior = list(
        A_offdiag = "normal(0, 2)",
        Q_sigma = "normal(0, 2)"
      ),
      scale = FALSE,
      chains = 1,
      refresh = 0,
      seed = seed,
      max_treedepth = 15
    )
    post <- extract_samples(fit)
    # posterior medians should recover fixed parameters within tolerance
    expect_equal(median(post$A[, 1, 1]), -0.5, tolerance = 1)
    expect_equal(median(post$A[, 2, 2]), -0.5, tolerance = 1)
    expect_equal(median(post$A[, 2, 1]),  3.0, tolerance = 1)
    expect_equal(median(post$Q[, 1, 1]),  1.5, tolerance = 1)
    expect_equal(median(post$Q[, 2, 2]),  1.5, tolerance = 1)
    expect_equal(median(post$b[, 1]), 0.0, tolerance = 1)
    expect_equal(median(post$b[, 2]), 0.0, tolerance = 1)
  })
}

#' @srrstats {G5.7} Standard errors for posterior estimates are smaller as
#'   number of observations increases
test_that("coev_fit() estimates smaller posterior SEs with more observations", {
  skip_if_not(run_extended_tests)
  # get smaller dataset
  withr::with_seed(1, {
    n_small <- 50
    tree_small <- ape::rcoal(n_small)
    d_small <- data.frame(
      id = tree_small$tip.label,
      x = rnorm(n_small),
      y = rnorm(n_small)
    )
  })
  # fit model to smaller dataset
  fit_small <-
    coev_fit(
      data = d_small,
      variables = list(
        x = "normal",
        y = "normal"
      ),
      id = "id",
      tree = tree_small,
      chains = 1,
      refresh = 0,
      seed = 1
    )
  s_small <- suppressWarnings(summary(fit_small))
  # get larger dataset
  withr::with_seed(1, {
    n_large <- 100
    tree_large <- ape::rcoal(n_large)
    d_large <- data.frame(
      id = tree_large$tip.label,
      x = rnorm(n_large),
      y = rnorm(n_large)
    )
  })
  # fit model to larger dataset
  fit_large <-
    coev_fit(
      data = d_large,
      variables = list(
        x = "normal",
        y = "normal"
      ),
      id = "id",
      tree = tree_large,
      chains = 1,
      refresh = 0,
      seed = 1
    )
  s_large <- suppressWarnings(summary(fit_large))
  # standard errors should be smaller with more observations
  expect_lt(s_large$auto[1, "Est.Error"], s_small$auto[1, "Est.Error"])
  expect_lt(s_large$auto[2, "Est.Error"], s_small$auto[2, "Est.Error"])
  expect_lt(s_large$cross[1, "Est.Error"], s_small$cross[1, "Est.Error"])
  expect_lt(s_large$cross[2, "Est.Error"], s_small$cross[2, "Est.Error"])
  expect_lt(s_large$sd_drift[1, "Est.Error"], s_small$sd_drift[1, "Est.Error"])
  expect_lt(s_large$sd_drift[2, "Est.Error"], s_small$sd_drift[2, "Est.Error"])
})
