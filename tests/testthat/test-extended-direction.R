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
  # check model converged
  expect_lt(max(suppressWarnings(fit$fit$summary())$rhat, na.rm = TRUE), 1.1)
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
