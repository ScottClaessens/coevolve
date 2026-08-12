#' @srrstats {G5.10} Flag extended tests
run_extended_tests <- identical(Sys.getenv("COEVOLVE_EXTENDED_TESTS"), "true")

#' @srrstats {G5.7} Standard errors for posterior estimates are smaller as
#'   number of observations increases
#' @srrstats {BS7.3} Algorithm scaling test
test_that("coev_fit() scales with increasing observations", {
  skip_if_not(run_extended_tests)
  # get dataset with repeated measures
  withr::with_seed(1, {
    n <- 10
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = tree$tip.label,
      x = rnorm(n),
      y = rnorm(n)
    )
    d <- rbind(d, d, d)
  })
  # fit model
  fit_small <-
    coev_fit(
      data = d,
      variables = list(
        x = "normal",
        y = "normal"
      ),
      id = "id",
      tree = tree,
      estimate_residual = FALSE,
      chains = 1,
      refresh = 0,
      seed = 1
    )
  s_small <- fit_small$fit$summary()
  # fit model with 5x the number of observations
  fit_large <-
    coev_fit(
      data = rbind(d, d, d, d, d),
      variables = list(
        x = "normal",
        y = "normal"
      ),
      id = "id",
      tree = tree,
      estimate_residual = FALSE,
      chains = 1,
      refresh = 0,
      seed = 1
    )
  s_large <- fit_large$fit$summary()
  # check models converged
  expect_lt(
    max(suppressWarnings(fit_small$fit$summary())$rhat, na.rm = TRUE),
    1.1
  )
  expect_lt(
    max(suppressWarnings(fit_large$fit$summary())$rhat, na.rm = TRUE),
    1.1
  )
  # standard errors should be smaller with more observations
  # autoregressive effects
  expect_lt(as.numeric(s_large[1, "sd"]), as.numeric(s_small[1, "sd"]))
  expect_lt(as.numeric(s_large[2, "sd"]), as.numeric(s_small[2, "sd"]))
  # cross effects
  expect_lt(as.numeric(s_large[3, "sd"]), as.numeric(s_small[3, "sd"]))
  expect_lt(as.numeric(s_large[4, "sd"]), as.numeric(s_small[4, "sd"]))
  # sd drift
  expect_lt(as.numeric(s_large[10, "sd"]), as.numeric(s_small[10, "sd"]))
  expect_lt(as.numeric(s_large[11, "sd"]), as.numeric(s_small[11, "sd"]))
  # larger data set takes longer to converge
  expect_lt(
    fit_small$fit$time()$total,
    fit_large$fit$time()$total
  )
})
