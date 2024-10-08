test_that("coev_pred_series() produces expected errors and output", {
  # load models
  m1 <- readRDS(test_path("fixtures", "coevfit_example1.rds"))
  m2 <- readRDS(test_path("fixtures", "coevfit_example2.rds"))
  m3 <- readRDS(test_path("fixtures", "coevfit_example3.rds"))
  m4 <- readRDS(test_path("fixtures", "coevfit_example4.rds"))
  m5 <- readRDS(test_path("fixtures", "coevfit_example5.rds"))
  m6 <- readRDS(test_path("fixtures", "coevfit_example6.rds"))
  m7 <- readRDS(test_path("fixtures", "coevfit_example7.rds"))
  m8 <- readRDS(test_path("fixtures", "coevfit_example8.rds"))
  m9 <- readRDS(test_path("fixtures", "coevfit_example9.rds"))
  m1 <- reload_fit(m1, filename = "coevfit_example1-1.csv")
  m2 <- reload_fit(m2, filename = "coevfit_example2-1.csv")
  m3 <- reload_fit(m3, filename = "coevfit_example3-1.csv")
  m4 <- reload_fit(m4, filename = "coevfit_example4-1.csv")
  m5 <- reload_fit(m5, filename = "coevfit_example5-1.csv")
  m6 <- reload_fit(m6, filename = "coevfit_example6-1.csv")
  m7 <- reload_fit(m7, filename = "coevfit_example7-1.csv")
  m8 <- reload_fit(m8, filename = "coevfit_example8-1.csv")
  m9 <- reload_fit(m9, filename = "coevfit_example9-1.csv")
  # expect the following errors
  expect_error(
    coev_pred_series(object = "fail"),
    "Argument 'object' must be a fitted coevolutionary model of class coevfit.",
    fixed = TRUE
  )
  expect_error(
    coev_pred_series(object = m1, eta_anc = c(1, "LCA")),
    paste0(
      "Argument 'eta_anc' must be numeric and equal in length to the number ",
      "of variables."
    ),
    fixed = TRUE
  )
  expect_error(
    coev_pred_series(object = m2, eta_anc = c(1, 2)),
    "Argument 'eta_anc' is not a named vector.",
    fixed = TRUE
  )
  expect_error(
    coev_pred_series(object = m2, eta_anc = c("a" = 1, "b" = 2)),
    paste0(
      "Argument 'eta_anc' has names different to the variables included ",
      "in the model."
    ),
    fixed = TRUE
  )
  expect_error(
    coev_pred_series(object = m1, tmax = -1),
    "Argument 'tmax' must be positive.",
    fixed = TRUE
  )
  expect_error(
    coev_pred_series(object = m1, ndraws = "fail"),
    "Argument 'ndraws' must be a single integer.",
    fixed = TRUE
  )
  expect_error(
    coev_pred_series(object = m1, ndraws = 0L),
    "Argument 'ndraws' must be between 1 and the total number of draws.",
    fixed = TRUE
  )
  expect_error(
    coev_pred_series(
      object = m1, ndraws = as.integer(nrow(m1$fit$draws()) + 1)
      ),
    "Argument 'ndraws' must be between 1 and the total number of draws.",
    fixed = TRUE
  )
  expect_error(
    coev_pred_series(object = m1, stochastic = "only"),
    "Argument 'stochastic' must be logical.",
    fixed = TRUE
  )
  # suppress warnings
  SW <- suppressWarnings
  # should run without error
  expect_no_error(SW(coev_pred_series(m1)))
  expect_no_error(SW(coev_pred_series(m2)))
  expect_no_error(SW(coev_pred_series(m3)))
  expect_no_error(SW(coev_pred_series(m4)))
  expect_no_error(SW(coev_pred_series(m5)))
  expect_no_error(SW(coev_pred_series(m6)))
  expect_no_error(SW(coev_pred_series(m7)))
  expect_no_error(SW(coev_pred_series(m8)))
  expect_no_error(SW(coev_pred_series(m9)))
  expect_no_error(SW(coev_pred_series(m1, stochastic = TRUE)))
  expect_no_error(SW(coev_pred_series(m2, stochastic = TRUE)))
  expect_no_error(SW(coev_pred_series(m3, stochastic = TRUE)))
  expect_no_error(SW(coev_pred_series(m4, stochastic = TRUE)))
  expect_no_error(SW(coev_pred_series(m5, stochastic = TRUE)))
  expect_no_error(SW(coev_pred_series(m6, stochastic = TRUE)))
  expect_no_error(SW(coev_pred_series(m7, stochastic = TRUE)))
  expect_no_error(SW(coev_pred_series(m8, stochastic = TRUE)))
  expect_no_error(SW(coev_pred_series(m9, stochastic = TRUE)))
  expect_no_error(SW(coev_pred_series(m1, ndraws = 1L)))
  expect_no_error(SW(coev_pred_series(m2, ndraws = 1L)))
  expect_no_error(SW(coev_pred_series(m3, ndraws = 1L)))
  expect_no_error(SW(coev_pred_series(m4, ndraws = 1L)))
  expect_no_error(SW(coev_pred_series(m5, ndraws = 1L)))
  expect_no_error(SW(coev_pred_series(m6, ndraws = 1L)))
  expect_no_error(SW(coev_pred_series(m7, ndraws = 1L)))
  expect_no_error(SW(coev_pred_series(m8, ndraws = 1L)))
  expect_no_error(SW(coev_pred_series(m9, ndraws = 1L)))
})
