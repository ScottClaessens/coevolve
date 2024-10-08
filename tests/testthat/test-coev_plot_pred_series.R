test_that("coev_plot_pred_series() produces expected errors and output", {
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
    coev_plot_pred_series(object = "fail"),
    "Argument 'object' must be a fitted coevolutionary model of class coevfit.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_pred_series(object = m1, prob = 95),
    "Argument 'prob' is not between 0 and 1.",
    fixed = TRUE
  )
  # suppress warnings
  SW <- suppressWarnings
  # should run without error and produce list of ggplot objects
  expect_no_error(SW(coev_plot_pred_series(m1)))
  expect_no_error(SW(coev_plot_pred_series(m2)))
  expect_no_error(SW(coev_plot_pred_series(m3)))
  expect_no_error(SW(coev_plot_pred_series(m4)))
  expect_no_error(SW(coev_plot_pred_series(m5)))
  expect_no_error(SW(coev_plot_pred_series(m6)))
  expect_no_error(SW(coev_plot_pred_series(m7)))
  expect_no_error(SW(coev_plot_pred_series(m8)))
  expect_no_error(SW(coev_plot_pred_series(m9)))
  expect_no_error(SW(coev_plot_pred_series(m1, stochastic = TRUE)))
  expect_no_error(SW(coev_plot_pred_series(m2, stochastic = TRUE)))
  expect_no_error(SW(coev_plot_pred_series(m3, stochastic = TRUE)))
  expect_no_error(SW(coev_plot_pred_series(m4, stochastic = TRUE)))
  expect_no_error(SW(coev_plot_pred_series(m5, stochastic = TRUE)))
  expect_no_error(SW(coev_plot_pred_series(m6, stochastic = TRUE)))
  expect_no_error(SW(coev_plot_pred_series(m7, stochastic = TRUE)))
  expect_no_error(SW(coev_plot_pred_series(m8, stochastic = TRUE)))
  expect_no_error(SW(coev_plot_pred_series(m9, stochastic = TRUE)))
  expect_no_error(SW(coev_plot_pred_series(m1, ndraws = 1L)))
  expect_no_error(SW(coev_plot_pred_series(m2, ndraws = 1L)))
  expect_no_error(SW(coev_plot_pred_series(m3, ndraws = 1L)))
  expect_no_error(SW(coev_plot_pred_series(m4, ndraws = 1L)))
  expect_no_error(SW(coev_plot_pred_series(m5, ndraws = 1L)))
  expect_no_error(SW(coev_plot_pred_series(m6, ndraws = 1L)))
  expect_no_error(SW(coev_plot_pred_series(m7, ndraws = 1L)))
  expect_no_error(SW(coev_plot_pred_series(m8, ndraws = 1L)))
  expect_no_error(SW(coev_plot_pred_series(m9, ndraws = 1L)))
  expect_true(methods::is(SW(coev_plot_pred_series(m1)), "ggplot"))
  expect_true(methods::is(SW(coev_plot_pred_series(m2)), "ggplot"))
  expect_true(methods::is(SW(coev_plot_pred_series(m3)), "ggplot"))
  expect_true(methods::is(SW(coev_plot_pred_series(m4)), "ggplot"))
  expect_true(methods::is(SW(coev_plot_pred_series(m5)), "ggplot"))
  expect_true(methods::is(SW(coev_plot_pred_series(m6)), "ggplot"))
  expect_true(methods::is(SW(coev_plot_pred_series(m7)), "ggplot"))
  expect_true(methods::is(SW(coev_plot_pred_series(m8)), "ggplot"))
  expect_true(methods::is(SW(coev_plot_pred_series(m9)), "ggplot"))
})
