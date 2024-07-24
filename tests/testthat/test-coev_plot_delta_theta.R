test_that("coev_plot_delta_theta() produces expected errors and output", {
  # load models
  m1 <- readRDS(test_path("fixtures", "coevfit_example1.rds"))
  m2 <- readRDS(test_path("fixtures", "coevfit_example2.rds"))
  m3 <- readRDS(test_path("fixtures", "coevfit_example3.rds"))
  m4 <- readRDS(test_path("fixtures", "coevfit_example4.rds"))
  m5 <- readRDS(test_path("fixtures", "coevfit_example5.rds"))
  m6 <- readRDS(test_path("fixtures", "coevfit_example6.rds"))
  m1 <- reload_fit(m1, filename = "coevfit_example1-1.csv")
  m2 <- reload_fit(m2, filename = "coevfit_example2-1.csv")
  m3 <- reload_fit(m3, filename = "coevfit_example3-1.csv")
  m4 <- reload_fit(m4, filename = "coevfit_example4-1.csv")
  m5 <- reload_fit(m5, filename = "coevfit_example5-1.csv")
  m6 <- reload_fit(m6, filename = "coevfit_example6-1.csv")
  # expect the following errors
  expect_error(
    coev_plot_delta_theta(object = "fail"),
    "Argument 'object' must be a fitted coevolutionary model of class coevfit."
  )
  expect_error(
    coev_plot_delta_theta(object = m1, variables = NA),
    "Argument 'variables' must be a character vector."
  )
  expect_error(
    coev_plot_delta_theta(object = m1, variables = "fail"),
    "Argument 'variables' must be of length > 1."
  )
  expect_error(
    coev_plot_delta_theta(object = m1, variables = c("x", "y", "fail")),
    "Some variables in 'variables' are not included in the fitted model."
  )
  expect_error(
    coev_plot_delta_theta(object = m1, variables = c("x", "y", "y")),
    "Argument 'variables' contains duplicates."
  )
  # suppress warnings
  SW <- suppressWarnings
  # should run without error and produce ggplot object
  expect_no_error(SW(coev_plot_delta_theta(m1)))
  expect_no_error(SW(coev_plot_delta_theta(m2)))
  expect_no_error(SW(coev_plot_delta_theta(m3)))
  expect_no_error(SW(coev_plot_delta_theta(m4)))
  expect_no_error(SW(coev_plot_delta_theta(m5)))
  expect_no_error(SW(coev_plot_delta_theta(m6)))
  expect_true(methods::is(SW(coev_plot_delta_theta(m1)), "ggplot"))
  expect_true(methods::is(SW(coev_plot_delta_theta(m2)), "ggplot"))
  expect_true(methods::is(SW(coev_plot_delta_theta(m3)), "ggplot"))
  expect_true(methods::is(SW(coev_plot_delta_theta(m4)), "ggplot"))
  expect_true(methods::is(SW(coev_plot_delta_theta(m5)), "ggplot"))
  expect_true(methods::is(SW(coev_plot_delta_theta(m6)), "ggplot"))
})
