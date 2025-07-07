test_that("coev_plot_delta_theta() produces expected errors and output", {
  # load models
  m01 <- readRDS(test_path("fixtures", "coevfit_example_01.rds"))
  m02 <- readRDS(test_path("fixtures", "coevfit_example_02.rds"))
  m03 <- readRDS(test_path("fixtures", "coevfit_example_03.rds"))
  m04 <- readRDS(test_path("fixtures", "coevfit_example_04.rds"))
  m05 <- readRDS(test_path("fixtures", "coevfit_example_05.rds"))
  m06 <- readRDS(test_path("fixtures", "coevfit_example_06.rds"))
  m07 <- readRDS(test_path("fixtures", "coevfit_example_07.rds"))
  m08 <- readRDS(test_path("fixtures", "coevfit_example_08.rds"))
  m09 <- readRDS(test_path("fixtures", "coevfit_example_09.rds"))
  m10 <- readRDS(test_path("fixtures", "coevfit_example_10.rds"))
  m01 <- reload_fit(m01, filename = "coevfit_example_01-1.csv")
  m02 <- reload_fit(m02, filename = "coevfit_example_02-1.csv")
  m03 <- reload_fit(m03, filename = "coevfit_example_03-1.csv")
  m04 <- reload_fit(m04, filename = "coevfit_example_04-1.csv")
  m05 <- reload_fit(m05, filename = "coevfit_example_05-1.csv")
  m06 <- reload_fit(m06, filename = "coevfit_example_06-1.csv")
  m07 <- reload_fit(m07, filename = "coevfit_example_07-1.csv")
  m08 <- reload_fit(m08, filename = "coevfit_example_08-1.csv")
  m09 <- reload_fit(m09, filename = "coevfit_example_09-1.csv")
  m10 <- reload_fit(m10, filename = "coevfit_example_10-1.csv")
  # expect the following errors
  #' @srrstats {G5.2, G5.2b} Test all error messages
  expect_error(
    coev_plot_delta_theta(object = "fail"),
    "Argument 'object' must be a fitted coevolutionary model of class coevfit.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_delta_theta(object = m01, variables = NA),
    "Argument 'variables' must be a character vector.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_delta_theta(object = m01, variables = "fail"),
    "Argument 'variables' must be of length > 1.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_delta_theta(object = m01, variables = c("x", "y", "fail")),
    "Some variables in 'variables' are not included in the fitted model.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_delta_theta(object = m01, variables = c("x", "y", "y")),
    "Argument 'variables' contains duplicates.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_delta_theta(object = m01, prob = "fail"),
    "Argument 'prob' must be numeric.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_delta_theta(object = m01, prob = c(0.5, 0.5)),
    "Argument 'prob' must be of length 1.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_delta_theta(object = m01, prob = -0.5),
    "Argument 'prob' must be between 0 and 1.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_delta_theta(object = m01, prob = 1.5),
    "Argument 'prob' must be between 0 and 1.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_delta_theta(object = m01, prob_outer = "fail"),
    "Argument 'prob_outer' must be numeric.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_delta_theta(object = m01, prob_outer = c(0.5, 0.5)),
    "Argument 'prob_outer' must be of length 1.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_delta_theta(object = m01, prob_outer = -0.5),
    "Argument 'prob_outer' must be between 0 and 1.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_delta_theta(object = m01, prob_outer = 1.5),
    "Argument 'prob_outer' must be between 0 and 1.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_delta_theta(object = m01, prob = 0.5, prob_outer = 0.4),
    "Argument 'prob_outer' must be greater than argument 'prob'.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_delta_theta(object = m01, limits = FALSE),
    "Argument 'limits' must be a numeric vector.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_delta_theta(object = m01, limits = 0),
    "Argument 'limits' must be of length 2.",
    fixed = TRUE
  )
  # suppress warnings
  sw <- suppressWarnings
  # should run without error and produce ggplot object
  expect_no_error(sw(coev_plot_delta_theta(m01)))
  expect_no_error(sw(coev_plot_delta_theta(m02)))
  expect_no_error(sw(coev_plot_delta_theta(m03)))
  expect_no_error(sw(coev_plot_delta_theta(m04)))
  expect_no_error(sw(coev_plot_delta_theta(m05)))
  expect_no_error(sw(coev_plot_delta_theta(m06)))
  expect_no_error(sw(coev_plot_delta_theta(m07)))
  expect_no_error(sw(coev_plot_delta_theta(m08)))
  expect_no_error(sw(coev_plot_delta_theta(m09)))
  expect_no_error(sw(coev_plot_delta_theta(m10)))
  expect_true(methods::is(sw(coev_plot_delta_theta(m01)), "ggplot"))
  expect_true(methods::is(sw(coev_plot_delta_theta(m02)), "ggplot"))
  expect_true(methods::is(sw(coev_plot_delta_theta(m03)), "ggplot"))
  expect_true(methods::is(sw(coev_plot_delta_theta(m04)), "ggplot"))
  expect_true(methods::is(sw(coev_plot_delta_theta(m05)), "ggplot"))
  expect_true(methods::is(sw(coev_plot_delta_theta(m06)), "ggplot"))
  expect_true(methods::is(sw(coev_plot_delta_theta(m07)), "ggplot"))
  expect_true(methods::is(sw(coev_plot_delta_theta(m08)), "ggplot"))
  expect_true(methods::is(sw(coev_plot_delta_theta(m09)), "ggplot"))
  expect_true(methods::is(sw(coev_plot_delta_theta(m10)), "ggplot"))
  # limits work as expected
  expect_no_error(sw(coev_plot_delta_theta(m01, limits = c(-5, 5))))
  expect_no_error(sw(coev_plot_delta_theta(m02, limits = c(-5, 5))))
  expect_no_error(sw(coev_plot_delta_theta(m03, limits = c(-5, 5))))
  expect_no_error(sw(coev_plot_delta_theta(m04, limits = c(-5, 5))))
  expect_no_error(sw(coev_plot_delta_theta(m05, limits = c(-5, 5))))
  expect_no_error(sw(coev_plot_delta_theta(m06, limits = c(-5, 5))))
  expect_no_error(sw(coev_plot_delta_theta(m07, limits = c(-5, 5))))
  expect_no_error(sw(coev_plot_delta_theta(m08, limits = c(-5, 5))))
  expect_no_error(sw(coev_plot_delta_theta(m09, limits = c(-5, 5))))
  expect_no_error(sw(coev_plot_delta_theta(m10, limits = c(-5, 5))))
  # prob and prob_outer work as expected
  expect_no_error(sw(coev_plot_delta_theta(m01, prob = 0.5, prob_outer = 0.89)))
  expect_no_error(sw(coev_plot_delta_theta(m02, prob = 0.5, prob_outer = 0.89)))
  expect_no_error(sw(coev_plot_delta_theta(m03, prob = 0.5, prob_outer = 0.89)))
  expect_no_error(sw(coev_plot_delta_theta(m04, prob = 0.5, prob_outer = 0.89)))
  expect_no_error(sw(coev_plot_delta_theta(m05, prob = 0.5, prob_outer = 0.89)))
  expect_no_error(sw(coev_plot_delta_theta(m06, prob = 0.5, prob_outer = 0.89)))
  expect_no_error(sw(coev_plot_delta_theta(m07, prob = 0.5, prob_outer = 0.89)))
  expect_no_error(sw(coev_plot_delta_theta(m08, prob = 0.5, prob_outer = 0.89)))
  expect_no_error(sw(coev_plot_delta_theta(m09, prob = 0.5, prob_outer = 0.89)))
  expect_no_error(sw(coev_plot_delta_theta(m10, prob = 0.5, prob_outer = 0.89)))
  # declaring variables works as expected
  expect_no_error(sw(coev_plot_delta_theta(m01, variables = c("x", "y"))))
})
