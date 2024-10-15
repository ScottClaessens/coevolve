test_that("coev_plot_delta_theta() produces expected errors and output", {
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
    coev_plot_delta_theta(object = "fail"),
    "Argument 'object' must be a fitted coevolutionary model of class coevfit.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_delta_theta(object = m1, variables = NA),
    "Argument 'variables' must be a character vector.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_delta_theta(object = m1, variables = "fail"),
    "Argument 'variables' must be of length > 1.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_delta_theta(object = m1, variables = c("x", "y", "fail")),
    "Some variables in 'variables' are not included in the fitted model.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_delta_theta(object = m1, variables = c("x", "y", "y")),
    "Argument 'variables' contains duplicates.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_delta_theta(object = m1, prob = "fail"),
    "Argument 'prob' must be numeric.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_delta_theta(object = m1, prob = c(0.5, 0.5)),
    "Argument 'prob' must be of length 1.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_delta_theta(object = m1, prob = -0.5),
    "Argument 'prob' must be between 0 and 1.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_delta_theta(object = m1, prob = 1.5),
    "Argument 'prob' must be between 0 and 1.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_delta_theta(object = m1, prob_outer = "fail"),
    "Argument 'prob_outer' must be numeric.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_delta_theta(object = m1, prob_outer = c(0.5, 0.5)),
    "Argument 'prob_outer' must be of length 1.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_delta_theta(object = m1, prob_outer = -0.5),
    "Argument 'prob_outer' must be between 0 and 1.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_delta_theta(object = m1, prob_outer = 1.5),
    "Argument 'prob_outer' must be between 0 and 1.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_delta_theta(object = m1, prob = 0.5, prob_outer = 0.4),
    "Argument 'prob_outer' must be greater than argument 'prob'.",
    fixed = TRUE
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
  expect_no_error(SW(coev_plot_delta_theta(m7)))
  expect_no_error(SW(coev_plot_delta_theta(m8)))
  expect_no_error(SW(coev_plot_delta_theta(m9)))
  expect_true(methods::is(SW(coev_plot_delta_theta(m1)), "ggplot"))
  expect_true(methods::is(SW(coev_plot_delta_theta(m2)), "ggplot"))
  expect_true(methods::is(SW(coev_plot_delta_theta(m3)), "ggplot"))
  expect_true(methods::is(SW(coev_plot_delta_theta(m4)), "ggplot"))
  expect_true(methods::is(SW(coev_plot_delta_theta(m5)), "ggplot"))
  expect_true(methods::is(SW(coev_plot_delta_theta(m6)), "ggplot"))
  expect_true(methods::is(SW(coev_plot_delta_theta(m7)), "ggplot"))
  expect_true(methods::is(SW(coev_plot_delta_theta(m8)), "ggplot"))
  expect_true(methods::is(SW(coev_plot_delta_theta(m9)), "ggplot"))
  # limits work as expected
  expect_no_error(SW(coev_plot_delta_theta(m1, limits = c(-5, 5))))
  expect_no_error(SW(coev_plot_delta_theta(m2, limits = c(-5, 5))))
  expect_no_error(SW(coev_plot_delta_theta(m3, limits = c(-5, 5))))
  expect_no_error(SW(coev_plot_delta_theta(m4, limits = c(-5, 5))))
  expect_no_error(SW(coev_plot_delta_theta(m5, limits = c(-5, 5))))
  expect_no_error(SW(coev_plot_delta_theta(m6, limits = c(-5, 5))))
  expect_no_error(SW(coev_plot_delta_theta(m7, limits = c(-5, 5))))
  expect_no_error(SW(coev_plot_delta_theta(m8, limits = c(-5, 5))))
  expect_no_error(SW(coev_plot_delta_theta(m9, limits = c(-5, 5))))
  # prob and prob_outer work as expected
  expect_no_error(SW(coev_plot_delta_theta(m1, prob = 0.5, prob_outer = 0.89)))
  expect_no_error(SW(coev_plot_delta_theta(m2, prob = 0.5, prob_outer = 0.89)))
  expect_no_error(SW(coev_plot_delta_theta(m3, prob = 0.5, prob_outer = 0.89)))
  expect_no_error(SW(coev_plot_delta_theta(m4, prob = 0.5, prob_outer = 0.89)))
  expect_no_error(SW(coev_plot_delta_theta(m5, prob = 0.5, prob_outer = 0.89)))
  expect_no_error(SW(coev_plot_delta_theta(m6, prob = 0.5, prob_outer = 0.89)))
  expect_no_error(SW(coev_plot_delta_theta(m7, prob = 0.5, prob_outer = 0.89)))
  expect_no_error(SW(coev_plot_delta_theta(m8, prob = 0.5, prob_outer = 0.89)))
  expect_no_error(SW(coev_plot_delta_theta(m9, prob = 0.5, prob_outer = 0.89)))
  # declaring variables works as expected
  expect_no_error(SW(coev_plot_delta_theta(m1, variables = c("x", "y"))))
})
