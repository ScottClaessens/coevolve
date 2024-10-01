test_that("coev_plot_predictive_check() produces expected errors and output", {
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
    coev_plot_predictive_check(object = "fail"),
    "Argument 'object' must be a fitted coevolutionary model of class coevfit."
  )
  expect_error(
    coev_plot_predictive_check(object = m1, variables = 0),
    paste0(
      "Argument 'variables' must be a character string or a ",
      "vector of character strings."
    )
  )
  expect_error(
    coev_plot_predictive_check(object = m1, variables = "fail"),
    paste0(
      "Argument 'variables' contains variable names that are not ",
      "included in the fitted model."
    )
  )
  expect_error(
    coev_plot_predictive_check(object = m1, ndraws = "fail"),
    "Argument 'ndraws' must be a single integer."
  )
  expect_error(
    coev_plot_predictive_check(object = m1, ndraws = 0L),
    "Argument 'ndraws' must be between 1 and the total number of draws."
  )
  expect_error(
    coev_plot_predictive_check(
      object = m1, ndraws = as.integer(nrow(m1$fit$draws()) + 1)
      ),
    "Argument 'ndraws' must be between 1 and the total number of draws."
  )
  # suppress warnings
  SW <- suppressWarnings
  # should run without error and produce list of ggplot objects
  expect_no_error(SW(coev_plot_predictive_check(m1)))
  expect_no_error(SW(coev_plot_predictive_check(m2)))
  expect_no_error(SW(coev_plot_predictive_check(m3)))
  expect_no_error(SW(coev_plot_predictive_check(m4)))
  expect_no_error(SW(coev_plot_predictive_check(m5)))
  expect_no_error(SW(coev_plot_predictive_check(m6)))
  expect_no_error(SW(coev_plot_predictive_check(m7)))
  expect_no_error(SW(coev_plot_predictive_check(m8)))
  expect_no_error(SW(coev_plot_predictive_check(m9)))
  expect_no_error(SW(coev_plot_predictive_check(m1, variables = "x")))
  expect_no_error(SW(coev_plot_predictive_check(m2, variables = "w")))
  expect_no_error(SW(coev_plot_predictive_check(m3, variables = "w")))
  expect_no_error(SW(coev_plot_predictive_check(m4, variables = "y")))
  expect_no_error(SW(coev_plot_predictive_check(m5, variables = "w")))
  expect_no_error(SW(coev_plot_predictive_check(m6, variables = "w")))
  expect_no_error(SW(coev_plot_predictive_check(m7, variables = "w")))
  expect_no_error(SW(coev_plot_predictive_check(m8, variables = "x")))
  expect_no_error(SW(coev_plot_predictive_check(m9, variables = "x")))
  expect_no_error(SW(coev_plot_predictive_check(m1, ndraws = 1L)))
  expect_no_error(SW(coev_plot_predictive_check(m2, ndraws = 1L)))
  expect_no_error(SW(coev_plot_predictive_check(m3, ndraws = 1L)))
  expect_no_error(SW(coev_plot_predictive_check(m4, ndraws = 1L)))
  expect_no_error(SW(coev_plot_predictive_check(m5, ndraws = 1L)))
  expect_no_error(SW(coev_plot_predictive_check(m6, ndraws = 1L)))
  expect_no_error(SW(coev_plot_predictive_check(m7, ndraws = 1L)))
  expect_no_error(SW(coev_plot_predictive_check(m8, ndraws = 1L)))
  expect_no_error(SW(coev_plot_predictive_check(m9, ndraws = 1L)))
  # should work with multiPhylo
  expect_no_error(SW(coev_plot_predictive_check(m8, tree_id = 2L)))
  # should work as expected with missing data
  # lower limit of plot should not be -9999.45
  expect_false(
    isTRUE(
      all.equal(
        ggplot2::layer_scales(
          SW(coev_plot_predictive_check(m6, variables = "x"))[[1]]
          )$x$range$range[1],
        -9999.45
        )
      )
    )
})
