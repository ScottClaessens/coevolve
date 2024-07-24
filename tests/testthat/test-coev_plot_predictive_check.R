test_that("coev_plot_predictive_check() produces expected errors and output", {
  # load model
  m <- readRDS(test_path("fixtures", "coevfit_example1.rds"))
  m <- reload_fit(m, filename = "coevfit_example1-1.csv")
  # expect the following errors
  expect_error(
    coev_plot_predictive_check(object = "fail"),
    "Argument 'object' must be a fitted coevolutionary model of class coevfit."
  )
  expect_error(
    coev_plot_predictive_check(object = m, variables = 0),
    paste0(
      "Argument 'variables' must be a character string or a ",
      "vector of character strings."
    )
  )
  expect_error(
    coev_plot_predictive_check(object = m, variables = "fail"),
    paste0(
      "Argument 'variables' contains variable names that are not ",
      "included in the fitted model."
    )
  )
  expect_error(
    coev_plot_predictive_check(object = m, ndraws = "fail"),
    "Argument 'ndraws' must be a single integer."
  )
  expect_error(
    coev_plot_predictive_check(object = m, ndraws = 0L),
    "Argument 'ndraws' must be between 1 and the total number of draws."
  )
  expect_error(
    coev_plot_predictive_check(
      object = m, ndraws = as.integer(nrow(m$fit$draws()) + 1)
      ),
    "Argument 'ndraws' must be between 1 and the total number of draws."
  )
  # suppress warnings
  SW <- suppressWarnings
  # should run without error and produce list of ggplot objects
  expect_no_error(SW(coev_plot_predictive_check(m)))
  expect_no_error(SW(coev_plot_predictive_check(m, variables = c("u","v"))))
  expect_no_error(SW(coev_plot_predictive_check(m, ndraws = 1L)))
})
