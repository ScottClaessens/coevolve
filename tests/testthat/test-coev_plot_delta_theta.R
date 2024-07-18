test_that("coev_plot_delta_theta() produces expected errors and output", {
  # load model
  m <- coevolve:::coevfit_example1
  # expect the following errors
  expect_error(
    coev_plot_delta_theta(object = "fail"),
    "Argument 'object' must be a fitted coevolutionary model of class coevfit."
  )
  expect_error(
    coev_plot_delta_theta(object = m, variables = NA),
    "Argument 'variables' must be a character vector."
  )
  expect_error(
    coev_plot_delta_theta(object = m, variables = "fail"),
    "Argument 'variables' must be of length > 1."
  )
  expect_error(
    coev_plot_delta_theta(object = m, variables = c("x", "y", "fail")),
    "Some variables in 'variables' are not included in the fitted model."
  )
  expect_error(
    coev_plot_delta_theta(object = m, variables = c("x", "y", "y")),
    "Argument 'variables' contains duplicates."
  )
  # suppress warnings
  SW <- suppressWarnings
  # should run without error and produce ggplot object
  expect_no_error(SW(coev_plot_delta_theta(m)))
  expect_true(methods::is(SW(coev_plot_delta_theta(m)), "ggplot"))
})
