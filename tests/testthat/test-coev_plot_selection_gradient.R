test_that("coev_plot_selection_gradient() produces expected errors and output", {
  # load model
  m <- coevolve:::coevfit_example1
  # expect the following errors
  expect_error(
    coev_plot_selection_gradient(object = "fail", var1 = "x", var2 = "y"),
    "Argument 'object' must be a fitted coevolutionary model of class coevfit."
  )
  expect_error(
    coev_plot_selection_gradient(object = m, var1 = 0:1, var2 = "y"),
    "Argument 'var1' must be a character string of length one."
  )
  expect_error(
    coev_plot_selection_gradient(object = m, var1 = "z", var2 = "y"),
    "Argument 'var1' must be a variable included in the fitted model."
  )
  expect_error(
    coev_plot_selection_gradient(object = m, var1 = "x", var2 = 0:1),
    "Argument 'var2' must be a character string of length one."
  )
  expect_error(
    coev_plot_selection_gradient(object = m, var1 = "x", var2 = "z"),
    "Argument 'var2' must be a variable included in the fitted model."
  )
  expect_error(
    coev_plot_selection_gradient(object = m, var1 = "x", var2 = "x"),
    "Argument 'var1' and 'var2' must refer to different variables."
  )
  expect_error(
    coev_plot_selection_gradient(
      object = m, var1 = "x", var2 = "y", contour = "hello"
      ),
    "Argument 'contour' must be logical."
  )
  # suppress warnings
  SW <- suppressWarnings
  # should run without error and produce ggplot object
  expect_no_error(
    SW(coev_plot_selection_gradient(m, var1 = "x", var2 = "y"))
    )
  expect_no_error(
    SW(coev_plot_selection_gradient(m, var1 = "x", var2 = "y", contour = TRUE))
    )
  expect_true(
    methods::is(
      SW(coev_plot_selection_gradient(m, var1 = "x", var2 = "y")),
      "ggplot"
      )
    )
})
