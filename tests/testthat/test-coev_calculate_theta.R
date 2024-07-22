test_that("coev_calculate_theta() produces expected errors and output", {
  # load model
  m <- readRDS(test_path("fixtures", "coevfit_example4.rds"))
  m <- reload_fit(m, filename = "coevfit_example4-1.csv")
  # expect the following errors
  expect_error(
    coev_calculate_theta(object = "fail"),
    "Argument 'object' must be a fitted coevolutionary model of class coevfit."
  )
  expect_error(
    coev_calculate_theta(object = m, intervention_values = list("fail")),
    "Argument 'intervention_values' is not a named list."
  )
  expect_error(
    coev_calculate_theta(object = m, intervention_values = list(fail = NA)),
    paste0(
      "At least one variable in 'intervention_values' is not included in ",
      "the fitted model."
    )
  )
  expect_error(
    coev_calculate_theta(object = m, intervention_values = list(y = NA)),
    paste0(
      "At least one coevolving variable in the model is not included in ",
      "argument 'intervention_values'."
    )
  )
  expect_error(
    coev_calculate_theta(
      object = m,
      intervention_values = list(y = NA, y = NA, z = NA)
    ),
    "Argument 'intervention_values' contains duplicated variable names."
  )
  expect_error(
    coev_calculate_theta(
      object = m,
      intervention_values = list(y = c(NA, NA), z = 0)
    ),
    "Values in 'intervention_values' must each be of length one."
  )
  expect_error(
    coev_calculate_theta(
      object = m,
      intervention_values = list(y = "fail", z = 0)
    ),
    "Values in 'intervention_values' must each be NA or numeric."
  )
  expect_error(
    coev_calculate_theta(
      object = m,
      intervention_values = list(y = 0, z = 0)
    ),
    paste0(
      "Argument 'intervention_values' must have at least one NA value ",
      "declaring a free variable. If all variables are held constant, the ",
      "system is already at equilibrium and there is nothing to compute."
    )
  )
  expect_error(
    coev_calculate_theta(
      object = m,
      intervention_values = list(y = NA, z = NA)
    ),
    paste0(
      "Argument 'intervention_values' must have at least one variable ",
      "held constant (i.e., all values are NA)."
    ),
    fixed = TRUE
  )
  # should run without error
  theta1 <-
    coev_calculate_theta(
      object = m,
      intervention_values = list(y = 0, z = NA)
      )
  theta2 <-
    coev_calculate_theta(
      object = m,
      intervention_values = list(y = NA, z = 0)
    )
  expect_no_error(theta1)
  expect_no_error(theta2)
  # output should be matrix of posterior draws
  expect_true(methods::is(theta1, "matrix"))
  expect_true(methods::is(theta2, "matrix"))
  # output column names should be equal to variable names
  expect_true(identical(colnames(theta1), names(m$variables)))
  expect_true(identical(colnames(theta2), names(m$variables)))
  # variables should be correctly held constant in output
  expect_true(all(theta1[,"y"] == 0))
  expect_true(all(theta2[,"z"] == 0))
})
