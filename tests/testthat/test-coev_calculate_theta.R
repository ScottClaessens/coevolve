test_that("coev_calculate_theta() produces expected errors and output", {
  # load model
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
    coev_calculate_theta(object = "fail"),
    "Argument 'object' must be a fitted coevolutionary model of class coevfit."
  )
  expect_error(
    coev_calculate_theta(object = m04, intervention_values = list("fail")),
    "Argument 'intervention_values' is not a named list."
  )
  expect_error(
    coev_calculate_theta(object = m04, intervention_values = list(fail = NA)),
    paste0(
      "At least one variable in 'intervention_values' is not included in ",
      "the fitted model."
    )
  )
  expect_error(
    coev_calculate_theta(object = m04, intervention_values = list(y = NA)),
    paste0(
      "All coevolving variables must be included in ",
      "argument 'intervention_values'."
    )
  )
  expect_error(
    coev_calculate_theta(
      object = m04,
      intervention_values = list(y = NA, y = NA, z = NA)
    ),
    "Argument 'intervention_values' contains duplicated variable names."
  )
  expect_error(
    coev_calculate_theta(
      object = m04,
      intervention_values = list(y = c(NA, NA), z = 0)
    ),
    "Values in 'intervention_values' must each be of length one."
  )
  expect_error(
    coev_calculate_theta(
      object = m04,
      intervention_values = list(y = "fail", z = 0)
    ),
    "Values in 'intervention_values' must each be NA or numeric."
  )
  expect_error(
    coev_calculate_theta(
      object = m04,
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
      object = m04,
      intervention_values = list(y = NA, z = NA)
    ),
    paste0(
      "Argument 'intervention_values' must have at least one variable ",
      "held constant (i.e., all values are NA)."
    ),
    fixed = TRUE
  )
  # should run without error
  theta01 <- coev_calculate_theta(m01, list(u = NA, v = 0, w = 0, x = 0, y = 0))
  theta02 <- coev_calculate_theta(m02, list(w = NA, x = 0))
  theta03 <- coev_calculate_theta(m03, list(w = NA, x = 0))
  theta04 <- coev_calculate_theta(m04, list(y = NA, z = 0))
  theta05 <- coev_calculate_theta(m05, list(w = NA, x = 0))
  theta06 <- coev_calculate_theta(m06, list(w = NA, x = 0))
  theta07 <- coev_calculate_theta(m07, list(w = NA, x = 0))
  theta08 <- coev_calculate_theta(m08, list(x = NA, y = 0))
  theta09 <- coev_calculate_theta(m09, list(x = NA, y = 0))
  theta10 <- coev_calculate_theta(m10, list(x = NA, y = 0))
  theta01_null <- coev_calculate_theta(m01, intervention_values = NULL)
  theta02_null <- coev_calculate_theta(m02, intervention_values = NULL)
  theta03_null <- coev_calculate_theta(m03, intervention_values = NULL)
  theta04_null <- coev_calculate_theta(m04, intervention_values = NULL)
  theta05_null <- coev_calculate_theta(m05, intervention_values = NULL)
  theta06_null <- coev_calculate_theta(m06, intervention_values = NULL)
  theta07_null <- coev_calculate_theta(m07, intervention_values = NULL)
  theta08_null <- coev_calculate_theta(m08, intervention_values = NULL)
  theta09_null <- coev_calculate_theta(m09, intervention_values = NULL)
  theta10_null <- coev_calculate_theta(m10, intervention_values = NULL)
  expect_no_error(theta01)
  expect_no_error(theta02)
  expect_no_error(theta03)
  expect_no_error(theta04)
  expect_no_error(theta05)
  expect_no_error(theta06)
  expect_no_error(theta07)
  expect_no_error(theta08)
  expect_no_error(theta09)
  expect_no_error(theta10)
  expect_no_error(theta01_null)
  expect_no_error(theta02_null)
  expect_no_error(theta03_null)
  expect_no_error(theta04_null)
  expect_no_error(theta05_null)
  expect_no_error(theta06_null)
  expect_no_error(theta07_null)
  expect_no_error(theta08_null)
  expect_no_error(theta09_null)
  expect_no_error(theta10_null)
  # output should be matrix of posterior draws
  expect_true(methods::is(theta01, "matrix"))
  expect_true(methods::is(theta02, "matrix"))
  expect_true(methods::is(theta03, "matrix"))
  expect_true(methods::is(theta04, "matrix"))
  expect_true(methods::is(theta05, "matrix"))
  expect_true(methods::is(theta06, "matrix"))
  expect_true(methods::is(theta07, "matrix"))
  expect_true(methods::is(theta08, "matrix"))
  expect_true(methods::is(theta09, "matrix"))
  expect_true(methods::is(theta10, "matrix"))
  expect_true(methods::is(theta01_null, "matrix"))
  expect_true(methods::is(theta02_null, "matrix"))
  expect_true(methods::is(theta03_null, "matrix"))
  expect_true(methods::is(theta04_null, "matrix"))
  expect_true(methods::is(theta05_null, "matrix"))
  expect_true(methods::is(theta06_null, "matrix"))
  expect_true(methods::is(theta07_null, "matrix"))
  expect_true(methods::is(theta08_null, "matrix"))
  expect_true(methods::is(theta09_null, "matrix"))
  expect_true(methods::is(theta10_null, "matrix"))
  # output column names should be equal to variable names
  expect_true(identical(colnames(theta01), names(m01$variables)))
  expect_true(identical(colnames(theta02), names(m02$variables)))
  expect_true(identical(colnames(theta03), names(m03$variables)))
  expect_true(identical(colnames(theta04), names(m04$variables)))
  expect_true(identical(colnames(theta05), names(m05$variables)))
  expect_true(identical(colnames(theta06), names(m06$variables)))
  expect_true(identical(colnames(theta07), names(m07$variables)))
  expect_true(identical(colnames(theta08), names(m08$variables)))
  expect_true(identical(colnames(theta09), names(m09$variables)))
  expect_true(identical(colnames(theta10), names(m10$variables)))
  expect_true(identical(colnames(theta01_null), names(m01$variables)))
  expect_true(identical(colnames(theta02_null), names(m02$variables)))
  expect_true(identical(colnames(theta03_null), names(m03$variables)))
  expect_true(identical(colnames(theta04_null), names(m04$variables)))
  expect_true(identical(colnames(theta05_null), names(m05$variables)))
  expect_true(identical(colnames(theta06_null), names(m06$variables)))
  expect_true(identical(colnames(theta07_null), names(m07$variables)))
  expect_true(identical(colnames(theta08_null), names(m08$variables)))
  expect_true(identical(colnames(theta09_null), names(m09$variables)))
  expect_true(identical(colnames(theta10_null), names(m10$variables)))
  # variables should be correctly held constant in output
  expect_true(all(theta01[, c("v", "w", "x", "y")] == 0))
  expect_true(all(theta02[, "x"] == 0))
  expect_true(all(theta03[, "x"] == 0))
  expect_true(all(theta04[, "z"] == 0))
  expect_true(all(theta05[, "x"] == 0))
  expect_true(all(theta06[, "x"] == 0))
  expect_true(all(theta07[, "x"] == 0))
  expect_true(all(theta08[, "y"] == 0))
  expect_true(all(theta09[, "y"] == 0))
  expect_true(all(theta10[, "y"] == 0))
  # free variables should not be constant in output
  expect_false(all(theta01_null[, "y"] == 0))
  expect_false(all(theta09_null[, "y"] == 0))
  # test for no NA, NaN, or Inf
  #' @srrstats {G5.3} Ensure no missing or undefined values
  expect_equal(sum(is.na(theta01)), 0)
  expect_equal(sum(is.nan(theta01)), 0)
  expect_equal(sum(is.infinite(theta01)), 0)
})
