test_that("coev_calculate_theta() produces expected errors and output", {
  # load model
  m1 <- readRDS(test_path("fixtures", "coevfit_example1.rds"))
  m2 <- readRDS(test_path("fixtures", "coevfit_example2.rds"))
  m3 <- readRDS(test_path("fixtures", "coevfit_example3.rds"))
  m4 <- readRDS(test_path("fixtures", "coevfit_example4.rds"))
  m5 <- readRDS(test_path("fixtures", "coevfit_example5.rds"))
  m6 <- readRDS(test_path("fixtures", "coevfit_example6.rds"))
  m7 <- readRDS(test_path("fixtures", "coevfit_example7.rds"))
  m8 <- readRDS(test_path("fixtures", "coevfit_example8.rds"))
  m1 <- reload_fit(m1, filename = "coevfit_example1-1.csv")
  m2 <- reload_fit(m2, filename = "coevfit_example2-1.csv")
  m3 <- reload_fit(m3, filename = "coevfit_example3-1.csv")
  m4 <- reload_fit(m4, filename = "coevfit_example4-1.csv")
  m5 <- reload_fit(m5, filename = "coevfit_example5-1.csv")
  m6 <- reload_fit(m6, filename = "coevfit_example6-1.csv")
  m7 <- reload_fit(m7, filename = "coevfit_example7-1.csv")
  m8 <- reload_fit(m8, filename = "coevfit_example8-1.csv")
  # expect the following errors
  expect_error(
    coev_calculate_theta(object = "fail"),
    "Argument 'object' must be a fitted coevolutionary model of class coevfit."
  )
  expect_error(
    coev_calculate_theta(object = m4, intervention_values = list("fail")),
    "Argument 'intervention_values' is not a named list."
  )
  expect_error(
    coev_calculate_theta(object = m4, intervention_values = list(fail = NA)),
    paste0(
      "At least one variable in 'intervention_values' is not included in ",
      "the fitted model."
    )
  )
  expect_error(
    coev_calculate_theta(object = m4, intervention_values = list(y = NA)),
    paste0(
      "All coevolving variables must be included in ",
      "argument 'intervention_values'."
    )
  )
  expect_error(
    coev_calculate_theta(
      object = m4,
      intervention_values = list(y = NA, y = NA, z = NA)
    ),
    "Argument 'intervention_values' contains duplicated variable names."
  )
  expect_error(
    coev_calculate_theta(
      object = m4,
      intervention_values = list(y = c(NA, NA), z = 0)
    ),
    "Values in 'intervention_values' must each be of length one."
  )
  expect_error(
    coev_calculate_theta(
      object = m4,
      intervention_values = list(y = "fail", z = 0)
    ),
    "Values in 'intervention_values' must each be NA or numeric."
  )
  expect_error(
    coev_calculate_theta(
      object = m4,
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
      object = m4,
      intervention_values = list(y = NA, z = NA)
    ),
    paste0(
      "Argument 'intervention_values' must have at least one variable ",
      "held constant (i.e., all values are NA)."
    ),
    fixed = TRUE
  )
  # should run without error
  theta1 <- coev_calculate_theta(m1, list(u = NA, v = 0, w = 0, x = 0, y = 0))
  theta2 <- coev_calculate_theta(m2, list(w = NA, x = 0))
  theta3 <- coev_calculate_theta(m3, list(w = NA, x = 0))
  theta4 <- coev_calculate_theta(m4, list(y = NA, z = 0))
  theta5 <- coev_calculate_theta(m5, list(w = NA, x = 0))
  theta6 <- coev_calculate_theta(m6, list(w = NA, x = 0))
  theta7 <- coev_calculate_theta(m7, list(w = NA, x = 0))
  theta8 <- coev_calculate_theta(m8, list(x = NA, y = 0))
  expect_no_error(theta1)
  expect_no_error(theta2)
  expect_no_error(theta3)
  expect_no_error(theta4)
  expect_no_error(theta5)
  expect_no_error(theta6)
  expect_no_error(theta7)
  expect_no_error(theta8)
  # output should be matrix of posterior draws
  expect_true(methods::is(theta1, "matrix"))
  expect_true(methods::is(theta2, "matrix"))
  expect_true(methods::is(theta3, "matrix"))
  expect_true(methods::is(theta4, "matrix"))
  expect_true(methods::is(theta5, "matrix"))
  expect_true(methods::is(theta6, "matrix"))
  expect_true(methods::is(theta7, "matrix"))
  expect_true(methods::is(theta8, "matrix"))
  # output column names should be equal to variable names
  expect_true(identical(colnames(theta1), names(m1$variables)))
  expect_true(identical(colnames(theta2), names(m2$variables)))
  expect_true(identical(colnames(theta3), names(m3$variables)))
  expect_true(identical(colnames(theta4), names(m4$variables)))
  expect_true(identical(colnames(theta5), names(m5$variables)))
  expect_true(identical(colnames(theta6), names(m6$variables)))
  expect_true(identical(colnames(theta7), names(m7$variables)))
  expect_true(identical(colnames(theta8), names(m8$variables)))
  # variables should be correctly held constant in output
  expect_true(all(theta1[,c("v","w","x","y")] == 0))
  expect_true(all(theta2[,"x"] == 0))
  expect_true(all(theta3[,"x"] == 0))
  expect_true(all(theta4[,"z"] == 0))
  expect_true(all(theta5[,"x"] == 0))
  expect_true(all(theta6[,"x"] == 0))
  expect_true(all(theta7[,"x"] == 0))
  expect_true(all(theta8[,"y"] == 0))
})
