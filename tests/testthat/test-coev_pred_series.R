test_that("coev_pred_series() produces expected errors and output", {
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
    coev_pred_series(object = "fail"),
    "Argument 'object' must be a fitted coevolutionary model of class coevfit.",
    fixed = TRUE
  )
  expect_error(
    coev_pred_series(object = m02, eta_anc = c(1, 2)),
    "Argument 'eta_anc' is not a named list",
    fixed = TRUE
  )
  expect_error(
    coev_pred_series(object = m02, eta_anc = list(fail = 0)),
    "At least one variable in 'eta_anc' is not included in the fitted model.",
    fixed = TRUE
  )
  expect_error(
    coev_pred_series(object = m02, eta_anc = list(x = 0)),
    "All coevolving variables must be included in argument 'eta_anc'.",
    fixed = TRUE
  )
  expect_error(
    coev_pred_series(
      object = m01,
      eta_anc = list(x = 0, y = 0, u = 0, v = 0, w = 0, w = 0) # repetition
    ),
    "Argument 'eta_anc' contains duplicated variable names.",
    fixed = TRUE
  )
  expect_error(
    coev_pred_series(
      object = m01,
      eta_anc = list(x = 0, y = 0, u = 0, v = 0, w = c(0, 0)) # not length 1
    ),
    "Values in 'eta_anc' must each be of length one.",
    fixed = TRUE
  )
  expect_error(
    coev_pred_series(
      object = m01,
      eta_anc = list(x = "LCA", y = 0, u = 0, v = 0, w = 0)
    ),
    paste0(
      "Values in 'eta_anc' must each be numeric."
    ),
    fixed = TRUE
  )
  expect_error(
    coev_pred_series(object = m01, tmax = "fail"),
    "Argument 'tmax' must be a single numeric value",
    fixed = TRUE
  )
  expect_error(
    coev_pred_series(object = m01, tmax = -1),
    "Argument 'tmax' must be positive.",
    fixed = TRUE
  )
  expect_error(
    coev_pred_series(object = m01, ndraws = "fail"),
    "Argument 'ndraws' must be numeric.",
    fixed = TRUE
  )
  expect_error(
    coev_pred_series(object = m01, ndraws = c(0, 0)),
    "Argument 'ndraws' must be a single integer",
    fixed = TRUE
  )
  expect_error(
    coev_pred_series(object = m01, ndraws = 0L),
    "Argument 'ndraws' must be between 1 and the total number of draws.",
    fixed = TRUE
  )
  expect_error(
    coev_pred_series(
      object = m01, ndraws = as.integer(nrow(m01$fit$draws()) + 1)
    ),
    "Argument 'ndraws' must be between 1 and the total number of draws.",
    fixed = TRUE
  )
  expect_error(
    coev_pred_series(object = m01, stochastic = "only"),
    "Argument 'stochastic' must be logical.",
    fixed = TRUE
  )
  expect_error(
    coev_pred_series(
      object = m01,
      intervention_values = list(u = NA, v = 0, w = 0, x = 0, y = 0),
      stochastic = TRUE
    ),
    "Argument 'stochastic' cannot be `TRUE` when intervention_values are set.",
    fixed = TRUE
  )
  expect_error(
    coev_pred_series(object = m01, intervention_values = "fail"),
    "Argument 'intervention_values' is not a named list.",
    fixed = TRUE
  )
  expect_error(
    coev_pred_series(object = m01, intervention_values = list(fail = 0)),
    paste0(
      "At least one variable in 'intervention_values' is not ",
      "included in the fitted model."
    ),
    fixed = TRUE
  )
  expect_error(
    coev_pred_series(
      object = m01,
      intervention_values = list(x = 0)
    ),
    paste0(
      "All coevolving variables must be included in ",
      "argument 'intervention_values'."
    ),
    fixed = TRUE
  )
  expect_error(
    coev_pred_series(
      object = m01,
      intervention_values = list(x = 0, y = 0, u = 0, v = 0, w = 0, w = 0)
    ),
    "Argument 'intervention_values' contains duplicated variable names.",
    fixed = TRUE
  )
  expect_error(
    coev_pred_series(
      object = m01,
      intervention_values = list(x = 0, y = 0, u = 0, v = 0, w = c(0, 0))
    ),
    "Values in 'intervention_values' must each be of length one.",
    fixed = TRUE
  )
  expect_error(
    coev_pred_series(
      object = m01,
      intervention_values = list(x = 0, y = 0, u = 0, v = 0, w = "fail")
    ),
    "Values in 'intervention_values' must each be NA or numeric.",
    fixed = TRUE
  )
  expect_error(
    coev_pred_series(
      object = m01,
      intervention_values = list(x = 0, y = 0, u = 0, v = 0, w = 0)
    ),
    paste0(
      "Argument 'intervention_values' must have at least one NA value ",
      "declaring a free variable. If all variables are held constant, the ",
      "system is already at equilibrium and there is nothing to compute."
    ),
    fixed = TRUE
  )
  expect_error(
    coev_pred_series(
      object = m01,
      intervention_values = list(x = NA, y = NA, u = NA, v = NA, w = NA)
    ),
    paste0(
      "Argument 'intervention_values' must have at least one variable ",
      "held constant (i.e., not all values are NA)."
    ),
    fixed = TRUE
  )
  # suppress warnings
  sw <- suppressWarnings
  # should run without error
  expect_no_error(sw(coev_pred_series(m01)))
  expect_no_error(sw(coev_pred_series(m02)))
  expect_no_error(sw(coev_pred_series(m03)))
  expect_no_error(sw(coev_pred_series(m04)))
  expect_no_error(sw(coev_pred_series(m05)))
  expect_no_error(sw(coev_pred_series(m06)))
  expect_no_error(sw(coev_pred_series(m07)))
  expect_no_error(sw(coev_pred_series(m08)))
  expect_no_error(sw(coev_pred_series(m09)))
  expect_no_error(sw(coev_pred_series(m10)))
  expect_no_error(sw(coev_pred_series(m01, stochastic = TRUE)))
  expect_no_error(sw(coev_pred_series(m02, stochastic = TRUE)))
  expect_no_error(sw(coev_pred_series(m03, stochastic = TRUE)))
  expect_no_error(sw(coev_pred_series(m04, stochastic = TRUE)))
  expect_no_error(sw(coev_pred_series(m05, stochastic = TRUE)))
  expect_no_error(sw(coev_pred_series(m06, stochastic = TRUE)))
  expect_no_error(sw(coev_pred_series(m07, stochastic = TRUE)))
  expect_no_error(sw(coev_pred_series(m08, stochastic = TRUE)))
  expect_no_error(sw(coev_pred_series(m09, stochastic = TRUE)))
  expect_no_error(sw(coev_pred_series(m10, stochastic = TRUE)))
  expect_no_error(sw(coev_pred_series(m01, ndraws = 1)))
  expect_no_error(sw(coev_pred_series(m02, ndraws = 1)))
  expect_no_error(sw(coev_pred_series(m03, ndraws = 1)))
  expect_no_error(sw(coev_pred_series(m04, ndraws = 1)))
  expect_no_error(sw(coev_pred_series(m05, ndraws = 1)))
  expect_no_error(sw(coev_pred_series(m06, ndraws = 1)))
  expect_no_error(sw(coev_pred_series(m07, ndraws = 1)))
  expect_no_error(sw(coev_pred_series(m08, ndraws = 1)))
  expect_no_error(sw(coev_pred_series(m09, ndraws = 1)))
  expect_no_error(sw(coev_pred_series(m10, ndraws = 1)))
  expect_no_error(
    sw(
      coev_pred_series(
        m01, intervention_values = list(u = NA, v = 0, w = 0, x = 0, y = 0)
      )
    )
  )
  expect_no_error(
    sw(coev_pred_series(m02, intervention_values = list(w = NA, x = 0)))
  )
  expect_no_error(
    sw(coev_pred_series(m03, intervention_values = list(w = NA, x = 0)))
  )
  expect_no_error(
    sw(coev_pred_series(m04, intervention_values = list(y = NA, z = 0)))
  )
  expect_no_error(
    sw(coev_pred_series(m05, intervention_values = list(w = NA, x = 0)))
  )
  expect_no_error(
    sw(coev_pred_series(m06, intervention_values = list(w = NA, x = 0)))
  )
  expect_no_error(
    sw(coev_pred_series(m07, intervention_values = list(w = NA, x = 0)))
  )
  expect_no_error(
    sw(coev_pred_series(m08, intervention_values = list(x = NA, y = 0)))
  )
  expect_no_error(
    sw(coev_pred_series(m09, intervention_values = list(x = NA, y = 0)))
  )
  expect_no_error(
    sw(coev_pred_series(m10, intervention_values = list(x = NA, y = 0)))
  )
  expect_no_error(
    sw(coev_pred_series(m01, eta_anc = list(u = 0, v = 0, w = 0, x = 0, y = 0)))
  )
  expect_no_error(sw(coev_pred_series(m02, eta_anc = list(w = 0, x = 0))))
  expect_no_error(sw(coev_pred_series(m03, eta_anc = list(w = 0, x = 0))))
  expect_no_error(sw(coev_pred_series(m04, eta_anc = list(y = 0, z = 0))))
  expect_no_error(sw(coev_pred_series(m05, eta_anc = list(w = 0, x = 0))))
  expect_no_error(sw(coev_pred_series(m06, eta_anc = list(w = 0, x = 0))))
  expect_no_error(sw(coev_pred_series(m07, eta_anc = list(w = 0, x = 0))))
  expect_no_error(sw(coev_pred_series(m08, eta_anc = list(x = 0, y = 0))))
  expect_no_error(sw(coev_pred_series(m09, eta_anc = list(x = 0, y = 0))))
  expect_no_error(sw(coev_pred_series(m10, eta_anc = list(x = 0, y = 0))))
  # test for no NA, NaN, or Inf
  #' @srrstats {G5.3} Ensure no missing or undefined values
  expect_equal(sum(is.na(sw(coev_pred_series(m01)))), 0)
  expect_equal(sum(is.nan(sw(coev_pred_series(m01)))), 0)
  expect_equal(sum(is.infinite(sw(coev_pred_series(m01)))), 0)
})
