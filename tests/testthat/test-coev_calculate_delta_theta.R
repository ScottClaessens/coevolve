test_that("coev_calculate_delta_theta() produces expected errors and output", {
  # load model
  m1 <- readRDS(test_path("fixtures", "coevfit_example1.rds"))
  m2 <- readRDS(test_path("fixtures", "coevfit_example2.rds"))
  m3 <- readRDS(test_path("fixtures", "coevfit_example3.rds"))
  m4 <- readRDS(test_path("fixtures", "coevfit_example4.rds"))
  m5 <- readRDS(test_path("fixtures", "coevfit_example5.rds"))
  m6 <- readRDS(test_path("fixtures", "coevfit_example6.rds"))
  m1 <- reload_fit(m1, filename = "coevfit_example1-1.csv")
  m2 <- reload_fit(m2, filename = "coevfit_example2-1.csv")
  m3 <- reload_fit(m3, filename = "coevfit_example3-1.csv")
  m4 <- reload_fit(m4, filename = "coevfit_example4-1.csv")
  m5 <- reload_fit(m5, filename = "coevfit_example5-1.csv")
  m6 <- reload_fit(m6, filename = "coevfit_example6-1.csv")
  # expect the following errors
  expect_error(
    coev_calculate_delta_theta(object = "fail"),
    "Argument 'object' must be a fitted coevolutionary model of class coevfit."
  )
  expect_error(
    coev_calculate_delta_theta(object = m1, response = c(NA, "fail")),
    "Argument 'response' must be a character string of length one."
  )
  expect_error(
    coev_calculate_delta_theta(object = m1, response = "fail"),
    "Argument 'response' must be a variable included in the fitted model."
  )
  expect_error(
    coev_calculate_delta_theta(object = m1, response = "x",
                               predictor = c(NA, "fail")),
    "Argument 'predictor' must be a character string of length one."
  )
  expect_error(
    coev_calculate_delta_theta(object = m1, response = "x", predictor = "fail"),
    "Argument 'predictor' must be a variable included in the fitted model."
  )
  expect_error(
    coev_calculate_delta_theta(object = m1, response = "x", predictor = "x"),
    "Argument 'response' and 'predictor' must refer to different variables."
  )
  # suppress warnings
  fun <- function(object, response, predictor) {
    suppressWarnings(coev_calculate_delta_theta(object, response, predictor))
  }
  # should run without error and produce draws_array objects
  expect_no_error(fun(m1, "x", "y"))
  expect_no_error(fun(m2, "w", "x"))
  expect_no_error(fun(m3, "w", "x"))
  expect_no_error(fun(m4, "y", "z"))
  expect_no_error(fun(m5, "w", "x"))
  expect_no_error(fun(m6, "w", "x"))
  expect_true(methods::is(fun(m1, "x", "y"), "draws_array"))
  expect_true(methods::is(fun(m2, "w", "x"), "draws_array"))
  expect_true(methods::is(fun(m3, "w", "x"), "draws_array"))
  expect_true(methods::is(fun(m4, "y", "z"), "draws_array"))
  expect_true(methods::is(fun(m5, "w", "x"), "draws_array"))
  expect_true(methods::is(fun(m6, "w", "x"), "draws_array"))
})

test_that("coev_calculate_delta_theta() works with repeated observations", {
  # load model
  m <- readRDS(test_path("fixtures", "coevfit_example7.rds"))
  m <- reload_fit(m, filename = "coevfit_example7-1.csv")
  # suppress warnings
  SW <- suppressWarnings
  # should run without error and produce draws_array object
  expect_no_error(
    SW(coev_calculate_delta_theta(m, response = "w", predictor = "x"))
  )
  expect_true(
    methods::is(
      SW(coev_calculate_delta_theta(m, response = "w", predictor = "x")),
      "draws_array"
    )
  )
})

test_that("coev_calculate_delta_theta() works with multiPhylo object", {
  # load model
  m <- readRDS(test_path("fixtures", "coevfit_example8.rds"))
  m <- reload_fit(m, filename = "coevfit_example8-1.csv")
  # suppress warnings
  SW <- suppressWarnings
  # should run without error and produce draws_array object
  expect_no_error(
    SW(coev_calculate_delta_theta(m, response = "x", predictor = "y"))
  )
  expect_true(
    methods::is(
      SW(coev_calculate_delta_theta(m, response = "x", predictor = "y")),
      "draws_array"
    )
  )
})

test_that("coev_calculate_delta_theta() works when Q offdiag == 0", {
  # load model
  m <- readRDS(test_path("fixtures", "coevfit_example9.rds"))
  m <- reload_fit(m, filename = "coevfit_example9-1.csv")
  # suppress warnings
  SW <- suppressWarnings
  # should run without error and produce draws_array object
  expect_no_error(
    SW(coev_calculate_delta_theta(m, response = "x", predictor = "y"))
  )
  expect_true(
    methods::is(
      SW(coev_calculate_delta_theta(m, response = "x", predictor = "y")),
      "draws_array"
    )
  )
})
