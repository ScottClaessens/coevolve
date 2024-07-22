test_that("coev_calculate_delta_theta() produces expected errors and output", {
  # load model
  m <- readRDS(test_path("fixtures", "coevfit_example1.rds"))
  m <- reload_fit(m, filename = "coevfit_example1-1.csv")
  # expect the following errors
  expect_error(
    coev_calculate_delta_theta(object = "fail"),
    "Argument 'object' must be a fitted coevolutionary model of class coevfit."
  )
  expect_error(
    coev_calculate_delta_theta(object = m, response = c(NA, "fail")),
    "Argument 'response' must be a character string of length one."
  )
  expect_error(
    coev_calculate_delta_theta(object = m, response = "fail"),
    "Argument 'response' must be a variable included in the fitted model."
  )
  expect_error(
    coev_calculate_delta_theta(object = m, response = "x",
                               predictor = c(NA, "fail")),
    "Argument 'predictor' must be a character string of length one."
  )
  expect_error(
    coev_calculate_delta_theta(object = m, response = "x", predictor = "fail"),
    "Argument 'predictor' must be a variable included in the fitted model."
  )
  expect_error(
    coev_calculate_delta_theta(object = m, response = "x", predictor = "x"),
    "Argument 'response' and 'predictor' must refer to different variables."
  )
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
