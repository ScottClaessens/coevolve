test_that("save_coevfit() produces expected errors and output", {
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
  expect_error(
    save_coevfit(object = "fail"),
    "Argument 'object' must be a fitted coevolutionary model of class coevfit.",
    fixed = TRUE
  )
  expect_error(
    save_coevfit(object = m01, file = 1),
    "Argument 'file' must be a string of length one.",
    fixed = TRUE
  )
  expect_error(
    save_coevfit(object = m01, file = c("fail", "fail")),
    "Argument 'file' must be a string of length one.",
    fixed = TRUE
  )
  # should run without error
  temp_file <- withr::local_tempfile(fileext = ".rds")
  expect_no_error(save_coevfit(object = m01, file = temp_file))
  expect_no_error(save_coevfit(object = m02, file = temp_file))
  expect_no_error(save_coevfit(object = m03, file = temp_file))
  expect_no_error(save_coevfit(object = m04, file = temp_file))
  expect_no_error(save_coevfit(object = m05, file = temp_file))
  expect_no_error(save_coevfit(object = m06, file = temp_file))
  expect_no_error(save_coevfit(object = m07, file = temp_file))
  expect_no_error(save_coevfit(object = m08, file = temp_file))
  expect_no_error(save_coevfit(object = m09, file = temp_file))
  expect_no_error(save_coevfit(object = m10, file = temp_file))
})
