test_that("coev_plot_flowfield() produces expected errors and output", {
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
    coev_plot_flowfield(object = "fail", var1 = "x", var2 = "y"),
    "Argument 'object' must be a fitted coevolutionary model of class coevfit."
  )
  expect_error(
    coev_plot_flowfield(object = m01, var1 = 0:1, var2 = "y"),
    "Argument 'var1' must be a character string of length one."
  )
  expect_error(
    coev_plot_flowfield(object = m01, var1 = "z", var2 = "y"),
    "Argument 'var1' must be a variable included in the fitted model."
  )
  expect_error(
    coev_plot_flowfield(object = m01, var1 = "x", var2 = 0:1),
    "Argument 'var2' must be a character string of length one."
  )
  expect_error(
    coev_plot_flowfield(object = m01, var1 = "x", var2 = "z"),
    "Argument 'var2' must be a variable included in the fitted model."
  )
  expect_error(
    coev_plot_flowfield(object = m01, var1 = "x", var2 = "x"),
    "Argument 'var1' and 'var2' must refer to different variables."
  )
  expect_error(
    coev_plot_flowfield(
      object = m01, var1 = "x", var2 = "y", nullclines = "hello"
    ),
    "Argument 'nullclines' must be logical."
  )
  expect_error(
    coev_plot_flowfield(
      object = m01, var1 = "x", var2 = "y", limits = "hello"
    ),
    "Argument 'limits' must be a numeric vector of length 2."
  )
  # should produce warning when at least three traits in the model
  expect_warning(
    coev_plot_flowfield(object = m01, var1 = "x", var2 = "y"),
    paste0(
      "Other traits were held constant at their median values to produce ",
      "this flowfield plot, which can potentially produce misleading ",
      "pictures of coevolutionary dynamics."
    )
  )
  # should run without error
  fun <- function(model, var1, var2,
                  nullclines = FALSE,
                  limits = c(-2.5, 2.5)) {
    suppressWarnings(
      coev_plot_flowfield(model, var1, var2, nullclines, limits)
    )
  }
  expect_no_error(fun(m01, "x", "y"))
  expect_no_error(fun(m02, "w", "x"))
  expect_no_error(fun(m03, "w", "x"))
  expect_no_error(fun(m04, "y", "z"))
  expect_no_error(fun(m05, "w", "x"))
  expect_no_error(fun(m06, "w", "x"))
  expect_no_error(fun(m07, "w", "x"))
  expect_no_error(fun(m08, "y", "x"))
  expect_no_error(fun(m09, "y", "x"))
  expect_no_error(fun(m10, "x", "y"))
  expect_no_error(fun(m01, "x", "y", nullclines = TRUE))
  expect_no_error(fun(m02, "w", "x", nullclines = TRUE))
  expect_no_error(fun(m03, "w", "x", nullclines = TRUE))
  expect_no_error(fun(m04, "y", "z", nullclines = TRUE))
  expect_no_error(fun(m05, "w", "x", nullclines = TRUE))
  expect_no_error(fun(m06, "w", "x", nullclines = TRUE))
  expect_no_error(fun(m07, "w", "x", nullclines = TRUE))
  expect_no_error(fun(m08, "y", "x", nullclines = TRUE))
  expect_no_error(fun(m09, "y", "x", nullclines = TRUE))
  expect_no_error(fun(m10, "x", "y", nullclines = TRUE))
  expect_no_error(fun(m01, "x", "y", limits = c(-3, 3)))
  expect_no_error(fun(m02, "w", "x", limits = c(-3, 3)))
  expect_no_error(fun(m03, "w", "x", limits = c(-3, 3)))
  expect_no_error(fun(m04, "y", "z", limits = c(-3, 3)))
  expect_no_error(fun(m05, "w", "x", limits = c(-3, 3)))
  expect_no_error(fun(m06, "w", "x", limits = c(-3, 3)))
  expect_no_error(fun(m07, "w", "x", limits = c(-3, 3)))
  expect_no_error(fun(m08, "y", "x", limits = c(-3, 3)))
  expect_no_error(fun(m09, "y", "x", limits = c(-3, 3)))
  expect_no_error(fun(m10, "x", "y", limits = c(-3, 3)))
})
