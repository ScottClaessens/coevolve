test_that("coev_plot_flowfield() produces expected errors and output", {
  # load models
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
    coev_plot_flowfield(object = "fail", var1 = "x", var2 = "y"),
    "Argument 'object' must be a fitted coevolutionary model of class coevfit."
  )
  expect_error(
    coev_plot_flowfield(object = m1, var1 = 0:1, var2 = "y"),
    "Argument 'var1' must be a character string of length one."
  )
  expect_error(
    coev_plot_flowfield(object = m1, var1 = "z", var2 = "y"),
    "Argument 'var1' must be a variable included in the fitted model."
  )
  expect_error(
    coev_plot_flowfield(object = m1, var1 = "x", var2 = 0:1),
    "Argument 'var2' must be a character string of length one."
  )
  expect_error(
    coev_plot_flowfield(object = m1, var1 = "x", var2 = "z"),
    "Argument 'var2' must be a variable included in the fitted model."
  )
  expect_error(
    coev_plot_flowfield(object = m1, var1 = "x", var2 = "x"),
    "Argument 'var1' and 'var2' must refer to different variables."
  )
  expect_error(
    coev_plot_flowfield(
      object = m1, var1 = "x", var2 = "y", nullclines = "hello"
      ),
    "Argument 'nullclines' must be logical."
  )
  expect_error(
    coev_plot_flowfield(
      object = m1, var1 = "x", var2 = "y", limits = "hello"
    ),
    "Argument 'limits' must be a numeric vector of length 2."
  )
  # should run without error
  fun <- function(model, var1, var2,
                  nullclines = FALSE,
                  limits = c(-2.5, 2.5)) {
    suppressWarnings(
      coev_plot_flowfield(model, var1, var2, nullclines, limits)
      )
  }
  expect_no_error(fun(m1, "x", "y"))
  expect_no_error(fun(m2, "w", "x"))
  expect_no_error(fun(m3, "w", "x"))
  expect_no_error(fun(m4, "y", "z"))
  expect_no_error(fun(m5, "w", "x"))
  expect_no_error(fun(m6, "w", "x"))
  expect_no_error(fun(m7, "w", "x"))
  expect_no_error(fun(m8, "y", "x"))
  expect_no_error(fun(m1, "x", "y", nullclines = TRUE))
  expect_no_error(fun(m2, "w", "x", nullclines = TRUE))
  expect_no_error(fun(m3, "w", "x", nullclines = TRUE))
  expect_no_error(fun(m4, "y", "z", nullclines = TRUE))
  expect_no_error(fun(m5, "w", "x", nullclines = TRUE))
  expect_no_error(fun(m6, "w", "x", nullclines = TRUE))
  expect_no_error(fun(m7, "w", "x", nullclines = TRUE))
  expect_no_error(fun(m8, "y", "x", nullclines = TRUE))
  expect_no_error(fun(m1, "x", "y", limits = c(-3, 3)))
  expect_no_error(fun(m2, "w", "x", limits = c(-3, 3)))
  expect_no_error(fun(m3, "w", "x", limits = c(-3, 3)))
  expect_no_error(fun(m4, "y", "z", limits = c(-3, 3)))
  expect_no_error(fun(m5, "w", "x", limits = c(-3, 3)))
  expect_no_error(fun(m6, "w", "x", limits = c(-3, 3)))
  expect_no_error(fun(m7, "w", "x", limits = c(-3, 3)))
  expect_no_error(fun(m8, "y", "x", limits = c(-3, 3)))
})
