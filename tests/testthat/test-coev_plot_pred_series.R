test_that("coev_plot_pred_series() produces expected errors and output", {
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
    coev_plot_pred_series(object = "fail"),
    paste0(
      "Argument 'object' must be a fitted coevolutionary model of class ",
      "'coevfit'."
    ),
    fixed = TRUE
  )
  expect_error(
    coev_plot_pred_series(object = m01, prob = 95),
    "Argument 'prob' must be a single numeric value between 0 and 1.",
    fixed = TRUE
  )
  # suppress warnings
  sw <- suppressWarnings
  # should run without error and produce list of ggplot objects
  expect_no_error(sw(coev_plot_pred_series(m01)))
  expect_no_error(sw(coev_plot_pred_series(m02)))
  expect_no_error(sw(coev_plot_pred_series(m03)))
  expect_no_error(sw(coev_plot_pred_series(m04)))
  expect_no_error(sw(coev_plot_pred_series(m05)))
  expect_no_error(sw(coev_plot_pred_series(m06)))
  expect_no_error(sw(coev_plot_pred_series(m07)))
  expect_no_error(sw(coev_plot_pred_series(m08)))
  expect_no_error(sw(coev_plot_pred_series(m09)))
  expect_no_error(sw(coev_plot_pred_series(m10)))
  expect_no_error(sw(coev_plot_pred_series(m01, stochastic = TRUE)))
  expect_no_error(sw(coev_plot_pred_series(m02, stochastic = TRUE)))
  expect_no_error(sw(coev_plot_pred_series(m03, stochastic = TRUE)))
  expect_no_error(sw(coev_plot_pred_series(m04, stochastic = TRUE)))
  expect_no_error(sw(coev_plot_pred_series(m05, stochastic = TRUE)))
  expect_no_error(sw(coev_plot_pred_series(m06, stochastic = TRUE)))
  expect_no_error(sw(coev_plot_pred_series(m07, stochastic = TRUE)))
  expect_no_error(sw(coev_plot_pred_series(m08, stochastic = TRUE)))
  expect_no_error(sw(coev_plot_pred_series(m09, stochastic = TRUE)))
  expect_no_error(sw(coev_plot_pred_series(m10, stochastic = TRUE)))
  expect_no_error(sw(coev_plot_pred_series(m01, ndraws = 1)))
  expect_no_error(sw(coev_plot_pred_series(m02, ndraws = 1)))
  expect_no_error(sw(coev_plot_pred_series(m03, ndraws = 1)))
  expect_no_error(sw(coev_plot_pred_series(m04, ndraws = 1)))
  expect_no_error(sw(coev_plot_pred_series(m05, ndraws = 1)))
  expect_no_error(sw(coev_plot_pred_series(m06, ndraws = 1)))
  expect_no_error(sw(coev_plot_pred_series(m07, ndraws = 1)))
  expect_no_error(sw(coev_plot_pred_series(m08, ndraws = 1)))
  expect_no_error(sw(coev_plot_pred_series(m09, ndraws = 1)))
  expect_no_error(sw(coev_plot_pred_series(m10, ndraws = 1)))
  expect_no_error(
    sw(
      coev_plot_pred_series(
        m01, intervention_values = list(u = NA, v = 0, w = 0, x = 0, y = 0)
      )
    )
  )
  expect_no_error(
    sw(coev_plot_pred_series(m02, intervention_values = list(w = NA, x = 0)))
  )
  expect_no_error(
    sw(coev_plot_pred_series(m03, intervention_values = list(w = NA, x = 0)))
  )
  expect_no_error(
    sw(coev_plot_pred_series(m04, intervention_values = list(y = NA, z = 0)))
  )
  expect_no_error(
    sw(coev_plot_pred_series(m05, intervention_values = list(w = NA, x = 0)))
  )
  expect_no_error(
    sw(coev_plot_pred_series(m06, intervention_values = list(w = NA, x = 0)))
  )
  expect_no_error(
    sw(coev_plot_pred_series(m07, intervention_values = list(w = NA, x = 0)))
  )
  expect_no_error(
    sw(coev_plot_pred_series(m08, intervention_values = list(x = NA, y = 0)))
  )
  expect_no_error(
    sw(coev_plot_pred_series(m09, intervention_values = list(x = NA, y = 0)))
  )
  expect_no_error(
    sw(coev_plot_pred_series(m10, intervention_values = list(x = NA, y = 0)))
  )
  expect_no_error(
    sw(
      coev_plot_pred_series(
        m01, eta_anc = list(u = 0, v = 0, w = 0, x = 0, y = 0)
      )
    )
  )
  expect_no_error(sw(coev_plot_pred_series(m02, eta_anc = list(w = 0, x = 0))))
  expect_no_error(sw(coev_plot_pred_series(m03, eta_anc = list(w = 0, x = 0))))
  expect_no_error(sw(coev_plot_pred_series(m04, eta_anc = list(y = 0, z = 0))))
  expect_no_error(sw(coev_plot_pred_series(m05, eta_anc = list(w = 0, x = 0))))
  expect_no_error(sw(coev_plot_pred_series(m06, eta_anc = list(w = 0, x = 0))))
  expect_no_error(sw(coev_plot_pred_series(m07, eta_anc = list(w = 0, x = 0))))
  expect_no_error(sw(coev_plot_pred_series(m08, eta_anc = list(x = 0, y = 0))))
  expect_no_error(sw(coev_plot_pred_series(m09, eta_anc = list(x = 0, y = 0))))
  expect_no_error(sw(coev_plot_pred_series(m10, eta_anc = list(x = 0, y = 0))))
  expect_true(methods::is(sw(coev_plot_pred_series(m01)), "ggplot"))
  expect_true(methods::is(sw(coev_plot_pred_series(m02)), "ggplot"))
  expect_true(methods::is(sw(coev_plot_pred_series(m03)), "ggplot"))
  expect_true(methods::is(sw(coev_plot_pred_series(m04)), "ggplot"))
  expect_true(methods::is(sw(coev_plot_pred_series(m05)), "ggplot"))
  expect_true(methods::is(sw(coev_plot_pred_series(m06)), "ggplot"))
  expect_true(methods::is(sw(coev_plot_pred_series(m07)), "ggplot"))
  expect_true(methods::is(sw(coev_plot_pred_series(m08)), "ggplot"))
  expect_true(methods::is(sw(coev_plot_pred_series(m09)), "ggplot"))
  expect_true(methods::is(sw(coev_plot_pred_series(m10)), "ggplot"))
})
