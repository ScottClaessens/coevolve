test_that("coev_plot_predictive_check() produces expected errors and output", {
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
    coev_plot_predictive_check(object = "fail"),
    "Argument 'object' must be a fitted coevolutionary model of class coevfit."
  )
  expect_error(
    coev_plot_predictive_check(object = m01, variables = 0),
    paste0(
      "Argument 'variables' must be a character string or a ",
      "vector of character strings."
    ),
    fixed = TRUE
  )
  expect_error(
    coev_plot_predictive_check(object = m01, variables = "fail"),
    paste0(
      "Argument 'variables' contains variable names that are not ",
      "included in the fitted model."
    ),
    fixed = TRUE
  )
  expect_error(
    coev_plot_predictive_check(object = m01, ndraws = "fail"),
    "Argument 'ndraws' must be numeric.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_predictive_check(object = m01, ndraws = c(1, 2)),
    "Argument 'ndraws' must be a single integer.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_predictive_check(object = m01, ndraws = 0),
    "Argument 'ndraws' must be between 1 and the total number of draws.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_predictive_check(object = m01, ndraws = nrow(m01$fit$draws()) + 1),
    "Argument 'ndraws' must be between 1 and the total number of draws.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_predictive_check(object = m01, tree_id = "fail"),
    "Argument 'tree_id' must be numeric.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_predictive_check(object = m01, tree_id = c(1, 2)),
    "Argument 'tree_id' must be a single integer.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_predictive_check(object = m01, tree_id = 0),
    "Argument 'tree_id' must be between 1 and the total number of trees.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_predictive_check(object = m01, tree_id = m01$stan_data$N_tree + 1),
    "Argument 'tree_id' must be between 1 and the total number of trees.",
    fixed = TRUE
  )
  # suppress warnings
  SW <- suppressWarnings
  # should run without error and produce list of ggplot objects
  expect_no_error(SW(coev_plot_predictive_check(m01)))
  expect_no_error(SW(coev_plot_predictive_check(m02)))
  expect_no_error(SW(coev_plot_predictive_check(m03)))
  expect_no_error(SW(coev_plot_predictive_check(m04)))
  expect_no_error(SW(coev_plot_predictive_check(m05)))
  expect_no_error(SW(coev_plot_predictive_check(m06)))
  expect_no_error(SW(coev_plot_predictive_check(m07)))
  expect_no_error(SW(coev_plot_predictive_check(m08)))
  expect_no_error(SW(coev_plot_predictive_check(m09)))
  expect_no_error(SW(coev_plot_predictive_check(m10)))
  expect_no_error(SW(coev_plot_predictive_check(m01, variables = "x")))
  expect_no_error(SW(coev_plot_predictive_check(m02, variables = "w")))
  expect_no_error(SW(coev_plot_predictive_check(m03, variables = "w")))
  expect_no_error(SW(coev_plot_predictive_check(m04, variables = "y")))
  expect_no_error(SW(coev_plot_predictive_check(m05, variables = "w")))
  expect_no_error(SW(coev_plot_predictive_check(m06, variables = "w")))
  expect_no_error(SW(coev_plot_predictive_check(m07, variables = "w")))
  expect_no_error(SW(coev_plot_predictive_check(m08, variables = "x")))
  expect_no_error(SW(coev_plot_predictive_check(m09, variables = "x")))
  expect_no_error(SW(coev_plot_predictive_check(m10, variables = "x")))
  expect_no_error(SW(coev_plot_predictive_check(m01, ndraws = 1)))
  expect_no_error(SW(coev_plot_predictive_check(m02, ndraws = 1)))
  expect_no_error(SW(coev_plot_predictive_check(m03, ndraws = 1)))
  expect_no_error(SW(coev_plot_predictive_check(m04, ndraws = 1)))
  expect_no_error(SW(coev_plot_predictive_check(m05, ndraws = 1)))
  expect_no_error(SW(coev_plot_predictive_check(m06, ndraws = 1)))
  expect_no_error(SW(coev_plot_predictive_check(m07, ndraws = 1)))
  expect_no_error(SW(coev_plot_predictive_check(m08, ndraws = 1)))
  expect_no_error(SW(coev_plot_predictive_check(m09, ndraws = 1)))
  expect_no_error(SW(coev_plot_predictive_check(m10, ndraws = 1)))
  # should work with multiPhylo
  expect_no_error(SW(coev_plot_predictive_check(m08, tree_id = 2)))
  # should work as expected with missing data
  # lower limit of plot should not be -9999.45
  expect_false(
    isTRUE(
      all.equal(
        ggplot2::layer_scales(
          SW(coev_plot_predictive_check(m06, variables = "x"))[[1]]
          )$x$range$range[1],
        -9999.45
        )
      )
    )
})
