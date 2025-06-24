test_that("coev_plot_trait_values() produces expected errors and output", {
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
    coev_plot_trait_values(object = "fail"),
    "Argument 'object' must be a fitted coevolutionary model of class coevfit.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_trait_values(object = m01, variables = "fail"),
    paste0(
      "Argument 'variables' must be a character vector ",
      "of at least length 2."
    ),
    fixed = TRUE
  )
  expect_error(
    coev_plot_trait_values(object = m01, variables = c("x", "fail")),
    paste0(
      "Argument 'variables' contains variable names that are not ",
      "included in the fitted model."
    ),
    fixed = TRUE
  )
  expect_error(
    coev_plot_trait_values(object = m01, ndraws = "fail"),
    "Argument 'ndraws' must be numeric.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_trait_values(object = m01, ndraws = c(1, 2)),
    "Argument 'ndraws' must be a single integer.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_trait_values(object = m01, ndraws = 0),
    "Argument 'ndraws' must be between 1 and the total number of draws.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_trait_values(object = m01, tree_id = "fail"),
    "Argument 'tree_id' must be numeric.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_trait_values(object = m01, tree_id = c(1, 2)),
    "Argument 'tree_id' must be a single integer.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_trait_values(object = m01, tree_id = 0),
    "Argument 'tree_id' must be between 1 and the total number of trees.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_trait_values(object = m01, tree_id = m01$stan_data$N_tree + 1),
    "Argument 'tree_id' must be between 1 and the total number of trees.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_trait_values(object = m01, tree_id = 1, xlim = 0),
    "Argument 'xlim' must be a numeric vector of length 2.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_trait_values(object = m01, tree_id = 1, xlim = c(0, 1), ylim = 0),
    "Argument 'ylim' must be a numeric vector of length 2.",
    fixed = TRUE
  )
  # should run without error and produce ggplot object
  fun <- function(model, variables = NULL, ndraws = 50, tree_id = NULL,
                  xlim = NULL, ylim = NULL) {
    suppressWarnings(
      coev_plot_trait_values(model, variables, ndraws, tree_id, xlim, ylim)
    )
  }
  expect_no_error(fun(m01))
  expect_no_error(fun(m02))
  expect_no_error(fun(m03))
  expect_no_error(fun(m04))
  expect_no_error(fun(m05))
  expect_no_error(fun(m06))
  expect_no_error(fun(m07))
  expect_no_error(fun(m08))
  expect_no_error(fun(m09))
  expect_no_error(fun(m10))
  expect_true(methods::is(fun(m01), "patchwork"))
  expect_true(methods::is(fun(m02), "patchwork"))
  expect_true(methods::is(fun(m03), "patchwork"))
  expect_true(methods::is(fun(m04), "patchwork"))
  expect_true(methods::is(fun(m05), "patchwork"))
  expect_true(methods::is(fun(m06), "patchwork"))
  expect_true(methods::is(fun(m07), "patchwork"))
  expect_true(methods::is(fun(m07), "patchwork"))
  expect_true(methods::is(fun(m08), "patchwork"))
  expect_true(methods::is(fun(m09), "patchwork"))
  expect_true(methods::is(fun(m10), "patchwork"))
  expect_no_error(fun(m01, variables = c("u", "v")))
  expect_no_error(fun(m02, variables = c("w", "x")))
  expect_no_error(fun(m03, variables = c("w", "x")))
  expect_no_error(fun(m04, variables = c("y", "z")))
  expect_no_error(fun(m05, variables = c("w", "x")))
  expect_no_error(fun(m06, variables = c("w", "x")))
  expect_no_error(fun(m07, variables = c("w", "x")))
  expect_no_error(fun(m08, variables = c("x", "y")))
  expect_no_error(fun(m09, variables = c("x", "y")))
  expect_no_error(fun(m10, variables = c("x", "y")))
  expect_no_error(fun(m01, ndraws = 10))
  expect_no_error(fun(m02, ndraws = 10))
  expect_no_error(fun(m03, ndraws = 10))
  expect_no_error(fun(m04, ndraws = 10))
  expect_no_error(fun(m05, ndraws = 10))
  expect_no_error(fun(m06, ndraws = 10))
  expect_no_error(fun(m07, ndraws = 10))
  expect_no_error(fun(m08, ndraws = 10))
  expect_no_error(fun(m09, ndraws = 10))
  expect_no_error(fun(m10, ndraws = 10))
  expect_no_error(fun(m01, tree_id = 1))
  expect_no_error(fun(m02, tree_id = 1))
  expect_no_error(fun(m03, tree_id = 1))
  expect_no_error(fun(m04, tree_id = 1))
  expect_no_error(fun(m05, tree_id = 1))
  expect_no_error(fun(m06, tree_id = 1))
  expect_no_error(fun(m07, tree_id = 1))
  expect_no_error(fun(m08, tree_id = 1))
  expect_no_error(fun(m09, tree_id = 1))
  expect_no_error(fun(m10, tree_id = 1))
  expect_no_error(fun(m01, xlim = c(-1, 1), ylim = c(-1, 1)))
  expect_no_error(fun(m02, xlim = c(-1, 1), ylim = c(-1, 1)))
  expect_no_error(fun(m03, xlim = c(-1, 1), ylim = c(-1, 1)))
  expect_no_error(fun(m04, xlim = c(-1, 1), ylim = c(-1, 1)))
  expect_no_error(fun(m05, xlim = c(-1, 1), ylim = c(-1, 1)))
  expect_no_error(fun(m06, xlim = c(-1, 1), ylim = c(-1, 1)))
  expect_no_error(fun(m07, xlim = c(-1, 1), ylim = c(-1, 1)))
  expect_no_error(fun(m08, xlim = c(-1, 1), ylim = c(-1, 1)))
  expect_no_error(fun(m09, xlim = c(-1, 1), ylim = c(-1, 1)))
  expect_no_error(fun(m10, xlim = c(-1, 1), ylim = c(-1, 1)))
})
