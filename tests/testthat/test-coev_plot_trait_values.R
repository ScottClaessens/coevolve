test_that("coev_plot_trait_values() produces expected errors and output", {
  # load models
  m1 <- readRDS(test_path("fixtures", "coevfit_example1.rds"))
  m2 <- readRDS(test_path("fixtures", "coevfit_example2.rds"))
  m3 <- readRDS(test_path("fixtures", "coevfit_example3.rds"))
  m4 <- readRDS(test_path("fixtures", "coevfit_example4.rds"))
  m5 <- readRDS(test_path("fixtures", "coevfit_example5.rds"))
  m6 <- readRDS(test_path("fixtures", "coevfit_example6.rds"))
  m7 <- readRDS(test_path("fixtures", "coevfit_example7.rds"))
  m8 <- readRDS(test_path("fixtures", "coevfit_example8.rds"))
  m9 <- readRDS(test_path("fixtures", "coevfit_example9.rds"))
  m1 <- reload_fit(m1, filename = "coevfit_example1-1.csv")
  m2 <- reload_fit(m2, filename = "coevfit_example2-1.csv")
  m3 <- reload_fit(m3, filename = "coevfit_example3-1.csv")
  m4 <- reload_fit(m4, filename = "coevfit_example4-1.csv")
  m5 <- reload_fit(m5, filename = "coevfit_example5-1.csv")
  m6 <- reload_fit(m6, filename = "coevfit_example6-1.csv")
  m7 <- reload_fit(m7, filename = "coevfit_example7-1.csv")
  m8 <- reload_fit(m8, filename = "coevfit_example8-1.csv")
  m9 <- reload_fit(m9, filename = "coevfit_example9-1.csv")
  # expect the following errors
  expect_error(
    coev_plot_trait_values(object = "fail"),
    "Argument 'object' must be a fitted coevolutionary model of class coevfit.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_trait_values(object = m1, variables = "fail"),
    paste0(
      "Argument 'variables' must be a character vector ",
      "of at least length 2."
    ),
    fixed = TRUE
  )
  expect_error(
    coev_plot_trait_values(object = m1, variables = c("x", "fail")),
    paste0(
      "Argument 'variables' contains variable names that are not ",
      "included in the fitted model."
    ),
    fixed = TRUE
  )
  expect_error(
    coev_plot_trait_values(object = m1, ndraws = "fail"),
    "Argument 'ndraws' must be numeric.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_trait_values(object = m1, ndraws = c(1, 2)),
    "Argument 'ndraws' must be a single integer.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_trait_values(object = m1, ndraws = 0),
    "Argument 'ndraws' must be between 1 and the total number of draws.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_trait_values(object = m1, tree_id = "fail"),
    "Argument 'tree_id' must be numeric.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_trait_values(object = m1, tree_id = c(1, 2)),
    "Argument 'tree_id' must be a single integer.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_trait_values(object = m1, tree_id = 0),
    "Argument 'tree_id' must be between 1 and the total number of trees.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_trait_values(object = m1, tree_id = m1$stan_data$N_tree + 1),
    "Argument 'tree_id' must be between 1 and the total number of trees.",
    fixed = TRUE
  )
  # should run without error and produce ggplot object
  fun <- function(model, variables = NULL, ndraws = 50, tree_id = NULL,
                  xlim = NULL, ylim = NULL) {
    suppressWarnings(
      coev_plot_trait_values(model, variables, ndraws, tree_id, xlim, ylim)
      )
  }
  expect_no_error(fun(m1))
  expect_no_error(fun(m2))
  expect_no_error(fun(m3))
  expect_no_error(fun(m4))
  expect_no_error(fun(m5))
  expect_no_error(fun(m6))
  expect_no_error(fun(m7))
  expect_no_error(fun(m8))
  expect_no_error(fun(m9))
  expect_true(methods::is(fun(m1), "patchwork"))
  expect_true(methods::is(fun(m2), "patchwork"))
  expect_true(methods::is(fun(m3), "patchwork"))
  expect_true(methods::is(fun(m4), "patchwork"))
  expect_true(methods::is(fun(m5), "patchwork"))
  expect_true(methods::is(fun(m6), "patchwork"))
  expect_true(methods::is(fun(m7), "patchwork"))
  expect_true(methods::is(fun(m7), "patchwork"))
  expect_true(methods::is(fun(m8), "patchwork"))
  expect_true(methods::is(fun(m9), "patchwork"))
  expect_no_error(fun(m1, variables = c("u","v")))
  expect_no_error(fun(m2, variables = c("w","x")))
  expect_no_error(fun(m3, variables = c("w","x")))
  expect_no_error(fun(m4, variables = c("y","z")))
  expect_no_error(fun(m5, variables = c("w","x")))
  expect_no_error(fun(m6, variables = c("w","x")))
  expect_no_error(fun(m7, variables = c("w","x")))
  expect_no_error(fun(m8, variables = c("x","y")))
  expect_no_error(fun(m9, variables = c("x","y")))
  expect_no_error(fun(m1, ndraws = 10))
  expect_no_error(fun(m2, ndraws = 10))
  expect_no_error(fun(m3, ndraws = 10))
  expect_no_error(fun(m4, ndraws = 10))
  expect_no_error(fun(m5, ndraws = 10))
  expect_no_error(fun(m6, ndraws = 10))
  expect_no_error(fun(m7, ndraws = 10))
  expect_no_error(fun(m8, ndraws = 10))
  expect_no_error(fun(m9, ndraws = 10))
  expect_no_error(fun(m1, tree_id = 1))
  expect_no_error(fun(m2, tree_id = 1))
  expect_no_error(fun(m3, tree_id = 1))
  expect_no_error(fun(m4, tree_id = 1))
  expect_no_error(fun(m5, tree_id = 1))
  expect_no_error(fun(m6, tree_id = 1))
  expect_no_error(fun(m7, tree_id = 1))
  expect_no_error(fun(m8, tree_id = 1))
  expect_no_error(fun(m9, tree_id = 1))
  expect_no_error(fun(m1, xlim = c(-1, 1), ylim = c(-1, 1)))
  expect_no_error(fun(m2, xlim = c(-1, 1), ylim = c(-1, 1)))
  expect_no_error(fun(m3, xlim = c(-1, 1), ylim = c(-1, 1)))
  expect_no_error(fun(m4, xlim = c(-1, 1), ylim = c(-1, 1)))
  expect_no_error(fun(m5, xlim = c(-1, 1), ylim = c(-1, 1)))
  expect_no_error(fun(m6, xlim = c(-1, 1), ylim = c(-1, 1)))
  expect_no_error(fun(m7, xlim = c(-1, 1), ylim = c(-1, 1)))
  expect_no_error(fun(m8, xlim = c(-1, 1), ylim = c(-1, 1)))
  expect_no_error(fun(m9, xlim = c(-1, 1), ylim = c(-1, 1)))
})
