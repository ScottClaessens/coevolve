test_that("coev_ancestral_states() produces expected errors", {
  # load single-tree model (5 distributions) and multi-tree model
  m01 <- readRDS(test_path("fixtures", "coevfit_example_01.rds"))
  m01 <- reload_fit(m01, filename = "coevfit_example_01-1.csv")
  # object must be coevfit
  expect_error(
    coev_ancestral_states(object = "fail"),
    "Argument 'object' must be a fitted coevolutionary model of class coevfit.",
    fixed = TRUE
  )
  # variables must be character
  expect_error(
    coev_ancestral_states(object = m01, variables = 123),
    "Argument 'variables' must be a character vector.",
    fixed = TRUE
  )
  # variables must be in model
  expect_error(
    coev_ancestral_states(object = m01, variables = "nonexistent"),
    paste0(
      "Argument 'variables' contains variable names that are not ",
      "included in the fitted model."
    ),
    fixed = TRUE
  )
  # nodes must be valid
  expect_error(
    coev_ancestral_states(object = m01, nodes = "bad_value"),
    "Argument 'nodes' must be \"internal\", \"all\", or an integer vector.",
    fixed = TRUE
  )
  # tree_id must be valid integer
  expect_error(
    coev_ancestral_states(object = m01, tree_id = "fail"),
    "Argument 'tree_id' must be a single integer.",
    fixed = TRUE
  )
  expect_error(
    coev_ancestral_states(object = m01, tree_id = 99),
    "Argument 'tree_id' must be between 1 and the number of trees.",
    fixed = TRUE
  )
  # scale must be valid
  expect_error(
    coev_ancestral_states(object = m01, scale = "bad"),
    "Argument 'scale' must be \"latent\" or \"response\".",
    fixed = TRUE
  )
  # summary must be logical
  expect_error(
    coev_ancestral_states(object = m01, summary = "fail"),
    "Argument 'summary' must be logical.",
    fixed = TRUE
  )
  # prob must be valid
  expect_error(
    coev_ancestral_states(object = m01, prob = 1.5),
    "Argument 'prob' must be a single numeric value between 0 and 1.",
    fixed = TRUE
  )
  expect_error(
    coev_ancestral_states(object = m01, prob = -0.1),
    "Argument 'prob' must be a single numeric value between 0 and 1.",
    fixed = TRUE
  )
})

test_that("coev_ancestral_states() returns correct structure (single tree, latent)", {
  # m01: single tree, 5 tips, 5 variables
  m01 <- readRDS(test_path("fixtures", "coevfit_example_01.rds"))
  m01 <- reload_fit(m01, filename = "coevfit_example_01-1.csv")
  sw <- suppressWarnings
  # default: internal nodes only, summary = TRUE, latent scale
  result <- sw(coev_ancestral_states(m01))
  # should be a data frame / tibble

  expect_true(is.data.frame(result))
  # should have the correct columns
  expect_true(all(c("node", "variable", "estimate", "lower", "upper")
                  %in% names(result)))
  # should have clade_pp = 1 for single tree
  expect_true("clade_pp" %in% names(result))
  expect_true(all(result$clade_pp == 1.0))
  # number of rows: internal nodes (N_tips - 1 = 4) * variables (5) = 20
  n_tips <- m01$stan_data$N_tips
  n_internal <- n_tips - 1
  n_vars <- length(m01$variables)
  expect_equal(nrow(result), n_internal * n_vars)
  # node IDs should be > N_tips (internal nodes in ape convention)
  expect_true(all(result$node > n_tips))
  # variable names should match the model
  expect_equal(sort(unique(result$variable)), sort(names(m01$variables)))
  # estimates should be finite
  expect_true(all(is.finite(result$estimate)))
  # lower should be <= estimate, estimate <= upper
  expect_true(all(result$lower <= result$estimate))
  expect_true(all(result$estimate <= result$upper))
  # ref_tree attribute should be set
  expect_true(!is.null(attr(result, "ref_tree")))
  expect_true(inherits(attr(result, "ref_tree"), "phylo"))
})

test_that("coev_ancestral_states() nodes = 'all' includes tips", {
  m01 <- readRDS(test_path("fixtures", "coevfit_example_01.rds"))
  m01 <- reload_fit(m01, filename = "coevfit_example_01-1.csv")
  sw <- suppressWarnings
  result <- sw(coev_ancestral_states(m01, nodes = "all"))
  n_tips <- m01$stan_data$N_tips
  n_seg <- m01$stan_data$N_seg
  n_vars <- length(m01$variables)
  # should include all nodes (tips + internal)
  expect_equal(nrow(result), n_seg * n_vars)
  # some node IDs should be <= N_tips (tips)
  expect_true(any(result$node <= n_tips))
  # some node IDs should be > N_tips (internal)
  expect_true(any(result$node > n_tips))
})

test_that("coev_ancestral_states() nodes = integer vector works", {
  m01 <- readRDS(test_path("fixtures", "coevfit_example_01.rds"))
  m01 <- reload_fit(m01, filename = "coevfit_example_01-1.csv")
  sw <- suppressWarnings
  n_tips <- m01$stan_data$N_tips
  # request specific internal nodes
  target_nodes <- c(n_tips + 1, n_tips + 2)
  result <- sw(coev_ancestral_states(m01, nodes = target_nodes))
  n_vars <- length(m01$variables)
  expect_equal(nrow(result), length(target_nodes) * n_vars)
  expect_equal(sort(unique(result$node)), sort(target_nodes))
})

test_that("coev_ancestral_states() summary = FALSE returns draws", {
  m01 <- readRDS(test_path("fixtures", "coevfit_example_01.rds"))
  m01 <- reload_fit(m01, filename = "coevfit_example_01-1.csv")
  sw <- suppressWarnings
  result <- sw(coev_ancestral_states(m01, summary = FALSE))
  # should be a list
  expect_true(is.list(result))
  expect_true("draws" %in% names(result))
  expect_true("ref_tree" %in% names(result))
  expect_true("variable_names" %in% names(result))
  # draws should be a 3D array: [draws, nodes, variables]
  expect_equal(length(dim(result$draws)), 3)
  n_tips <- m01$stan_data$N_tips
  n_internal <- n_tips - 1
  n_vars <- length(m01$variables)
  expect_equal(dim(result$draws)[2], n_internal) # internal nodes only
  expect_equal(dim(result$draws)[3], n_vars)
  # all values should be finite
  expect_true(all(is.finite(result$draws)))
})

test_that("coev_ancestral_states() variables argument subsets correctly", {
  m01 <- readRDS(test_path("fixtures", "coevfit_example_01.rds"))
  m01 <- reload_fit(m01, filename = "coevfit_example_01-1.csv")
  sw <- suppressWarnings
  result <- sw(coev_ancestral_states(m01, variables = c("u", "v")))
  expect_equal(sort(unique(result$variable)), sort(c("u", "v")))
  n_tips <- m01$stan_data$N_tips
  n_internal <- n_tips - 1
  expect_equal(nrow(result), n_internal * 2)
})

test_that("coev_ancestral_states() prob argument changes CI width", {
  m01 <- readRDS(test_path("fixtures", "coevfit_example_01.rds"))
  m01 <- reload_fit(m01, filename = "coevfit_example_01-1.csv")
  sw <- suppressWarnings
  narrow <- sw(coev_ancestral_states(m01, prob = 0.5))
  wide <- sw(coev_ancestral_states(m01, prob = 0.99))
  # wider prob should give wider intervals
  ci_width_narrow <- narrow$upper - narrow$lower
  ci_width_wide <- wide$upper - wide$lower
  expect_true(all(ci_width_wide >= ci_width_narrow - 1e-10))
})

test_that("coev_ancestral_states() works across all fixture models", {
  sw <- suppressWarnings
  models <- list(
    m01 = list(rds = "coevfit_example_01.rds", csv = "coevfit_example_01-1.csv"),
    m02 = list(rds = "coevfit_example_02.rds", csv = "coevfit_example_02-1.csv"),
    m04 = list(rds = "coevfit_example_04.rds", csv = "coevfit_example_04-1.csv"),
    m05 = list(rds = "coevfit_example_05.rds", csv = "coevfit_example_05-1.csv"),
    m06 = list(rds = "coevfit_example_06.rds", csv = "coevfit_example_06-1.csv"),
    m09 = list(rds = "coevfit_example_09.rds", csv = "coevfit_example_09-1.csv"),
    m10 = list(rds = "coevfit_example_10.rds", csv = "coevfit_example_10-1.csv")
  )
  for (name in names(models)) {
    m <- readRDS(test_path("fixtures", models[[name]]$rds))
    m <- reload_fit(m, filename = models[[name]]$csv)
    expect_no_error(sw(coev_ancestral_states(m)), label = name)
    result <- sw(coev_ancestral_states(m))
    expect_true(is.data.frame(result), label = name)
    expect_true(all(is.finite(result$estimate)), label = name)
  }
})
