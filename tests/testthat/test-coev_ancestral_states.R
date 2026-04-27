test_that("coev_ancestral_states() produces expected errors", {
  m01 <- readRDS(test_path("fixtures", "coevfit_example_01.rds"))
  m01 <- reload_fit(m01, filename = "coevfit_example_01-1.csv")
  expect_error(
    coev_ancestral_states(object = "fail"),
    "Argument 'object' must be a fitted coevolutionary model of class coevfit.",
    fixed = TRUE
  )
  expect_error(
    coev_ancestral_states(object = m01, variables = 123),
    "Argument 'variables' must be a character vector.",
    fixed = TRUE
  )
  expect_error(
    coev_ancestral_states(object = m01, variables = "nonexistent"),
    paste0(
      "Argument 'variables' contains variable names that are not ",
      "included in the fitted model."
    ),
    fixed = TRUE
  )
  expect_error(
    coev_ancestral_states(object = m01, nodes = "bad_value"),
    "Argument 'nodes' must be",
    fixed = TRUE
  )
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
  expect_error(
    coev_ancestral_states(object = m01, scale = "bad"),
    "Argument 'scale' must be",
    fixed = TRUE
  )
  expect_error(
    coev_ancestral_states(object = m01, summary = "fail"),
    "Argument 'summary' must be logical.",
    fixed = TRUE
  )
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

test_that("coev_ancestral_states() returns correct structure", {
  m01 <- readRDS(test_path("fixtures", "coevfit_example_01.rds"))
  m01 <- reload_fit(m01, filename = "coevfit_example_01-1.csv")
  sw <- suppressWarnings
  result <- sw(coev_ancestral_states(m01))

  expect_true(is.data.frame(result))
  expect_true(all(
    c("node", "variable", "estimate", "lower", "upper") %in%
      names(result)
  ))
  expect_true("clade_pp" %in% names(result))
  expect_true(all(result$clade_pp == 1.0))
  n_tips <- m01$stan_data$N_tips
  n_internal <- n_tips - 1
  n_vars <- length(m01$variables)
  expect_equal(nrow(result), n_internal * n_vars)
  expect_true(all(result$node > n_tips))
  expect_equal(
    sort(unique(result$variable)),
    sort(names(m01$variables))
  )
  expect_true(all(is.finite(result$estimate)))
  expect_true(all(result$lower <= result$estimate))
  expect_true(all(result$estimate <= result$upper))
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
  expect_equal(nrow(result), n_seg * n_vars)
  expect_true(any(result$node <= n_tips))
  expect_true(any(result$node > n_tips))
})

test_that("coev_ancestral_states() nodes = integer vector works", {
  m01 <- readRDS(test_path("fixtures", "coevfit_example_01.rds"))
  m01 <- reload_fit(m01, filename = "coevfit_example_01-1.csv")
  sw <- suppressWarnings
  n_tips <- m01$stan_data$N_tips
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
  expect_true(is.list(result))
  expect_true("draws" %in% names(result))
  expect_true("ref_tree" %in% names(result))
  expect_true("variable_names" %in% names(result))
  expect_equal(length(dim(result$draws)), 3)
  n_tips <- m01$stan_data$N_tips
  n_internal <- n_tips - 1
  n_vars <- length(m01$variables)
  expect_equal(dim(result$draws)[2], n_internal)
  expect_equal(dim(result$draws)[3], n_vars)
  expect_true(all(is.finite(result$draws)))
})

test_that("coev_ancestral_states() variables subsets correctly", {
  m01 <- readRDS(test_path("fixtures", "coevfit_example_01.rds"))
  m01 <- reload_fit(m01, filename = "coevfit_example_01-1.csv")
  sw <- suppressWarnings
  result <- sw(coev_ancestral_states(m01, variables = c("u", "v")))
  expect_equal(sort(unique(result$variable)), sort(c("u", "v")))
  n_tips <- m01$stan_data$N_tips
  n_internal <- n_tips - 1
  expect_equal(nrow(result), n_internal * 2)
})

test_that("coev_ancestral_states() prob changes CI width", {
  m01 <- readRDS(test_path("fixtures", "coevfit_example_01.rds"))
  m01 <- reload_fit(m01, filename = "coevfit_example_01-1.csv")
  sw <- suppressWarnings
  narrow <- sw(coev_ancestral_states(m01, prob = 0.5))
  wide <- sw(coev_ancestral_states(m01, prob = 0.99))
  ci_width_narrow <- narrow$upper - narrow$lower
  ci_width_wide <- wide$upper - wide$lower
  expect_true(all(ci_width_wide >= ci_width_narrow - 1e-10))
})

test_that("coev_ancestral_states() bernoulli response is [0,1]", {
  m09 <- readRDS(test_path("fixtures", "coevfit_example_09.rds"))
  m09 <- reload_fit(m09, filename = "coevfit_example_09-1.csv")
  sw <- suppressWarnings
  result <- sw(coev_ancestral_states(m09, scale = "response"))
  expect_true(all(result$estimate >= 0 & result$estimate <= 1))
  expect_true(all(result$lower >= 0))
  expect_true(all(result$upper <= 1))
})

test_that("coev_ancestral_states() poisson response is positive", {
  m04 <- readRDS(test_path("fixtures", "coevfit_example_04.rds"))
  m04 <- reload_fit(m04, filename = "coevfit_example_04-1.csv")
  sw <- suppressWarnings
  result <- sw(coev_ancestral_states(m04, scale = "response"))
  expect_true(all(result$estimate > 0))
  expect_true(all(result$lower >= 0))
})

test_that("coev_ancestral_states() normal response = latent", {
  m10 <- readRDS(test_path("fixtures", "coevfit_example_10.rds"))
  m10 <- reload_fit(m10, filename = "coevfit_example_10-1.csv")
  sw <- suppressWarnings
  latent <- sw(coev_ancestral_states(m10, scale = "latent"))
  response <- sw(coev_ancestral_states(m10, scale = "response"))
  expect_equal(latent$estimate, response$estimate, tolerance = 1e-10)
  expect_equal(latent$lower, response$lower, tolerance = 1e-10)
  expect_equal(latent$upper, response$upper, tolerance = 1e-10)
})

test_that("coev_ancestral_states() gamma response is positive", {
  m01 <- readRDS(test_path("fixtures", "coevfit_example_01.rds"))
  m01 <- reload_fit(m01, filename = "coevfit_example_01-1.csv")
  sw <- suppressWarnings
  result <- sw(coev_ancestral_states(
    m01, variables = "u", scale = "response"
  ))
  expect_true(all(result$estimate > 0))
  expect_true(all(result$lower > 0))
})

test_that("coev_ancestral_states() ordinal response has probs", {
  m02 <- readRDS(test_path("fixtures", "coevfit_example_02.rds"))
  m02 <- reload_fit(m02, filename = "coevfit_example_02-1.csv")
  sw <- suppressWarnings
  result <- sw(coev_ancestral_states(
    m02, variables = "x", scale = "response"
  ))
  expect_true(any(grepl("^prob_", names(result))))
  prob_cols <- grep("^prob_", names(result), value = TRUE)
  for (col in prob_cols) {
    expect_true(all(result[[col]] >= 0 & result[[col]] <= 1))
  }
  # medians of individual category probabilities need not sum to

  # exactly 1 (median doesn't distribute over sums), but should
  # be close
  prob_sums <- rowSums(result[, prob_cols])
  expect_true(all(prob_sums > 0.5 & prob_sums < 1.5))
})

test_that("coev_ancestral_states() multi-tree tree_id works", {
  m08 <- readRDS(test_path("fixtures", "coevfit_example_08.rds"))
  m08 <- reload_fit(m08, filename = "coevfit_example_08-1.csv")
  sw <- suppressWarnings
  result1 <- sw(coev_ancestral_states(m08, tree_id = 1))
  expect_true(is.data.frame(result1))
  expect_true(all(result1$clade_pp == 1.0))
  result2 <- sw(coev_ancestral_states(m08, tree_id = 2))
  expect_true(is.data.frame(result2))
})

test_that("coev_ancestral_states() multi-tree default = tree 1", {
  m08 <- readRDS(test_path("fixtures", "coevfit_example_08.rds"))
  m08 <- reload_fit(m08, filename = "coevfit_example_08-1.csv")
  sw <- suppressWarnings
  result_default <- sw(coev_ancestral_states(m08))
  result_t1 <- sw(coev_ancestral_states(m08, tree_id = 1))
  expect_equal(result_default$estimate, result_t1$estimate)
})

test_that("coev_ancestral_states() works across fixture models", {
  sw <- suppressWarnings
  fixtures <- c(
    "01", "02", "04", "05", "06", "08", "09", "10"
  )
  for (id in fixtures) {
    rds <- paste0("coevfit_example_", id, ".rds")
    csv <- paste0("coevfit_example_", id, "-1.csv")
    m <- readRDS(test_path("fixtures", rds))
    m <- reload_fit(m, filename = csv)
    expect_no_error(sw(coev_ancestral_states(m)))
    result <- sw(coev_ancestral_states(m))
    expect_true(is.data.frame(result))
    expect_true(all(is.finite(result$estimate)))
    expect_no_error(
      sw(coev_ancestral_states(m, scale = "response"))
    )
  }
})
