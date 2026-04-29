test_that("coev_plot_node_labels() produces expected errors and output", {
  m01 <- readRDS(test_path("fixtures", "coevfit_example_01.rds"))
  m01 <- reload_fit(m01, filename = "coevfit_example_01-1.csv")
  # expect errors for invalid inputs
  expect_error(
    coev_plot_node_labels(object = "fail"),
    "Argument 'object' must be a fitted coevolutionary model of class coevfit.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_node_labels(object = m01, tip_labels = "fail"),
    "Argument 'tip_labels' must be TRUE or FALSE.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_node_labels(object = m01, node_cex = "fail"),
    "Argument 'node_cex' must be a single positive number.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_node_labels(object = m01, node_cex = -1),
    "Argument 'node_cex' must be a single positive number.",
    fixed = TRUE
  )
  expect_error(
    coev_plot_node_labels(object = m01, node_cex = c(0.5, 0.6)),
    "Argument 'node_cex' must be a single positive number.",
    fixed = TRUE
  )
  # should run without error and return phylo invisibly
  expect_no_error(result <- coev_plot_node_labels(m01))
  expect_true(inherits(result, "phylo"))
  # tip_labels = TRUE should also work
  expect_no_error(coev_plot_node_labels(m01, tip_labels = TRUE))
  # custom node_cex and node_col
  expect_no_error(
    coev_plot_node_labels(m01, node_cex = 0.4, node_col = "blue")
  )
})

test_that("coev_plot_node_labels() handles multiPhylo trees", {
  m08 <- readRDS(test_path("fixtures", "coevfit_example_08.rds"))
  m08 <- reload_fit(m08, filename = "coevfit_example_08-1.csv")
  expect_true(inherits(m08$tree, "multiPhylo"))
  expect_no_error(result <- coev_plot_node_labels(m08))
  expect_true(inherits(result, "phylo"))
})
