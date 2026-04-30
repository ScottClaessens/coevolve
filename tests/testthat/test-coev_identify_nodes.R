test_that("coev_identify_nodes() validates inputs", {
  expect_error(
    coev_identify_nodes(tree = "fail"),
    "Argument 'tree' must be a 'phylo' object.",
    fixed = TRUE
  )
  tr <- ape::rcoal(5)
  expect_error(
    coev_identify_nodes(tree = tr, n = "fail"),
    "Argument 'n' must be a single positive integer.",
    fixed = TRUE
  )
  expect_error(
    coev_identify_nodes(tree = tr, n = 0),
    "Argument 'n' must be a single positive integer.",
    fixed = TRUE
  )
  expect_error(
    coev_identify_nodes(tree = tr, n = c(1, 2)),
    "Argument 'n' must be a single positive integer.",
    fixed = TRUE
  )
  expect_error(
    coev_identify_nodes(tree = tr, n = 1.5),
    "Argument 'n' must be a single positive integer.",
    fixed = TRUE
  )
})
