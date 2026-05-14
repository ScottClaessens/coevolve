#' @srrstats {G5.10} Flag extended tests
run_extended_tests <- identical(Sys.getenv("COEVOLVE_EXTENDED_TESTS"), "true")

#' @srrstats {BS7.0, BS7.1} Prior recovery tests
test_that("coev_fit() recovers prior distribution", {
  skip_if_not(run_extended_tests)
  # get dummy dataset
  withr::with_seed(1, {
    n <- 5
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = tree$tip.label,
      x = rnorm(n),
      y = rnorm(n)
    )
  })
  # fit model with prior_only = TRUE
  fit <-
    suppressWarnings(
      coev_fit(
        data = d,
        variables = list(
          x = "normal",
          y = "normal"
        ),
        id = "id",
        tree = tree,
        prior = list(
          A_offdiag = "normal(2, 0.5)",
          b = "normal(-1, 0.5)"
        ),
        prior_only = TRUE,
        chains = 1,
        refresh = 0,
        seed = 1
      )
    )
  # estimates should recover prior within tolerance
  prior <- extract_samples(fit)
  expect_equal(median(prior$A[, 1, 2]), 2, tolerance = 0.1)
  expect_equal(median(prior$A[, 2, 1]), 2, tolerance = 0.1)
  expect_equal(median(prior$b[, 1]), -1, tolerance = 0.1)
  expect_equal(median(prior$b[, 2]), -1, tolerance = 0.1)
  expect_equal(sd(prior$A[, 1, 2]), 0.5, tolerance = 0.1)
  expect_equal(sd(prior$A[, 2, 1]), 0.5, tolerance = 0.1)
  expect_equal(sd(prior$b[, 1]), 0.5, tolerance = 0.1)
  expect_equal(sd(prior$b[, 2]), 0.5, tolerance = 0.1)
})
