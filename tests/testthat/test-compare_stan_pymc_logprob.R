# Shared skip helpers --------------------------------------------------------

.skip_compare <- function() {
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("posterior")
  skip_if_not(coevolve:::check_pymc_available(), message = "PyMC not available")
  stan_ok <- tryCatch(
    !is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)),
    error = function(e) FALSE
  )
  skip_if_not(stan_ok, message = "CmdStan not configured")
}

# Test 1: prior_only=TRUE, all-normal -- terminal_drift has flat prior in Stan,
# no terminal_drift in PyMC; no adjustment needed. Both densities should agree
# closely after the tip-edge z_drift correction.
test_that("compare_stan_pymc_logprob: prior_only all-normal density parity", {
  .skip_compare()

  withr::with_seed(42, {
    n <- 4L
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = tree$tip.label,
      x = rnorm(n),
      y = rnorm(n)
    )
  })

  out <- compare_stan_pymc_logprob(
    data = d,
    variables = list(x = "normal", y = "normal"),
    id = "id",
    tree = tree,
    prior_only = TRUE,
    estimate_correlated_drift = FALSE,
    iter_warmup = 120L,
    iter_sampling = 1L,
    chains = 1L,
    seed = 7L,
    tol_warn = Inf
  )

  expect_true(is.finite(out$logprob_stan))
  expect_true(is.finite(out$logprob_stan_adj))
  expect_true(is.finite(out$logprob_pymc))
  expect_true(out$n_params > 5L)
  expect_lt(out$abs_diff, 1.0)
})

# Test 2: prior_only=FALSE, mixed variables (normal + ordered_logistic) with a
# missing normal observation. This forces needs_terminal_drift=TRUE in both Stan
# and PyMC, so both backends have identical parameter structure. The terminal_drift
# adjustment should be zero; only tip-edge z_drift correction applies.
test_that("compare_stan_pymc_logprob: likelihood mixed vars + missing data parity", {
  .skip_compare()

  withr::with_seed(99, {
    n <- 4L
    tree <- ape::rcoal(n)
    x_vals <- rnorm(n)
    x_vals[1] <- NA  # one missing normal obs forces needs_terminal_drift=TRUE
    d <- data.frame(
      id = tree$tip.label,
      x = x_vals,
      y = sample(1:3, n, replace = TRUE)
    )
  })

  out <- compare_stan_pymc_logprob(
    data = d,
    variables = list(x = "normal", y = "ordered_logistic"),
    id = "id",
    tree = tree,
    prior_only = FALSE,
    estimate_correlated_drift = FALSE,
    iter_warmup = 150L,
    iter_sampling = 1L,
    chains = 1L,
    seed = 11L,
    tol_warn = Inf
  )

  expect_true(is.finite(out$logprob_stan))
  expect_true(is.finite(out$logprob_stan_adj))
  expect_true(is.finite(out$logprob_pymc))
  # terminal_drift adjustment must be zero: both backends have terminal_drift
  expect_equal(out$stan_terminal_drift_log_prior, 0.0)
  expect_lt(out$abs_diff, 1.0)
})
