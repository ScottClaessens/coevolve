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

# Shared 4-tip tree + data fixture ------------------------------------------
# Returns list(tree, d_normal, d_binary, d_count, d_ordered, d_positive, d_repeated)
.make_fixture <- function(seed = 1L) {
  withr::with_seed(seed, {
    n <- 4L
    tree <- ape::rcoal(n)
    tips <- tree$tip.label
    list(
      tree = tree,
      d_normal   = data.frame(id = tips, x = rnorm(n), y = rnorm(n)),
      d_binary   = data.frame(id = tips, x = rnorm(n), y = sample(0:1, n, replace = TRUE)),
      d_count    = data.frame(id = tips, x = rnorm(n),
                              # Ensure overdispersion: var(y) > mean(y)
                              y = c(0L, 0L, 5L, 10L)),
      d_ordered  = data.frame(
        id = tips, x = rnorm(n),
        y  = as.ordered(factor(c(1L, 2L, 3L, 2L), levels = 1:3))
      ),
      d_positive = data.frame(id = tips, x = rnorm(n), y = abs(rnorm(n)) + 0.1),
      d_repeated = data.frame(id = rep(tips, each = 2), x = rnorm(2 * n), y = rnorm(2 * n))
    )
  })
}

.cmp <- function(data, variables, tree, prior_only = TRUE,
                 estimate_correlated_drift = FALSE, seed = 7L, ...) {
  compare_stan_pymc_logprob(
    data = data, variables = variables, id = "id", tree = tree,
    prior_only = prior_only,
    estimate_correlated_drift = estimate_correlated_drift,
    seed = seed, tol_warn = Inf, ...
  )
}

# tol = 0.03: accommodates matrix exponential numerical precision (Padé
# 12-term + 8 squarings vs Stan's native matrix_exp) while still catching
# real parameterization bugs (LKJ Jacobian discrepancy was 0.15-0.51).
.expect_parity <- function(out, tol = 0.03) {
  expect_true(is.finite(out$logprob_stan))
  expect_true(is.finite(out$logprob_stan_adj))
  expect_true(is.finite(out$logprob_pymc))
  expect_lt(out$abs_diff, tol)
}

# Test 1: prior_only=TRUE, all-normal, diagonal Q
test_that("compare_stan_pymc_logprob: prior_only all-normal density parity", {
  .skip_compare()
  fx <- .make_fixture(42L)
  out <- .cmp(fx$d_normal, list(x = "normal", y = "normal"), fx$tree,
              prior_only = TRUE)
  .expect_parity(out)
  expect_true(out$n_params > 5L)
})

# Test 2: likelihood normal+ordered_logistic + missing data (forces terminal_drift in both)
test_that("compare_stan_pymc_logprob: likelihood mixed vars + missing data parity", {
  .skip_compare()
  fx <- .make_fixture(99L)
  fx$d_ordered$x[1] <- NA  # missing normal obs forces needs_terminal_drift=TRUE
  out <- .cmp(fx$d_ordered, list(x = "normal", y = "ordered_logistic"), fx$tree,
              prior_only = FALSE, seed = 11L)
  .expect_parity(out)
  # terminal_drift present in both; adjustment must be zero
  expect_equal(out$stan_terminal_drift_log_prior, 0.0)
})

# Test 3: prior_only=TRUE, correlated drift (LKJ Cholesky inverse transform)
test_that("compare_stan_pymc_logprob: prior_only correlated drift parity", {
  .skip_compare()
  fx <- .make_fixture(42L)
  out <- .cmp(fx$d_normal, list(x = "normal", y = "normal"), fx$tree,
              prior_only = TRUE, estimate_correlated_drift = TRUE, seed = 13L)
  .expect_parity(out)
})

# Test 4: likelihood bernoulli_logit (prior_only=FALSE so terminal_drift is
# constrained by std_normal() in both Stan and PyMC)
test_that("compare_stan_pymc_logprob: likelihood bernoulli_logit parity", {
  .skip_compare()
  fx <- .make_fixture(42L)
  out <- .cmp(fx$d_binary, list(x = "normal", y = "bernoulli_logit"), fx$tree,
              prior_only = FALSE, seed = 17L)
  .expect_parity(out)
  expect_equal(out$stan_terminal_drift_log_prior, 0.0)
})

# Test 5: likelihood bernoulli_logit (non-normal response; terminal_drift in both)
test_that("compare_stan_pymc_logprob: likelihood bernoulli_logit parity (seed 19)", {
  .skip_compare()
  fx <- .make_fixture(42L)
  out <- .cmp(fx$d_binary, list(x = "normal", y = "bernoulli_logit"), fx$tree,
              prior_only = FALSE, seed = 19L)
  .expect_parity(out)
  expect_equal(out$stan_terminal_drift_log_prior, 0.0)
})

# Test 6: likelihood poisson_softplus (prior_only=FALSE so terminal_drift constrained)
test_that("compare_stan_pymc_logprob: likelihood poisson_softplus parity", {
  .skip_compare()
  fx <- .make_fixture(42L)
  out <- .cmp(fx$d_count, list(x = "normal", y = "poisson_softplus"), fx$tree,
              prior_only = FALSE, seed = 23L)
  .expect_parity(out)
})

# Test 7: likelihood gamma_log (prior_only=FALSE)
test_that("compare_stan_pymc_logprob: likelihood gamma_log parity", {
  .skip_compare()
  fx <- .make_fixture(42L)
  out <- .cmp(fx$d_positive, list(x = "normal", y = "gamma_log"), fx$tree,
              prior_only = FALSE, seed = 29L)
  .expect_parity(out)
})

# Test 8: likelihood correlated drift + bernoulli (prior_only=FALSE; J=2, non-all-normal)
test_that("compare_stan_pymc_logprob: likelihood correlated drift + bernoulli parity", {
  .skip_compare()
  fx <- .make_fixture(42L)
  out <- .cmp(fx$d_binary, list(x = "normal", y = "bernoulli_logit"), fx$tree,
              prior_only = FALSE, estimate_correlated_drift = TRUE, seed = 31L)
  .expect_parity(out)
})

# Test 9: likelihood negative_binomial_softplus (prior_only=FALSE)
test_that("compare_stan_pymc_logprob: likelihood negative_binomial_softplus parity", {
  .skip_compare()
  fx <- .make_fixture(42L)
  out <- .cmp(fx$d_count, list(x = "normal", y = "negative_binomial_softplus"), fx$tree,
              prior_only = FALSE, seed = 37L)
  .expect_parity(out)
})

# Test 10: likelihood repeated measures (2 obs per tip).  prior_only=FALSE so
# terminal_drift ~ std_normal() is active in both Stan and PyMC.
test_that("compare_stan_pymc_logprob: likelihood repeated measures parity", {
  .skip_compare()
  fx <- .make_fixture(42L)
  out <- .cmp(fx$d_repeated, list(x = "normal", y = "normal"), fx$tree,
              prior_only = FALSE, seed = 41L)
  .expect_parity(out)
})

# Test 11: prior_only=TRUE, measurement_error (SE provided as column in data)
test_that("compare_stan_pymc_logprob: prior_only measurement_error parity", {
  .skip_compare()
  fx <- .make_fixture(42L)
  d <- fx$d_normal
  d$x_se <- 0.1
  d$y_se <- 0.1
  out <- .cmp(d, list(x = "normal", y = "normal"), fx$tree,
              prior_only = TRUE, measurement_error = list(x = "x_se", y = "y_se"),
              seed = 43L)
  .expect_parity(out)
})

# Test 12: prior_only=TRUE, effects_mat restriction
test_that("compare_stan_pymc_logprob: prior_only effects_mat restriction parity", {
  .skip_compare()
  fx <- .make_fixture(42L)
  em <- matrix(c(TRUE, FALSE, TRUE, TRUE), 2, 2,
               dimnames = list(c("x", "y"), c("x", "y")))
  out <- .cmp(fx$d_normal, list(x = "normal", y = "normal"), fx$tree,
              prior_only = TRUE, effects_mat = em, seed = 47L)
  .expect_parity(out)
})

# Test 13: repeated measures + measurement_error (per-obs SE)
test_that("compare_stan_pymc_logprob: repeated + ME parity", {
  .skip_compare()
  fx <- .make_fixture(42L)
  d <- fx$d_repeated
  d$x_se <- runif(nrow(d), 0.05, 0.2)
  d$y_se <- runif(nrow(d), 0.05, 0.2)
  out <- .cmp(d, list(x = "normal", y = "normal"), fx$tree,
              prior_only = FALSE,
              measurement_error = list(x = "x_se", y = "y_se"),
              seed = 53L)
  .expect_parity(out)
})

# Test 14: all-non-normal (exercises obs_lp=zeros path)
test_that("compare_stan_pymc_logprob: all-non-normal parity", {
  .skip_compare()
  fx <- .make_fixture(42L)
  d <- data.frame(
    id = fx$tree$tip.label,
    x  = sample(0:1, 4, replace = TRUE),
    y  = as.ordered(factor(c(1L, 2L, 3L, 2L), levels = 1:3))
  )
  vars <- list(x = "bernoulli_logit", y = "ordered_logistic")
  out <- .cmp(d, vars, fx$tree,
              prior_only = FALSE, seed = 59L)
  .expect_parity(out)
})
