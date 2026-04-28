# Test that Stan and JAX log-densities agree at the same unconstrained
# point across all model configurations. This is the core correctness
# check for the JAX backend — analogous to the transpailer approach.
#
# A small smoke subset (ordered logistic prior-only, normal prior-only,
# normal full-likelihood) runs by default. The full 19-config suite is
# gated behind COEVOLVE_EXTENDED_TESTS=true to keep the default test
# suite under 10 min; each Stan compile in this suite takes ~30–60s.

skip_if_not(
  coevolve:::check_jax_available(),
  message = "JAX not available - skipping logp comparison tests"
)

#' @srrstats {G5.10} Flag extended tests
run_extended_tests <- identical(Sys.getenv("COEVOLVE_EXTENDED_TESTS"), "true")

# Helper: run compare_stan_jax_logprob and assert gradient agreement.
# Prior-only tests hit machine precision; likelihood tests accumulate
# floating-point errors through the tree traversal and MVN evaluation,
# so we use a looser tolerance (1e-2) for those.
expect_logp_agreement <- function(..., grad_tol = 1e-4,
                                  offset_sd_tol = grad_tol) {
  result <- coevolve:::compare_stan_jax_logprob(
    ..., n_points = 3L, seed = 1L, grad_tol = grad_tol
  )
  testthat::expect_lt(
    result$offset_sd, offset_sd_tol,
    label = "logp offset should be constant (sd ≈ 0)"
  )
  testthat::expect_lt(
    result$max_grad_diff, grad_tol,
    label = "gradient discrepancy"
  )
}

test_that("logp agrees: ordered logistic (authority)", {
  expect_logp_agreement(
    data = authority$data,
    variables = list(
      political_authority = "ordered_logistic",
      religious_authority = "ordered_logistic"
    ),
    id = "language",
    tree = authority$phylogeny,
    prior = list(A_offdiag = "normal(0, 2)")
  )
})

test_that("logp agrees: normal (simulated)", {
  withr::with_seed(1, {
    n <- 10
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = tree$tip.label,
      x = rnorm(n),
      y = rnorm(n)
    )
  })
  expect_logp_agreement(
    data = d,
    variables = list(x = "normal", y = "normal"),
    id = "id",
    tree = tree
  )
})

test_that("logp agrees: bernoulli_logit (simulated)", {
  skip_if_not(run_extended_tests)
  withr::with_seed(2, {
    n <- 10
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = tree$tip.label,
      x = rnorm(n),
      y = rbinom(n, 1, 0.5)
    )
  })
  expect_logp_agreement(
    data = d,
    variables = list(x = "normal", y = "bernoulli_logit"),
    id = "id",
    tree = tree
  )
})

test_that("logp agrees: poisson_softplus (simulated)", {
  skip_if_not(run_extended_tests)
  withr::with_seed(3, {
    n <- 10
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = tree$tip.label,
      x = rnorm(n),
      y = rpois(n, 5)
    )
  })
  expect_logp_agreement(
    data = d,
    variables = list(x = "normal", y = "poisson_softplus"),
    id = "id",
    tree = tree
  )
})

test_that("logp agrees: negative_binomial_softplus (simulated)", {
  skip_if_not(run_extended_tests)
  withr::with_seed(4, {
    n <- 10
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = tree$tip.label,
      x = rnorm(n),
      y = as.integer(rnbinom(n, size = 5, mu = 10))
    )
  })
  expect_logp_agreement(
    data = d,
    variables = list(x = "normal", y = "negative_binomial_softplus"),
    id = "id",
    tree = tree
  )
})

test_that("logp agrees: exact GP spatial control", {
  skip_if_not(run_extended_tests)
  expect_logp_agreement(
    data = authority$data,
    variables = list(
      political_authority = "ordered_logistic",
      religious_authority = "ordered_logistic"
    ),
    id = "language",
    tree = authority$phylogeny,
    lon_lat = authority$coordinates,
    prior = list(A_offdiag = "normal(0, 2)")
  )
})

test_that("logp agrees: HSGP spatial control", {
  skip_if_not(run_extended_tests)
  expect_logp_agreement(
    data = authority$data,
    variables = list(
      political_authority = "ordered_logistic",
      religious_authority = "ordered_logistic"
    ),
    id = "language",
    tree = authority$phylogeny,
    lon_lat = authority$coordinates,
    dist_k = 3,
    prior = list(A_offdiag = "normal(0, 2)")
  )
})

test_that("logp agrees: repeated measures (normal)", {
  skip_if_not(run_extended_tests)
  expect_logp_agreement(
    data = repeated$data,
    variables = list(x = "normal", y = "normal"),
    id = "species",
    tree = repeated$phylogeny
  )
})

test_that("logp agrees: multiphylo (2 trees)", {
  skip_if_not(run_extended_tests)
  tree2 <- c(authority$phylogeny, authority$phylogeny)
  class(tree2) <- "multiPhylo"
  expect_logp_agreement(
    data = authority$data,
    variables = list(
      political_authority = "ordered_logistic",
      religious_authority = "ordered_logistic"
    ),
    id = "language",
    tree = tree2,
    prior = list(A_offdiag = "normal(0, 2)")
  )
})

test_that("logp agrees: measurement error (normal with SE)", {
  skip_if_not(run_extended_tests)
  withr::with_seed(5, {
    n <- 10
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = tree$tip.label,
      x = rnorm(n),
      y = rnorm(n),
      x_se = rexp(n, 5),
      y_se = rexp(n, 5)
    )
  })
  expect_logp_agreement(
    data = d,
    variables = list(x = "normal", y = "normal"),
    id = "id",
    tree = tree,
    measurement_error = list(x = "x_se", y = "y_se")
  )
})

test_that("logp agrees: correlated drift disabled", {
  skip_if_not(run_extended_tests)
  withr::with_seed(6, {
    n <- 10
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = tree$tip.label,
      x = rnorm(n),
      y = rnorm(n)
    )
  })
  expect_logp_agreement(
    data = d,
    variables = list(x = "normal", y = "normal"),
    id = "id",
    tree = tree,
    estimate_correlated_drift = FALSE
  )
})

test_that("logp agrees: with effects_mat restricting cross-effects", {
  skip_if_not(run_extended_tests)
  withr::with_seed(7, {
    n <- 10
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = tree$tip.label,
      x = rnorm(n),
      y = rnorm(n)
    )
  })
  em <- matrix(c(TRUE, FALSE, TRUE, TRUE), nrow = 2,
               dimnames = list(c("x", "y"), c("x", "y")))
  expect_logp_agreement(
    data = d,
    variables = list(x = "normal", y = "normal"),
    id = "id",
    tree = tree,
    effects_mat = em
  )
})

# ------------------------------------------------------------------
# prior_only=FALSE: test the full likelihood
#
# Likelihood tests use a looser tolerance (1e-2) because the MVN
# evaluation through the tree traversal accumulates floating-point
# errors from matrix_exp, Cholesky, and einsum operations — these
# differ between Stan's C++ and JAX's XLA implementations at the
# ~1e-3 level even though the computations are mathematically
# equivalent.
# ------------------------------------------------------------------

test_that("logp agrees WITH likelihood: ordered logistic", {
  skip_if_not(run_extended_tests)
  expect_logp_agreement(
    data = authority$data,
    variables = list(
      political_authority = "ordered_logistic",
      religious_authority = "ordered_logistic"
    ),
    id = "language",
    tree = authority$phylogeny,
    prior = list(A_offdiag = "normal(0, 2)"),
    prior_only = FALSE,
    grad_tol = 1e-2
  )
})

test_that("logp agrees WITH likelihood: normal", {
  withr::with_seed(1, {
    n <- 10
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = tree$tip.label,
      x = rnorm(n),
      y = rnorm(n)
    )
  })
  expect_logp_agreement(
    data = d,
    variables = list(x = "normal", y = "normal"),
    id = "id",
    tree = tree,
    prior_only = FALSE,
    grad_tol = 1e-2
  )
})

test_that("logp agrees WITH likelihood: mixed normal + bernoulli", {
  skip_if_not(run_extended_tests)
  withr::with_seed(2, {
    n <- 10
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = tree$tip.label,
      x = rnorm(n),
      y = rbinom(n, 1, 0.5)
    )
  })
  expect_logp_agreement(
    data = d,
    variables = list(x = "normal", y = "bernoulli_logit"),
    id = "id",
    tree = tree,
    prior_only = FALSE,
    grad_tol = 1e-2
  )
})

test_that("logp agrees WITH likelihood: poisson_softplus", {
  skip_if_not(run_extended_tests)
  withr::with_seed(3, {
    n <- 10
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = tree$tip.label,
      x = rnorm(n),
      y = rpois(n, 5)
    )
  })
  expect_logp_agreement(
    data = d,
    variables = list(x = "normal", y = "poisson_softplus"),
    id = "id",
    tree = tree,
    prior_only = FALSE,
    grad_tol = 1e-2
  )
})

test_that("logp agrees WITH likelihood: exact GP", {
  skip_if_not(run_extended_tests)
  expect_logp_agreement(
    data = authority$data,
    variables = list(
      political_authority = "ordered_logistic",
      religious_authority = "ordered_logistic"
    ),
    id = "language",
    tree = authority$phylogeny,
    lon_lat = authority$coordinates,
    prior = list(A_offdiag = "normal(0, 2)"),
    prior_only = FALSE,
    grad_tol = 1e-2
  )
})

test_that("logp agrees WITH likelihood: repeated measures", {
  skip_if_not(run_extended_tests)
  expect_logp_agreement(
    data = repeated$data,
    variables = list(x = "normal", y = "normal"),
    id = "species",
    tree = repeated$phylogeny,
    prior_only = FALSE,
    grad_tol = 1e-2
  )
})

test_that("logp agrees WITH likelihood: multiphylo", {
  skip_if_not(run_extended_tests)
  tree2 <- c(authority$phylogeny, authority$phylogeny)
  class(tree2) <- "multiPhylo"
  expect_logp_agreement(
    data = authority$data,
    variables = list(
      political_authority = "ordered_logistic",
      religious_authority = "ordered_logistic"
    ),
    id = "language",
    tree = tree2,
    prior = list(A_offdiag = "normal(0, 2)"),
    prior_only = FALSE,
    grad_tol = 1e-2
  )
})
