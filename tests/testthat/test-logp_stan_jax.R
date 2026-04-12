# Test that Stan and JAX log-densities agree at the same unconstrained
# point across all model configurations. This is the core correctness
# check for the JAX backend — analogous to the transpailer approach.

skip_if_not(
  coevolve:::check_jax_available(),
  message = "JAX not available - skipping logp comparison tests"
)

# Helper: run compare_stan_jax_logprob and assert gradient agreement
expect_logp_agreement <- function(..., grad_tol = 1e-4) {
  result <- compare_stan_jax_logprob(..., n_points = 3L, seed = 1L,
                                     grad_tol = grad_tol)
  expect_lt(
    result$offset_sd, 1e-4,
    label = "logp offset should be constant (sd ≈ 0)"
  )
  expect_lt(
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
  withr::with_seed(4, {
    n <- 10
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = tree$tip.label,
      x = rnorm(n),
      y = rnbinom(n, size = 5, mu = 10)
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
  expect_logp_agreement(
    data = repeated$data,
    variables = list(x = "normal", y = "normal"),
    id = "species",
    tree = repeated$phylogeny
  )
})

test_that("logp agrees: multiphylo (2 trees)", {
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
  withr::with_seed(7, {
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
    effects_mat = matrix(c(TRUE, FALSE, TRUE, TRUE), nrow = 2)
  )
})
