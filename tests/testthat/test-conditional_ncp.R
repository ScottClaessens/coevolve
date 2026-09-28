# Tests for the conditionally non-centred terminal drift (#124).
#
# In models with Gaussian variables and no repeated observations, the
# latent terminal drift of non-Gaussian variables (and of missing Gaussian
# values) is parameterised as a standard normal innovation conditional on
# the observed Gaussian residuals. This is a change of variables only: the
# equivalence tests below check that the log density matches the centred
# parameterisation up to the Jacobian and that the pointwise log likelihood
# is unchanged. A smoke subset runs by default; the remaining configurations
# are gated behind COEVOLVE_EXTENDED_TESTS=true because each configuration
# compiles two Stan models.

#' @srrstats {G5.10} Flag extended tests
run_extended_tests <- identical(Sys.getenv("COEVOLVE_EXTENDED_TESTS"), "true")

sim_mixed_data <- function(n = 12, seed = 1) {
  withr::with_seed(seed, {
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = tree$tip.label,
      x = rnorm(n),
      w = rnorm(n),
      z = rbinom(n, 1, 0.5),
      o = ordered(sample(1:3, n, replace = TRUE)),
      p = rpois(n, 3),
      x_se = rexp(n, 5),
      lon = runif(n, -10, 10),
      lat = runif(n, -10, 10)
    )
  })
  list(tree = tree, data = d)
}

expect_ncp_equivalent <- function(...) {
  res <- compare_ncp_centred(...)
  testthat::expect_false(identical(res$code_ncp, res$code_centred))
  testthat::expect_equal(
    res$lp_ncp, res$lp_centred_plus_jacobian,
    tolerance = 1e-8, label = "non-centred log density"
  )
  testthat::expect_equal(
    res$log_lik_ncp, res$log_lik_centred,
    tolerance = 1e-10, label = "non-centred pointwise log likelihood"
  )
}

test_that("use_conditional_ncp() selects the intended configurations", {
  sim <- sim_mixed_data()
  d <- sim$data
  d_miss <- d
  d_miss$x[1] <- NA
  d_dup <- rbind(d, d)
  mixed <- c("normal", "bernoulli_logit")
  expect_true(use_conditional_ncp(d, c("x", "z"), mixed, "id"))
  expect_true(use_conditional_ncp(d_miss, "x", "normal", "id"))
  expect_false(use_conditional_ncp(d, "x", "normal", "id"))
  expect_false(use_conditional_ncp(d, c("x", "w"), c("normal", "normal"), "id"))
  expect_false(use_conditional_ncp(d, "z", "bernoulli_logit", "id"))
  expect_false(use_conditional_ncp(d_dup, c("x", "z"), mixed, "id"))
})

test_that("Stan code outside the non-centred scope is unchanged", {
  sim <- sim_mixed_data()
  d_dup <- rbind(sim$data, sim$data)
  configs <- list(
    list(data = sim$data, variables = list(x = "normal", w = "normal")),
    list(data = sim$data, variables = list(z = "bernoulli_logit",
                                           o = "ordered_logistic")),
    list(data = d_dup, variables = list(x = "normal", z = "bernoulli_logit")),
    list(data = d_dup, variables = list(x = "normal", z = "bernoulli_logit"),
         estimate_residual = FALSE)
  )
  for (cfg in configs) {
    args <- c(cfg, list(id = "id", tree = sim$tree, log_lik = TRUE))
    code <- do.call(coev_make_stancode, args)
    code_centred <- with_mocked_bindings(
      do.call(coev_make_stancode, args),
      use_conditional_ncp = function(...) FALSE
    )
    expect_identical(code, code_centred)
    expect_false(grepl("ncp_terminal_drift", code, fixed = TRUE))
  }
})

test_that("non-centred Stan code adds the Jacobian and realised drift", {
  sim <- sim_mixed_data()
  code <- coev_make_stancode(
    data = sim$data,
    variables = list(x = "normal", z = "bernoulli_logit"),
    id = "id",
    tree = sim$tree,
    log_lik = TRUE
  )
  expect_match(code, "vector ncp_terminal_drift(", fixed = TRUE)
  expect_match(code, "target += ncp_log_det(", fixed = TRUE)
  # both model and generated quantities use the realised drift
  expect_match(code, "tdrift = ncp_terminal_drift(tdrift,", fixed = TRUE)
  gq <- sub(".*generated quantities\\{", "", code)
  expect_match(gq, "tdrifts = ncp_terminal_drift(tdrifts,", fixed = TRUE)
})

test_that("non-centred model is equivalent: normal + bernoulli", {
  skip_on_cran()
  sim <- sim_mixed_data()
  expect_ncp_equivalent(
    data = sim$data,
    variables = list(x = "normal", z = "bernoulli_logit"),
    id = "id",
    tree = sim$tree
  )
})

test_that("non-centred model is equivalent: missing data, non-normal first", {
  skip_on_cran()
  sim <- sim_mixed_data()
  d <- sim$data
  d$x[c(2, 5)] <- NA
  d$w[c(5, 7)] <- NA
  d$z[c(3, 7)] <- NA
  expect_ncp_equivalent(
    data = d,
    variables = list(z = "bernoulli_logit", x = "normal", w = "normal"),
    id = "id",
    tree = sim$tree
  )
})

test_that("non-centred model is equivalent: effects_mat, no correlated drift", {
  skip_on_cran()
  skip_if_not(run_extended_tests)
  sim <- sim_mixed_data()
  vars <- c("x", "w", "z")
  em <- matrix(TRUE, 3, 3, dimnames = list(vars, vars))
  em["x", "z"] <- em["z", "x"] <- FALSE
  expect_ncp_equivalent(
    data = sim$data,
    variables = list(x = "normal", w = "normal", z = "bernoulli_logit"),
    id = "id",
    tree = sim$tree,
    effects_mat = em,
    estimate_correlated_drift = FALSE
  )
})

test_that("non-centred model is equivalent: several response types", {
  skip_on_cran()
  skip_if_not(run_extended_tests)
  sim <- sim_mixed_data()
  d <- sim$data
  d$o[4] <- NA
  expect_ncp_equivalent(
    data = d,
    variables = list(
      o = "ordered_logistic",
      x = "normal",
      p = "poisson_softplus"
    ),
    id = "id",
    tree = sim$tree
  )
})

test_that("non-centred model is equivalent: measurement error", {
  skip_on_cran()
  skip_if_not(run_extended_tests)
  sim <- sim_mixed_data()
  d <- sim$data
  d$x[3] <- NA
  expect_ncp_equivalent(
    data = d,
    variables = list(x = "normal", z = "bernoulli_logit"),
    id = "id",
    tree = sim$tree,
    measurement_error = list(x = "x_se")
  )
})

test_that("non-centred model is equivalent: multiPhylo", {
  skip_on_cran()
  skip_if_not(run_extended_tests)
  sim <- sim_mixed_data()
  d <- sim$data
  d$x[2] <- NA
  tree2 <- withr::with_seed(2, ape::rcoal(12, tip.label = sim$tree$tip.label))
  trees <- c(sim$tree, tree2)
  expect_ncp_equivalent(
    data = d,
    variables = list(x = "normal", z = "bernoulli_logit"),
    id = "id",
    tree = trees
  )
})

test_that("non-centred model is equivalent: normal only with missing data", {
  skip_on_cran()
  skip_if_not(run_extended_tests)
  sim <- sim_mixed_data()
  d <- sim$data
  d$x[c(1, 4)] <- NA
  d$w[c(4, 9)] <- NA
  expect_ncp_equivalent(
    data = d,
    variables = list(x = "normal", w = "normal"),
    id = "id",
    tree = sim$tree
  )
})

test_that("non-centred model is equivalent: exact Gaussian process", {
  skip_on_cran()
  skip_if_not(run_extended_tests)
  sim <- sim_mixed_data()
  lon_lat <- data.frame(
    id = sim$data$id,
    longitude = sim$data$lon,
    latitude = sim$data$lat
  )
  expect_ncp_equivalent(
    data = sim$data,
    variables = list(x = "normal", z = "bernoulli_logit"),
    id = "id",
    tree = sim$tree,
    lon_lat = lon_lat
  )
})
