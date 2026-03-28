test_that("check_jax_available() returns logical", {
  result <- coevolve:::check_jax_available()
  expect_type(result, "logical")
  expect_length(result, 1L)
})

test_that("normalize_nuts_sampler accepts stan and nutpie", {
  expect_equal(
    coevolve:::normalize_nuts_sampler("stan"), "stan"
  )
  expect_equal(
    coevolve:::normalize_nuts_sampler("nutpie"), "nutpie"
  )
  expect_equal(
    coevolve:::normalize_nuts_sampler("NUTPIE"), "nutpie"
  )
  expect_equal(
    coevolve:::normalize_nuts_sampler("Stan"), "stan"
  )
})

test_that("normalize_nuts_sampler rejects invalid values", {
  expect_error(
    coevolve:::normalize_nuts_sampler("invalid"),
    "one of: stan, nutpie",
    fixed = TRUE
  )
  expect_error(
    coevolve:::normalize_nuts_sampler("pymc"),
    "one of: stan, nutpie",
    fixed = TRUE
  )
  expect_error(
    coevolve:::normalize_nuts_sampler(123),
    "Argument 'nuts_sampler' must be a single",
    fixed = TRUE
  )
  expect_error(
    coevolve:::normalize_nuts_sampler(NA_character_),
    "Argument 'nuts_sampler' must be a single",
    fixed = TRUE
  )
  expect_error(
    coevolve:::normalize_nuts_sampler(c("stan", "nutpie")),
    "Argument 'nuts_sampler' must be a single",
    fixed = TRUE
  )
})

test_that("coev_fit errors with nuts_sampler = 'invalid'", {
  withr::with_seed(1, {
    n <- 3
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = tree$tip.label,
      x = rnorm(n),
      y = rnorm(n)
    )
  })
  expect_error(
    coev_fit(
      data = d,
      variables = list(x = "normal", y = "normal"),
      id = "id",
      tree = tree,
      nuts_sampler = "invalid"
    ),
    "one of: stan, nutpie",
    fixed = TRUE
  )
})

test_that("stop_if_jax_not_available gives informative error", {
  # temporarily override check_jax_available in namespace
  ns <- asNamespace("coevolve")
  orig <- ns$check_jax_available
  unlockBinding("check_jax_available", ns)
  on.exit({
    assign("check_jax_available", orig, envir = ns)
    lockBinding("check_jax_available", ns)
  })
  assign(
    "check_jax_available",
    function() FALSE,
    envir = ns
  )
  expect_error(
    coevolve:::stop_if_jax_not_available(),
    "JAX backend is not available",
    fixed = TRUE
  )
})

test_that("coev_fit with nuts_sampler='jax' returns coevfit", {
  skip_if_not(
    coevolve:::check_jax_available(),
    message = "JAX not available - skipping JAX tests"
  )
  fit <- coev_fit(
    data = authority$data,
    variables = list(
      political_authority = "ordered_logistic",
      religious_authority = "ordered_logistic"
    ),
    id = "language",
    tree = authority$phylogeny,
    nuts_sampler = "nutpie",
    chains = 1,
    iter_warmup = 50,
    iter_sampling = 50,
    seed = 1
  )
  expect_s3_class(fit, "coevfit")
  expect_true(!is.null(fit$fit))
  expect_true(!is.null(fit$stan_data))
  expect_equal(fit$nuts_sampler, "nutpie")
  expect_equal(fit$backend, "nutpie")
})

test_that("jax_fit class is set on fit$fit", {
  skip_if_not(
    coevolve:::check_jax_available(),
    message = "JAX not available - skipping JAX tests"
  )
  fit <- coev_fit(
    data = authority$data,
    variables = list(
      political_authority = "ordered_logistic",
      religious_authority = "ordered_logistic"
    ),
    id = "language",
    tree = authority$phylogeny,
    nuts_sampler = "nutpie",
    chains = 1,
    iter_warmup = 50,
    iter_sampling = 50,
    seed = 1
  )
  expect_s3_class(fit$fit, "jax_fit")
})

test_that("draws extraction from JAX fit works", {
  skip_if_not(
    coevolve:::check_jax_available(),
    message = "JAX not available - skipping JAX tests"
  )
  fit <- coev_fit(
    data = authority$data,
    variables = list(
      political_authority = "ordered_logistic",
      religious_authority = "ordered_logistic"
    ),
    id = "language",
    tree = authority$phylogeny,
    nuts_sampler = "nutpie",
    chains = 1,
    iter_warmup = 50,
    iter_sampling = 50,
    seed = 1
  )
  draws <- fit$fit$draws()
  expect_s3_class(draws, "draws_array")
  expect_equal(posterior::nchains(draws), 1L)
  expect_equal(posterior::niterations(draws), 50L)
  # key parameters should be present
  vars <- posterior::variables(draws)
  expect_true(any(grepl("^A\\[", vars)))
  expect_true(any(grepl("^b\\[", vars)))
})

test_that("summary of JAX fit prints without error", {
  skip_if_not(
    coevolve:::check_jax_available(),
    message = "JAX not available - skipping JAX tests"
  )
  fit <- coev_fit(
    data = authority$data,
    variables = list(
      political_authority = "ordered_logistic",
      religious_authority = "ordered_logistic"
    ),
    id = "language",
    tree = authority$phylogeny,
    nuts_sampler = "nutpie",
    chains = 1,
    iter_warmup = 50,
    iter_sampling = 50,
    seed = 1
  )
  sw <- suppressWarnings
  expect_no_error(sw(summary(fit)))
  expect_no_error(sw(print(fit)))
  expect_no_error(sw(print(summary(fit))))
})

test_that("JAX fit respects sampling arguments", {
  skip_if_not(
    coevolve:::check_jax_available(),
    message = "JAX not available - skipping JAX tests"
  )
  n_chains <- 2L
  n_warmup <- 30L
  n_sampling <- 40L
  seed_val <- 42L
  fit <- coev_fit(
    data = authority$data,
    variables = list(
      political_authority = "ordered_logistic",
      religious_authority = "ordered_logistic"
    ),
    id = "language",
    tree = authority$phylogeny,
    nuts_sampler = "nutpie",
    chains = n_chains,
    iter_warmup = n_warmup,
    iter_sampling = n_sampling,
    seed = seed_val
  )
  wrapper <- fit$fit
  expect_equal(wrapper$chains, n_chains)
  expect_equal(wrapper$iter_warmup, n_warmup)
  expect_equal(wrapper$iter_sampling, n_sampling)
  expect_equal(wrapper$seed, seed_val)
  draws <- wrapper$draws()
  expect_equal(posterior::nchains(draws), n_chains)
  expect_equal(
    posterior::niterations(draws), n_sampling
  )
})

test_that("JAX fit metadata method works", {
  skip_if_not(
    coevolve:::check_jax_available(),
    message = "JAX not available - skipping JAX tests"
  )
  fit <- coev_fit(
    data = authority$data,
    variables = list(
      political_authority = "ordered_logistic",
      religious_authority = "ordered_logistic"
    ),
    id = "language",
    tree = authority$phylogeny,
    nuts_sampler = "nutpie",
    chains = 1,
    iter_warmup = 50,
    iter_sampling = 50,
    seed = 1
  )
  meta <- fit$fit$metadata()
  expect_type(meta, "list")
  expect_equal(meta$num_chains, 1L)
  expect_equal(meta$iter_sampling, 50L)
  expect_equal(meta$iter_warmup, 50L)
  expect_true(length(meta$stan_variables) > 0)
})

test_that("JAX fit extract_samples works", {
  skip_if_not(
    coevolve:::check_jax_available(),
    message = "JAX not available - skipping JAX tests"
  )
  fit <- coev_fit(
    data = authority$data,
    variables = list(
      political_authority = "ordered_logistic",
      religious_authority = "ordered_logistic"
    ),
    id = "language",
    tree = authority$phylogeny,
    nuts_sampler = "nutpie",
    chains = 1,
    iter_warmup = 50,
    iter_sampling = 50,
    seed = 1
  )
  sw <- suppressWarnings
  expect_no_error(sw(samples <- extract_samples(fit)))
  expect_type(samples, "list")
  expect_true("A" %in% names(samples))
  expect_true("b" %in% names(samples))
})

test_that("JAX fit works with simple normal model", {
  skip_if_not(
    coevolve:::check_jax_available(),
    message = "JAX not available - skipping JAX tests"
  )
  withr::with_seed(1, {
    n <- 5
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = tree$tip.label,
      x = rnorm(n),
      y = rnorm(n)
    )
  })
  fit <- coev_fit(
    data = d,
    variables = list(
      x = "normal",
      y = "normal"
    ),
    id = "id",
    tree = tree,
    nuts_sampler = "nutpie",
    chains = 1,
    iter_warmup = 50,
    iter_sampling = 50,
    seed = 1
  )
  expect_s3_class(fit, "coevfit")
  expect_s3_class(fit$fit, "jax_fit")
  sw <- suppressWarnings
  expect_no_error(sw(summary(fit)))
})

test_that("JAX fit with exact GP returns coevfit with GP params", {
  skip_if_not(
    coevolve:::check_jax_available(),
    message = "JAX not available - skipping JAX tests"
  )
  fit <- coev_fit(
    data = authority$data,
    variables = list(
      political_authority = "ordered_logistic",
      religious_authority = "ordered_logistic"
    ),
    id = "language",
    tree = authority$phylogeny,
    lon_lat = authority$coordinates,
    prior = list(A_offdiag = "normal(0, 2)"),
    nuts_sampler = "nutpie",
    chains = 1,
    iter_warmup = 50,
    iter_sampling = 50,
    seed = 1
  )
  expect_s3_class(fit, "coevfit")
  expect_s3_class(fit$fit, "jax_fit")
  vars <- posterior::variables(fit$fit$draws())
  expect_true(any(grepl("^sigma_dist", vars)))
  expect_true(any(grepl("^rho_dist", vars)))
})

test_that("JAX fit with HSGP returns coevfit with GP params", {
  skip_if_not(
    coevolve:::check_jax_available(),
    message = "JAX not available - skipping JAX tests"
  )
  fit <- coev_fit(
    data = authority$data,
    variables = list(
      political_authority = "ordered_logistic",
      religious_authority = "ordered_logistic"
    ),
    id = "language",
    tree = authority$phylogeny,
    lon_lat = authority$coordinates,
    dist_k = 3,
    prior = list(A_offdiag = "normal(0, 2)"),
    nuts_sampler = "nutpie",
    chains = 1,
    iter_warmup = 50,
    iter_sampling = 50,
    seed = 1
  )
  expect_s3_class(fit, "coevfit")
  expect_s3_class(fit$fit, "jax_fit")
  vars <- posterior::variables(fit$fit$draws())
  expect_true(any(grepl("^sigma_dist", vars)))
  expect_true(any(grepl("^rho_dist", vars)))
})
