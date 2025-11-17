test_that("check_nutpie_available() detects nutpie installation", {
  # this test checks if we can detect nutpie availability
  # initially will fail until we implement check_nutpie_available()
  skip_if_not(coevolve:::check_nutpie_available(),
              message = "nutpie not available - skipping nutpie tests")
  # if we get here, nutpie should be available
  expect_true(coevolve:::check_nutpie_available())
})

test_that("nutpie can compile simple Stan model", {
  skip_if_not(coevolve:::check_nutpie_available(),
              message = "nutpie not available - skipping nutpie tests")
  # simple Stan model for testing
  stan_code <- "
  data {
    int<lower=0> N;
    vector[N] y;
  }
  parameters {
    real mu;
    real<lower=0> sigma;
  }
  model {
    mu ~ normal(0, 1);
    sigma ~ exponential(1);
    y ~ normal(mu, sigma);
  }
  "
  expect_no_error({
    compiled <- nutpie_compile_stan_model(stan_code)
  })
  expect_true(!is.null(compiled))
})

test_that("nutpie can sample from simple model", {
  skip_if_not(coevolve:::check_nutpie_available(),
              message = "nutpie not available - skipping nutpie tests")
  stan_code <- "
  data {
    int<lower=0> N;
    vector[N] y;
  }
  parameters {
    real mu;
    real<lower=0> sigma;
  }
  model {
    mu ~ normal(0, 1);
    sigma ~ exponential(1);
    y ~ normal(mu, sigma);
  }
  "
  withr::with_seed(1, {
    data_list <- list(N = 10L, y = rnorm(10))
  })
  expect_no_error({
    trace <- nutpie_sample(stan_code, data_list,
                           num_chains = 2L,
                           num_samples = 100L,
                           num_warmup = 50L,
                           seed = 12345L)
  })
  expect_true(!is.null(trace))
  # convert to draws_array
  draws <- convert_nutpie_draws(trace)
  # verify structure
  expect_s3_class(draws, "draws_array")
  expect_equal(posterior::ndraws(draws), 200L) # 2 chains * 100 samples
  expect_equal(posterior::nchains(draws), 2L)
  expect_true("mu" %in% posterior::variables(draws))
  expect_true("sigma" %in% posterior::variables(draws))
})

test_that("coev_fit() works with backend = 'nutpie'", {
  skip_if_not(coevolve:::check_nutpie_available(),
              message = "nutpie not available - skipping nutpie tests")
  # simple coevolutionary model
  withr::with_seed(1, {
    n <- 5
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = tree$tip.label,
      x = rnorm(n),
      y = rnorm(n)
    )
  })
  # fit with nutpie
  expect_no_error({
    fit <- coev_fit(
      data = d,
      variables = list(x = "normal", y = "normal"),
      id = "id",
      tree = tree,
      backend = "nutpie",
      chains = 2,
      iter_sampling = 100,
      iter_warmup = 50,
      seed = 12345,
      refresh = 0
    )
  })
  # verify coevfit structure
  expect_s3_class(fit, "coevfit")
  expect_true(!is.null(fit$fit))
  expect_true(!is.null(fit$stan_code))
  expect_true(!is.null(fit$stan_data))
  # summary method should work
  expect_no_error({
    s <- summary(fit)
  })
  expect_s3_class(s, "coevsummary")
  # summary.coevfit returns a list with elements like auto, cross, etc.
  # Each element is a data.frame with columns Estimate, Rhat, Bulk_ESS
  # Check that the summary structure is correct
  expect_true("auto" %in% names(s))
  if (nrow(s$auto) > 0) {
    expect_true("Estimate" %in% names(s$auto))
    expect_true("Rhat" %in% names(s$auto))
    expect_true("Bulk_ESS" %in% names(s$auto))
  }
  # extract samples method should work
  expect_no_error({
    samples <- extract_samples(fit)
  })
  expect_type(samples, "list")
  expect_true("A" %in% names(samples))
  expect_true("b" %in% names(samples))
  # plot method should work
  expect_no_error({
    p <- plot(fit, plot = FALSE)
  })
  expect_true(!is.null(p))
  # coev_ functions should run without error
  sw <- suppressWarnings
  expect_no_error(sw(coev_calculate_theta(fit, list(x = NA, y = 0))))
  expect_no_error(sw(coev_calculate_delta_theta(fit, "x", "y")))
  expect_no_error(sw(coev_plot_delta_theta(fit)))
  expect_no_error(sw(coev_plot_flowfield(fit, "x", "y")))
  expect_no_error(sw(coev_plot_pred_series(fit)))
  expect_no_error(sw(coev_plot_predictive_check(fit)))
  expect_no_error(sw(coev_plot_selection_gradient(fit, "x", "y")))
  expect_no_error(sw(coev_plot_trait_values(fit)))
  expect_no_error(sw(coev_pred_series(fit)))
})

test_that("nutpie and cmdstanr produce similar results", {
  skip_if_not(coevolve:::check_nutpie_available(),
              message = "nutpie not available - skipping nutpie tests")
  withr::with_seed(1, {
    n <- 5
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = tree$tip.label,
      x = rnorm(n),
      y = rnorm(n)
    )
  })
  # fit with both samplers using same seed
  fit_cmdstanr <- coev_fit(
    data = d,
    variables = list(x = "normal", y = "normal"),
    id = "id",
    tree = tree,
    backend = "cmdstanr",
    chains = 2,
    iter_sampling = 100,
    iter_warmup = 50,
    seed = 12345,
    refresh = 0
  )
  fit_nutpie <- coev_fit(
    data = d,
    variables = list(x = "normal", y = "normal"),
    id = "id",
    tree = tree,
    backend = "nutpie",
    chains = 2,
    iter_sampling = 100,
    iter_warmup = 50,
    seed = 12345,
    refresh = 0
  )
  # extract posterior means for key parameters
  samples_cmdstanr <- extract_samples(fit_cmdstanr)
  samples_nutpie <- extract_samples(fit_nutpie)
  # compare means (should be similar, not identical)
  # note: different samplers will produce different samples
  # we're just checking that both produce reasonable results
  expect_true(is.numeric(samples_cmdstanr$A))
  expect_true(is.numeric(samples_nutpie$A))
  # both should have similar means (within reasonable tolerance)
  mean_cmdstanr <- mean(samples_cmdstanr$A)
  mean_nutpie <- mean(samples_nutpie$A)
  # allow for some difference due to different samplers
  expect_lt(abs(mean_cmdstanr - mean_nutpie), 1.0)
})

test_that("nutpie handles errors gracefully", {
  skip_if_not(coevolve:::check_nutpie_available(),
              message = "nutpie not available - skipping nutpie tests")
  # invalid Stan code should produce informative error
  invalid_code <- "invalid stan code here"
  expect_error(
    nutpie_compile_stan_model(invalid_code),
    regexp = ".*"
  )
  # invalid data should produce informative error
  stan_code <- "
  data {
    int<lower=0> N;
    vector[N] y;
  }
  parameters {
    real mu;
  }
  model {
    mu ~ normal(0, 1);
    y ~ normal(mu, 1);
  }
  "
  invalid_data <- list(N = 10L, y = "not a vector")
  expect_error(
    nutpie_sample(stan_code, invalid_data,
                  num_chains = 1L,
                  num_samples = 10L,
                  num_warmup = 5L),
    regexp = ".*"
  )
})

test_that("coev_fit() converts parallel_chains to cores for nutpie", {
  skip_if_not(coevolve:::check_nutpie_available(),
              message = "nutpie not available - skipping nutpie tests")
  withr::with_seed(1, {
    n <- 5
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = tree$tip.label,
      x = rnorm(n),
      y = rnorm(n)
    )
  })
  # Test that parallel_chains argument is accepted and converted to cores
  # If conversion fails, nutpie would error on unknown argument
  expect_no_error({
    fit <- coev_fit(
      data = d,
      variables = list(x = "normal", y = "normal"),
      id = "id",
      tree = tree,
      backend = "nutpie",
      chains = 2,
      iter_sampling = 50,
      iter_warmup = 25,
      seed = 12345,
      parallel_chains = 2,  # Should be converted to cores=2 for nutpie
      refresh = 0
    )
  })
  # Verify fit succeeded
  expect_s3_class(fit, "coevfit")
  expect_true(!is.null(fit$fit))
})

test_that("coev_fit() errors when backend = 'nutpie' but nutpie unavailable", {
  # mock check_nutpie_available to return FALSE
  # this test ensures proper error message when nutpie not available
  # note: this test may need to be adjusted based on implementation
  # for now, we expect an informative error
  withr::with_seed(1, {
    n <- 3
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = tree$tip.label,
      x = rnorm(n),
      y = rnorm(n)
    )
  })
  # if nutpie is actually available, skip this test
  skip_if(coevolve:::check_nutpie_available(),
          message = "nutpie is available - skipping unavailable test")
  expect_error(
    coev_fit(
      data = d,
      variables = list(x = "normal", y = "normal"),
      id = "id",
      tree = tree,
      backend = "nutpie",
      chains = 1,
      iter_sampling = 10,
      iter_warmup = 5,
      seed = 1,
      refresh = 0
    ),
    regexp = "nutpie.*not.*available|nutpie.*not.*installed"
  )
})
