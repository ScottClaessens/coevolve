# Tests for nutpie integration with coevolve
#
# These tests follow test-driven development principles:
# 1. Write tests first (they will initially fail)
# 2. Implement minimal code to make tests pass
# 3. Refactor and expand tests

test_that("check_nutpie_available() detects nutpie installation", {
  # This test checks if we can detect nutpie availability
  # Initially will fail until we implement check_nutpie_available()
  
  skip_if_not(check_nutpie_available(), 
              message = "nutpie not available - skipping nutpie tests")
  
  # If we get here, nutpie should be available
  expect_true(check_nutpie_available())
})

test_that("nutpie can compile simple Stan model", {
  skip_if_not(check_nutpie_available(), 
              message = "nutpie not available - skipping nutpie tests")
  
  # Simple Stan model for testing
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
  
  # This will fail until we implement nutpie compilation
  expect_no_error({
    compiled <- nutpie_compile_stan_model(stan_code)
  })
  
  expect_true(!is.null(compiled))
})

test_that("nutpie can sample from simple model", {
  skip_if_not(check_nutpie_available(), 
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
  
  data_list <- list(N = 10L, y = rnorm(10))
  
  # This will fail until we implement nutpie sampling
  expect_no_error({
    trace <- nutpie_sample(stan_code, data_list, 
                          num_chains = 2L, 
                          num_samples = 100L,
                          num_warmup = 50L,
                          seed = 12345L)
  })
  
  expect_true(!is.null(trace))
})

test_that("nutpie draws can be converted to draws_array", {
  skip_if_not(check_nutpie_available(), 
              message = "nutpie not available - skipping nutpie tests")
  
  # Fit a simple model with nutpie
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
  
  data_list <- list(N = 10L, y = rnorm(10))
  
  trace <- nutpie_sample(stan_code, data_list,
                         num_chains = 2L,
                         num_samples = 100L,
                         num_warmup = 50L,
                         seed = 12345L)
  
  # Convert to draws_array
  draws <- convert_nutpie_draws(trace)
  
  # Verify structure
  expect_s3_class(draws, "draws_array")
  expect_equal(posterior::ndraws(draws), 200L) # 2 chains * 100 samples
  expect_equal(posterior::nchains(draws), 2L)
  expect_true("mu" %in% posterior::variables(draws))
  expect_true("sigma" %in% posterior::variables(draws))
})

test_that("coev_fit() works with sampler = 'nutpie'", {
  skip_if_not(check_nutpie_available(), 
              message = "nutpie not available - skipping nutpie tests")
  
  # Simple coevolutionary model
  withr::with_seed(1, {
    n <- 5
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = tree$tip.label,
      x = rnorm(n),
      y = rnorm(n)
    )
  })
  
  # Fit with nutpie
  expect_no_error({
    fit <- coev_fit(
      data = d,
      variables = list(x = "normal", y = "normal"),
      id = "id",
      tree = tree,
      sampler = "nutpie",
      chains = 2,
      iter_sampling = 100,
      iter_warmup = 50,
      seed = 12345,
      refresh = 0
    )
  })
  
  # Verify coevfit structure
  expect_s3_class(fit, "coevfit")
  expect_true(!is.null(fit$fit))
  expect_true(!is.null(fit$stan_code))
  expect_true(!is.null(fit$stan_data))
})

test_that("summary() works with nutpie-fitted models", {
  skip_if_not(check_nutpie_available(), 
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
  
  fit <- coev_fit(
    data = d,
    variables = list(x = "normal", y = "normal"),
    id = "id",
    tree = tree,
    sampler = "nutpie",
    chains = 2,
    iter_sampling = 100,
    iter_warmup = 50,
    seed = 12345,
    refresh = 0
  )
  
  # Summary should work
  expect_no_error({
    s <- summary(fit)
  })
  
  expect_s3_class(s, "coevsummary")
  expect_true("Estimate" %in% names(s))
  expect_true("Rhat" %in% names(s))
  expect_true("Bulk_ESS" %in% names(s))
})

test_that("extract_samples() works with nutpie-fitted models", {
  skip_if_not(check_nutpie_available(), 
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
  
  fit <- coev_fit(
    data = d,
    variables = list(x = "normal", y = "normal"),
    id = "id",
    tree = tree,
    sampler = "nutpie",
    chains = 2,
    iter_sampling = 100,
    iter_warmup = 50,
    seed = 12345,
    refresh = 0
  )
  
  # Extract samples should work
  expect_no_error({
    samples <- extract_samples(fit)
  })
  
  expect_type(samples, "list")
  expect_true("A" %in% names(samples))
  expect_true("b" %in% names(samples))
})

test_that("plot() works with nutpie-fitted models", {
  skip_if_not(check_nutpie_available(), 
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
  
  fit <- coev_fit(
    data = d,
    variables = list(x = "normal", y = "normal"),
    id = "id",
    tree = tree,
    sampler = "nutpie",
    chains = 2,
    iter_sampling = 100,
    iter_warmup = 50,
    seed = 12345,
    refresh = 0
  )
  
  # Plot should work
  expect_no_error({
    p <- plot(fit, plot = FALSE)
  })
  
  expect_true(!is.null(p))
})

test_that("nutpie and cmdstanr produce similar results", {
  skip_if_not(check_nutpie_available(), 
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
  
  # Fit with both samplers using same seed
  fit_cmdstanr <- coev_fit(
    data = d,
    variables = list(x = "normal", y = "normal"),
    id = "id",
    tree = tree,
    sampler = "cmdstanr",
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
    sampler = "nutpie",
    chains = 2,
    iter_sampling = 100,
    iter_warmup = 50,
    seed = 12345,
    refresh = 0
  )
  
  # Extract posterior means for key parameters
  samples_cmdstanr <- extract_samples(fit_cmdstanr)
  samples_nutpie <- extract_samples(fit_nutpie)
  
  # Compare means (should be similar, not identical)
  # Note: Different samplers will produce different samples
  # We're just checking that both produce reasonable results
  expect_true(is.numeric(samples_cmdstanr$A))
  expect_true(is.numeric(samples_nutpie$A))
  
  # Both should have similar means (within reasonable tolerance)
  mean_cmdstanr <- mean(samples_cmdstanr$A)
  mean_nutpie <- mean(samples_nutpie$A)
  
  # Allow for some difference due to different samplers
  expect_lt(abs(mean_cmdstanr - mean_nutpie), 1.0)
})

test_that("nutpie handles errors gracefully", {
  skip_if_not(check_nutpie_available(), 
              message = "nutpie not available - skipping nutpie tests")
  
  # Invalid Stan code should produce informative error
  invalid_code <- "invalid stan code here"
  
  expect_error(
    nutpie_compile_stan_model(invalid_code),
    regexp = ".*"
  )
  
  # Invalid data should produce informative error
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

test_that("coev_fit() defaults to cmdstanr when sampler not specified", {
  # Should work without specifying sampler (backward compatibility)
  withr::with_seed(1, {
    n <- 3
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = tree$tip.label,
      x = rnorm(n),
      y = rnorm(n)
    )
  })
  
  expect_no_error({
    fit <- coev_fit(
      data = d,
      variables = list(x = "normal", y = "normal"),
      id = "id",
      tree = tree,
      chains = 1,
      iter_sampling = 10,
      iter_warmup = 5,
      seed = 1,
      refresh = 0
    )
  })
  
  expect_s3_class(fit, "coevfit")
})

test_that("coev_fit() errors when sampler = 'nutpie' but nutpie unavailable", {
  # Mock check_nutpie_available to return FALSE
  # This test ensures proper error message when nutpie not available
  
  # Note: This test may need to be adjusted based on implementation
  # For now, we expect an informative error
  withr::with_seed(1, {
    n <- 3
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = tree$tip.label,
      x = rnorm(n),
      y = rnorm(n)
    )
  })
  
  # If nutpie is actually available, skip this test
  skip_if(check_nutpie_available(),
          message = "nutpie is available - skipping unavailable test")
  
  expect_error(
    coev_fit(
      data = d,
      variables = list(x = "normal", y = "normal"),
      id = "id",
      tree = tree,
      sampler = "nutpie",
      chains = 1,
      iter_sampling = 10,
      iter_warmup = 5,
      seed = 1,
      refresh = 0
    ),
    regexp = "nutpie.*not.*available|nutpie.*not.*installed"
  )
})

