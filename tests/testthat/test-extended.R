#' @srrstats {G5.10} Flag extended tests
run_extended_tests <- identical(Sys.getenv("COEVOLVE_EXTENDED_TESTS"), "true")

test_that("coev_fit() fits test fixture models without error", {
  skip_if_not(run_extended_tests)
  withr::with_seed(1234, {
    # simulate data
    n <- 5
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = tree$tip.label,
      u = rgamma(n, shape = 1, rate = 1),
      v = rnorm(n),
      w = rbinom(n, size = 1, prob = 0.5),
      x = ordered(sample(1:4, size = n, replace = TRUE)),
      y = as.integer(rnbinom(n, mu = 3, size = 1)),
      z = as.integer(rnbinom(n, mu = 3, size = 1))
    )
    # stan arguments
    warmup <- 50
    iter <- 50
    chains <- 1
    # fit model without distance matrix
    coevfit_example_01 <-
      coev_fit(
        data = d,
        variables = list(
          u = "gamma_log",
          v = "normal",
          w = "bernoulli_logit",
          x = "ordered_logistic",
          y = "poisson_softplus"
        ),
        id = "id",
        tree = tree,
        chains = chains,
        iter_warmup = warmup,
        iter_sampling = iter,
        adapt_delta = 0.99,
        seed = 12345
      )
    # fit model with distance matrix
    dist_mat <- as.matrix(dist(rnorm(n)))
    rownames(dist_mat) <- colnames(dist_mat) <- d$id
    coevfit_example_02 <-
      coev_fit(
        data = d,
        variables = list(
          w = "bernoulli_logit",
          x = "ordered_logistic"
        ),
        id = "id",
        tree = tree,
        dist_mat = dist_mat,
        chains = chains,
        iter_warmup = warmup,
        iter_sampling = iter,
        adapt_delta = 0.99,
        seed = 12345
      )
    # fit prior only model
    coevfit_example_03 <-
      coev_fit(
        data = d,
        variables = list(
          w = "bernoulli_logit",
          x = "ordered_logistic"
        ),
        id = "id",
        tree = tree,
        chains = chains,
        iter_warmup = warmup,
        iter_sampling = iter,
        adapt_delta = 0.99,
        seed = 12345,
        prior_only = TRUE
      )
    # fit negative binomial model
    coevfit_example_04 <-
      coev_fit(
        data = d,
        variables = list(
          y = "poisson_softplus",
          z = "negative_binomial_softplus"
        ),
        id = "id",
        tree = tree,
        chains = chains,
        iter_warmup = warmup,
        iter_sampling = iter,
        adapt_delta = 0.99,
        seed = 12345
      )
    # fit model with effects matrix
    effects_mat <- matrix(
      c(TRUE, TRUE,
        FALSE, TRUE),
      byrow = TRUE,
      nrow = 2,
      ncol = 2,
      dimnames = list(c("w", "x"), c("w", "x"))
    )
    coevfit_example_05 <-
      coev_fit(
        data = d,
        variables = list(
          w = "bernoulli_logit",
          x = "ordered_logistic"
        ),
        id = "id",
        tree = tree,
        effects_mat = effects_mat,
        chains = chains,
        iter_warmup = warmup,
        iter_sampling = iter,
        adapt_delta = 0.99,
        seed = 12345
      )
    # fit model with missing data
    d$w[c(1, 2)] <- NA
    d$x[c(1, 3)] <- NA
    coevfit_example_06 <-
      coev_fit(
        data = d,
        variables = list(
          w = "bernoulli_logit",
          x = "ordered_logistic"
        ),
        id = "id",
        tree = tree,
        chains = chains,
        iter_warmup = warmup,
        iter_sampling = iter,
        adapt_delta = 0.99,
        seed = 12345
      )
    # fit model with repeated observations
    d <- data.frame(
      id = rep(tree$tip.label, each = 3),
      w = rbinom(n * 3, size = 1, prob = 0.5),
      x = ordered(sample(1:4, size = n * 3, replace = TRUE))
    )
    coevfit_example_07 <-
      coev_fit(
        data = d,
        variables = list(
          w = "bernoulli_logit",
          x = "ordered_logistic"
        ),
        id = "id",
        tree = tree,
        chains = chains,
        iter_warmup = warmup,
        iter_sampling = iter,
        adapt_delta = 0.99,
        seed = 12345
      )
    # fit model with multiPhylo object
    tree <- c(tree, ape::rcoal(n))
    d <- data.frame(
      id = tree[[1]]$tip.label,
      x = rbinom(n, 1, 0.5),
      y = rbinom(n, 1, 0.5)
    )
    coevfit_example_08 <-
      coev_fit(
        data = d,
        variables = list(
          x = "bernoulli_logit",
          y = "bernoulli_logit"
        ),
        id = "id",
        tree = tree,
        chains = chains,
        iter_warmup = warmup,
        iter_sampling = iter,
        adapt_delta = 0.99,
        seed = 12345
      )
    # set Q off diagonals to zero
    tree <- tree[[1]]
    d <- data.frame(
      id = tree$tip.label,
      x = rbinom(n, 1, 0.5),
      y = rbinom(n, 1, 0.5)
    )
    coevfit_example_09 <-
      coev_fit(
        data = d,
        variables = list(
          x = "bernoulli_logit",
          y = "bernoulli_logit"
        ),
        id = "id",
        tree = tree,
        estimate_correlated_drift = FALSE,
        chains = chains,
        iter_warmup = warmup,
        iter_sampling = iter,
        adapt_delta = 0.99,
        seed = 12345
      )
    # fit model with measurement error
    d <- data.frame(
      id = tree$tip.label,
      x = rnorm(n),
      y = rnorm(n),
      x_se = rexp(n, 5),
      y_se = rexp(n, 5)
    )
    coevfit_example_10 <-
      coev_fit(
        data = d,
        variables = list(
          x = "normal",
          y = "normal"
        ),
        id = "id",
        tree = tree,
        measurement_error = list(
          x = "x_se",
          y = "y_se"
        ),
        chains = chains,
        iter_warmup = warmup,
        iter_sampling = iter,
        adapt_delta = 0.99,
        seed = 12345
      )
  })
  # did all models fit without error?
  sw <- suppressWarnings
  expect_no_error(sw(coevfit_example_01))
  expect_no_error(sw(coevfit_example_02))
  expect_no_error(sw(coevfit_example_03))
  expect_no_error(sw(coevfit_example_04))
  expect_no_error(sw(coevfit_example_05))
  expect_no_error(sw(coevfit_example_06))
  expect_no_error(sw(coevfit_example_07))
  expect_no_error(sw(coevfit_example_08))
  expect_no_error(sw(coevfit_example_09))
  expect_no_error(sw(coevfit_example_10))
})

#' @srrstats {G5.4, G5.4a, G5.5} Extended correctness test with fixed seed
test_that("coev_fit() estimates correct direction of selection", {
  skip_if_not(run_extended_tests)
  # simulate data with fixed seed, where x -> y but not vice versa
  withr::with_seed(1, {
    sim <-
      coev_simulate_coevolution(
        n = 100,
        variables = c("x", "y"),
        selection_matrix = matrix(
          c(0.9, 0.0,
            0.9, 0.9),
          nrow = 2,
          byrow = TRUE,
          dimnames = list(c("x", "y"), c("x", "y"))
        ),
        drift = c("x" = 0.01, "y" = 0.01),
        prob_split = 0.05
      )
  })
  # fit model
  fit <-
    coev_fit(
      data = sim$data,
      variables = list(
        x = "normal",
        y = "normal"
      ),
      id = "species",
      tree = sim$tree,
      estimate_correlated_drift = FALSE,
      chains = 1,
      seed = 1,
      refresh = 0
    )
  # check model converged
  expect_lt(max(suppressWarnings(fit$fit$summary())$rhat, na.rm = TRUE), 1.1)
  # get delta theta values
  delta_theta_x_to_y <-
    coev_calculate_delta_theta(
      object = fit,
      predictor = "x",
      response = "y"
    )
  delta_theta_y_to_x <-
    coev_calculate_delta_theta(
      object = fit,
      predictor = "y",
      response = "x"
    )
  # for x -> y, 95% CI for delta theta should be greater than zero
  expect_true(
    quantile(delta_theta_x_to_y, 0.025) > 0 &&
      quantile(delta_theta_x_to_y, 0.975) > 0
  )
  # for y -> x, 95% CI for delta theta should include zero
  expect_true(
    quantile(delta_theta_y_to_x, 0.025) < 0 &&
      quantile(delta_theta_y_to_x, 0.975) > 0
  )
})

#' @srrstats {G5.6, G5.6a, G5.6b, G5.9, G5.9a, G5.9b, BS7.2} Extended parameter
#'   recovery tests with multiple fixed seeds for data simulation and cmdstanr
#' @srrstats {BS7.4, BS7.4a} Predicted values are on same scale as input data
for (seed in 1:3) {
  test_that(paste0("coev_fit() recovers parameters (seed = ", seed, ")"), {
    skip_if_not(run_extended_tests)
    # get dummy data
    withr::with_seed(seed, {
      n <- 100
      tree <- ape::rcoal(n)
      d <- data.frame(
        id = tree$tip.label,
        x = rnorm(n),
        y = rnorm(n)
      )
    })
    # get stan data list with prior_only
    sdata <- coev_make_standata(
      data = d,
      variables = list(
        x = "normal",
        y = "normal"
      ),
      id = "id",
      tree = tree,
      estimate_correlated_drift = FALSE,
      prior_only = TRUE
    )
    # manually fix parameters by editing the stan code
    scode <- coev_make_stancode(
      data = d,
      variables = list(
        x = "normal",
        y = "normal"
      ),
      id = "id",
      tree = tree,
      estimate_correlated_drift = FALSE
    )
    scode_fixed <- manually_fix_parameters(scode)
    # simulate dataset from model with fixed parameters
    sim <-
      suppressWarnings(
        cmdstanr::cmdstan_model(
          stan_file = cmdstanr::write_stan_file(scode_fixed)
        )$sample(
          data = sdata,
          chains = 1,
          refresh = 0,
          seed = seed,
          iter_warmup = 50,
          iter_sampling = 1
        )
      )
    draws <- posterior::as_draws_rvars(sim)
    d_sim <- data.frame(
      id = tree$tip.label,
      x = posterior::draws_of(draws$yrep)[1, 1, 1:n, 1],
      y = posterior::draws_of(draws$yrep)[1, 1, 1:n, 2]
    )
    # fit model to simulated data
    effects_mat <- matrix(TRUE, nrow = 2, ncol = 2,
                          dimnames = list(c("x", "y"), c("x", "y")))
    effects_mat[1, 2] <- FALSE
    fit <-
      suppressWarnings(
        coev_fit(
          data = d_sim,
          variables = list(
            x = "normal",
            y = "normal"
          ),
          id = "id",
          tree = tree,
          estimate_correlated_drift = FALSE,
          effects_mat = effects_mat,
          prior = list(
            A_offdiag = "normal(0, 2)",
            Q_sigma = "normal(0, 2)"
          ),
          scale = FALSE,
          chains = 1,
          refresh = 0,
          seed = seed,
          max_treedepth = 15
        )
      )
    # check model converged
    expect_lt(max(suppressWarnings(fit$fit$summary())$rhat, na.rm = TRUE), 1.1)
    # posterior medians should recover fixed parameters within tolerance
    post <- extract_samples(fit)
    expect_equal(median(post$A[, 1, 1]), -0.5, tolerance = 1)
    expect_equal(median(post$A[, 2, 2]), -0.5, tolerance = 1)
    expect_equal(median(post$A[, 2, 1]),  1.0, tolerance = 1)
    expect_equal(median(post$Q[, 1, 1]),  1.5, tolerance = 1)
    expect_equal(median(post$Q[, 2, 2]),  1.5, tolerance = 1)
    expect_equal(median(post$b[, 1]), 0.0, tolerance = 1)
    expect_equal(median(post$b[, 2]), 0.0, tolerance = 1)
    # median predicted values on same scale as input data
    pred <- apply(post$yrep, c(3, 4), median)
    expect_equal(median(d_sim$x), median(pred[, 1]), tolerance = 0.5)
    expect_equal(median(d_sim$y), median(pred[, 2]), tolerance = 0.5)
    expect_equal(sd(d_sim$x), sd(pred[, 1]), tolerance = 0.5)
    expect_equal(sd(d_sim$y), sd(pred[, 2]), tolerance = 0.5)
  })
}

#' @srrstats {G5.7} Standard errors for posterior estimates are smaller as
#'   number of observations increases
#' @srrstats {BS7.3} Algorithm scaling test
test_that("coev_fit() scales with increasing observations", {
  skip_if_not(run_extended_tests)
  # get dataset with repeated measures
  withr::with_seed(1, {
    n <- 10
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = tree$tip.label,
      x = rnorm(n),
      y = rnorm(n)
    )
    d <- rbind(d, d, d)
  })
  # fit model
  fit_small <-
    coev_fit(
      data = d,
      variables = list(
        x = "normal",
        y = "normal"
      ),
      id = "id",
      tree = tree,
      estimate_residual = FALSE,
      chains = 1,
      refresh = 0,
      seed = 1
    )
  s_small <- fit_small$fit$summary()
  # fit model with 5x the number of observations
  fit_large <-
    coev_fit(
      data = rbind(d, d, d, d, d),
      variables = list(
        x = "normal",
        y = "normal"
      ),
      id = "id",
      tree = tree,
      estimate_residual = FALSE,
      chains = 1,
      refresh = 0,
      seed = 1
    )
  s_large <- fit_large$fit$summary()
  # check models converged
  expect_lt(
    max(suppressWarnings(fit_small$fit$summary())$rhat, na.rm = TRUE),
    1.1
  )
  expect_lt(
    max(suppressWarnings(fit_large$fit$summary())$rhat, na.rm = TRUE),
    1.1
  )
  # standard errors should be smaller with more observations
  # autoregressive effects
  expect_lt(as.numeric(s_large[1, "sd"]), as.numeric(s_small[1, "sd"]))
  expect_lt(as.numeric(s_large[2, "sd"]), as.numeric(s_small[2, "sd"]))
  # cross effects
  expect_lt(as.numeric(s_large[3, "sd"]), as.numeric(s_small[3, "sd"]))
  expect_lt(as.numeric(s_large[4, "sd"]), as.numeric(s_small[4, "sd"]))
  # sd drift
  expect_lt(as.numeric(s_large[10, "sd"]), as.numeric(s_small[10, "sd"]))
  expect_lt(as.numeric(s_large[11, "sd"]), as.numeric(s_small[11, "sd"]))
  # larger data set takes longer to converge
  expect_lt(
    fit_small$fit$time()$total,
    fit_large$fit$time()$total
  )
})

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
