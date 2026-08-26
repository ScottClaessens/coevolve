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

test_that("summary.jax_fit forwards quantile formulas as columns", {
  draws <- posterior::as_draws_array(
    array(
      stats::rnorm(4 * 100 * 2),
      dim = c(100, 4, 2),
      dimnames = list(
        iteration = NULL,
        chain = NULL,
        variable = c("a", "b")
      )
    )
  )
  fit <- coevolve:::create_jax_wrapper(
    trace_result = NULL,
    draws_array = draws,
    stan_variables = c("a", "b"),
    iter_sampling = 100L,
    iter_warmup = 100L,
    chains = 4L
  )
  s <- fit$summary(
    NULL,
    Estimate = "mean",
    `Est.Error` = "sd",
    ~quantile(.x, probs = c(0.025, 0.975)),
    Rhat = "rhat"
  )
  expect_true(all(c("Estimate", "Est.Error", "2.5%", "97.5%", "Rhat")
                  %in% names(s)))
})

test_that("coev_fit with nuts_sampler='nutpie' returns coevfit", {
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

test_that("Plotting functions work with JAX fit", {
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
    prior = list(A_offdiag = "normal(0, 2)"),
    nuts_sampler = "nutpie",
    chains = 1,
    iter_warmup = 50,
    iter_sampling = 50,
    seed = 1
  )
  expect_no_error(coev_plot_trait_values(fit))
  expect_no_error(coev_plot_delta_theta(fit))
  expect_no_error(coev_plot_pred_series(fit))
  expect_no_error(
    coev_plot_flowfield(
      fit, var1 = "political_authority", var2 = "religious_authority"
    )
  )
  expect_no_error(
    coev_plot_selection_gradient(
      fit, var1 = "political_authority", var2 = "religious_authority"
    )
  )
})

test_that("JAX fit generates yrep posterior predictions", {
  skip_if_not(
    coevolve:::check_jax_available(),
    message = "JAX not available - skipping JAX tests"
  )
  fit_args <- list(
    data = authority$data,
    variables = list(
      political_authority = "ordered_logistic",
      religious_authority = "ordered_logistic"
    ),
    id = "language",
    tree = authority$phylogeny,
    prior = list(A_offdiag = "normal(0, 2)"),
    nuts_sampler = "nutpie",
    chains = 2,
    iter_warmup = 100,
    iter_sampling = 100,
    seed = 1
  )
  fit <- do.call(coev_fit, fit_args)
  sd <- fit$stan_data
  # yrep is returned with the same dimensions as the Stan backend
  yrep <- posterior::draws_of(
    posterior::as_draws_rvars(fit$fit$draws(variables = "yrep"))$yrep
  )
  expect_equal(
    dim(yrep), c(fit$nsamples, sd$N_tree, sd$N_obs, sd$J)
  )
  # ordinal predictions fall within the observed category range
  expect_true(all(yrep == round(yrep)))
  expect_true(all(yrep >= 1 & yrep <= max(sd$y)))
  # predictions vary across draws and across chains
  expect_gt(stats::sd(yrep[, 1, 1, 1]), 0)
  draws_arr <- fit$fit$draws(variables = "yrep")
  expect_false(
    identical(as.vector(draws_arr[1, 1, ]), as.vector(draws_arr[1, 2, ]))
  )
  # predictive checks now work for nutpie fits (#118)
  sw <- suppressWarnings
  expect_no_error(pp <- sw(coev_plot_predictive_check(fit)))
  expect_named(pp, names(fit$variables))
  expect_no_error(sw(coev_plot_predictive_check(fit, ndraws = 1)))
  # same seed reproduces the same predictions
  expect_equal(
    posterior::draws_of(
      posterior::as_draws_rvars(
        do.call(coev_fit, fit_args)$fit$draws(variables = "yrep")
      )$yrep
    ),
    yrep
  )
})

test_that("JAX yrep works with repeated observations and normal variables", {
  skip_if_not(
    coevolve:::check_jax_available(),
    message = "JAX not available - skipping JAX tests"
  )
  fit <- coev_fit(
    data = repeated$data,
    variables = list(
      x = "normal",
      y = "normal"
    ),
    id = "species",
    tree = repeated$phylogeny,
    nuts_sampler = "nutpie",
    chains = 1,
    iter_warmup = 100,
    iter_sampling = 100,
    seed = 1
  )
  sd <- fit$stan_data
  yrep <- posterior::draws_of(
    posterior::as_draws_rvars(fit$fit$draws(variables = "yrep"))$yrep
  )
  expect_equal(
    dim(yrep), c(fit$nsamples, sd$N_tree, sd$N_obs, sd$J)
  )
  expect_true(all(is.finite(yrep)))
  expect_no_error(suppressWarnings(coev_plot_predictive_check(fit)))
})

test_that("JAX fit returns log_lik when log_lik = TRUE", {
  skip_if_not(
    coevolve:::check_jax_available(),
    message = "JAX not available - skipping JAX tests"
  )
  fit_args <- list(
    data = authority$data,
    variables = list(
      political_authority = "ordered_logistic",
      religious_authority = "ordered_logistic"
    ),
    id = "language",
    tree = authority$phylogeny,
    prior = list(A_offdiag = "normal(0, 2)"),
    nuts_sampler = "nutpie",
    chains = 1,
    iter_warmup = 100,
    iter_sampling = 100,
    seed = 1
  )
  # log_lik is off by default, matching the Stan backend
  fit <- do.call(coev_fit, fit_args)
  expect_false(
    any(grepl("^log_lik", posterior::variables(fit$fit$draws())))
  )
  fit <- do.call(coev_fit, c(fit_args, list(log_lik = TRUE)))
  sd <- fit$stan_data
  log_lik <- posterior::draws_of(
    posterior::as_draws_rvars(fit$fit$draws(variables = "log_lik"))$log_lik
  )
  expect_equal(dim(log_lik), c(fit$nsamples, sd$N_obs * sd$J))
  expect_true(all(is.finite(log_lik)))
  expect_true(all(log_lik <= 0))
  expect_no_error(suppressWarnings(summary(fit)))
  skip_if_not_installed("loo")
  expect_no_error(suppressWarnings(loo::loo(log_lik)))
})

test_that("JAX log_lik handles normal variables and missing data", {
  skip_if_not(
    coevolve:::check_jax_available(),
    message = "JAX not available - skipping JAX tests"
  )
  withr::with_seed(1, {
    n <- 15
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = tree$tip.label,
      x = rnorm(n),
      y = rnorm(n)
    )
  })
  d$x[1:3] <- NA
  fit <- coev_fit(
    data = d,
    variables = list(x = "normal", y = "normal"),
    id = "id",
    tree = tree,
    nuts_sampler = "nutpie",
    log_lik = TRUE,
    chains = 1,
    iter_warmup = 100,
    iter_sampling = 100,
    seed = 1
  )
  sd <- fit$stan_data
  log_lik <- posterior::draws_of(
    posterior::as_draws_rvars(fit$fit$draws(variables = "log_lik"))$log_lik
  )
  expect_equal(dim(log_lik), c(fit$nsamples, sd$N_obs * sd$J))
  expect_true(all(is.finite(log_lik)))
  # missing cells contribute nothing (log_sum_exp of zeros over one tree)
  miss_ids <- which(as.vector(sd$miss) == 1)
  expect_true(all(log_lik[, miss_ids] == 0))
  expect_true(all(log_lik[, -miss_ids] != 0))
})

test_that("JAX fit works with prior_only for normal variables", {
  skip_if_not(
    coevolve:::check_jax_available(),
    message = "JAX not available - skipping JAX tests"
  )
  withr::with_seed(3, {
    n <- 8
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = tree$tip.label,
      x = rnorm(n),
      y = rnorm(n)
    )
  })
  # terminal_drift has no prior of its own outside the likelihood block, so
  # prior-only sampling used to be improper and failed to initialise (#118)
  fit <- coev_fit(
    data = d,
    variables = list(x = "normal", y = "normal"),
    id = "id",
    tree = tree,
    nuts_sampler = "nutpie",
    prior_only = TRUE,
    chains = 2,
    iter_warmup = 300,
    iter_sampling = 300,
    seed = 1
  )
  expect_s3_class(fit, "coevfit")
  rhat <- posterior::summarise_draws(
    posterior::subset_draws(fit$fit$draws(), variable = "eta_anc"),
    "rhat"
  )$rhat
  expect_true(all(rhat < 1.05))
  expect_no_error(suppressWarnings(coev_plot_predictive_check(fit)))
})

test_that("coev_ancestral_states() works with JAX fit", {
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
    prior = list(A_offdiag = "normal(0, 2)"),
    nuts_sampler = "nutpie",
    chains = 1,
    iter_warmup = 50,
    iter_sampling = 50,
    seed = 1
  )
  # latent scale summary
  asr <- coev_ancestral_states(fit)
  expect_true(is.data.frame(asr))
  expect_true(all(
    c("node", "variable", "estimate", "lower", "upper") %in% names(asr)
  ))
  expect_true(all(is.finite(asr$estimate)))
  # response scale summary (exercises cutpoint extraction)
  asr_resp <- coev_ancestral_states(fit, scale = "response")
  expect_true("category" %in% names(asr_resp))
  expect_true(all(asr_resp$estimate >= 0 & asr_resp$estimate <= 1))
  # raw draws
  asr_draws <- coev_ancestral_states(fit, summary = FALSE)
  expect_true(is.array(asr_draws$draws))
  expect_equal(length(dim(asr_draws$draws)), 3)
})
