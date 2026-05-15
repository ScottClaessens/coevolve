#' @srrstats {G5.10} Flag extended tests
run_extended_tests <- identical(Sys.getenv("COEVOLVE_EXTENDED_TESTS"), "true")

#' @srrstats {G5.6, G5.6a, G5.6b, G5.9, G5.9a, G5.9b, BS7.2} Extended parameter
#'   recovery tests with multiple fixed seeds for data simulation and cmdstanr
#' @srrstats {BS7.4, BS7.4a} Predicted values are on same scale as input data
for (seed in 1:2) {
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
