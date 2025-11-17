test_that("coev_fit() produces expected errors", {
  # simulate data
  withr::with_seed(1, {
    n <- 20
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = tree$tip.label,
      v = as.integer(rnbinom(n, mu = 4, size = 1)),
      w = rnorm(n),
      x = rbinom(n, size = 1, prob = 0.5),
      y = ordered(sample(1:4, size = n, replace = TRUE)),
      z = rpois(n, 3)
    )
  })
  # expect the following errors
  #' @srrstats {G5.2, G5.2b} Test all error messages
  expect_error(
    coev_fit(
      data = log, # should not accept functions
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree
    ),
    "Argument 'data' must be coercible to a data.frame.",
    fixed = TRUE
  )
  #' @srrstats {G5.8, G5.8a} Test for zero-length data
  expect_error(
    coev_fit(
      data = data.frame(), # empty
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree
    ),
    "Argument 'data' does not contain observations.",
    fixed = TRUE
  )
  expect_error(
    coev_fit(
      data = d,
      variables = "testing", # not a named list
      id = "id",
      tree = tree
    ),
    "Argument 'variables' is not a named list.",
    fixed = TRUE
  )
  expect_error(
    coev_fit(
      data = d,
      variables = list(
        # incorrect variable names (a/b instead of x/y)
        a = "bernoulli_logit",
        b = "ordered_logistic"
      ),
      id = "id",
      tree = tree
    ),
    "Some variable names are not valid column names in the data.",
    fixed = TRUE
  )
  expect_error(
    coev_fit(
      data = d,
      variables = list(
        x = "poisson", # incorrect response distribution
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree
    ),
    paste0(
      "Response distributions other than 'bernoulli_logit', ",
      "'ordered_logistic', 'poisson_softplus', ",
      "'negative_binomial_softplus', 'normal', and 'gamma_log' ",
      "are not yet supported."
    ),
    fixed = TRUE
  )
  expect_error(
    coev_fit(
      data = d,
      variables = list(
        # only x declared
        x = "bernoulli_logit"
      ),
      id = "id",
      tree = tree
    ),
    "Must be at least two coevolving variables.",
    fixed = TRUE
  )
  expect_error(
    coev_fit(
      data = dplyr::tibble(
        id = tree$tip.label,
        x = list("test"),
        y = list("test")
      ),
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree
    ),
    "Data must not contain list columns.",
    fixed = TRUE
  )
  expect_error(
    coev_fit(
      data = data.frame(
        id = tree$tip.label,
        x = Inf,
        y = -Inf
      ),
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree
    ),
    "Data must not contain undefined values (i.e., Inf or -Inf).",
    fixed = TRUE
  )
  #' @srrstats {G5.8, G5.8c} Test for columns with only NA values
  expect_error(
    coev_fit(
      data = data.frame(
        id = tree$tip.label,
        x = NA,
        y = NA
      ),
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree
    ),
    "Data must not contain columns with only NA values.",
    fixed = TRUE
  )
  #' @srrstats {G5.8, G5.8b, G5.8d} Tests for data of unsupported types
  expect_error(
    coev_fit(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "bernoulli_logit" # not binary integer
      ),
      id = "id",
      tree = tree
    ),
    paste0(
      "Variables following the 'bernoulli_logit' response distribution ",
      "must be integers with values of 0/1 in the data. Try using the ",
      "as.integer() function to convert variables to integers."
    ),
    fixed = TRUE
  )
  expect_error(
    coev_fit(
      data = d,
      variables = list(
        x = "ordered_logistic", # not ordered factor
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree
    ),
    paste0(
      "Variables following the 'ordered_logistic' response distribution ",
      "must be ordered factors in the data. Try using the as.ordered() ",
      "function to convert variables to ordered factors."
    ),
    fixed = TRUE
  )
  expect_error(
    coev_fit(
      data = d,
      variables = list(
        y = "poisson_softplus", # not integer >= 0
        z = "poisson_softplus"
      ),
      id = "id",
      tree = tree
    ),
    paste0(
      "Variables following the 'poisson_softplus' response distribution ",
      "must be integers greater than or equal to zero in the data. Try ",
      "using the as.integer() function to convert variables to integers."
    ),
    fixed = TRUE
  )
  expect_error(
    coev_fit(
      data = d,
      variables = list(
        y = "negative_binomial_softplus", # not integer >= 0
        z = "negative_binomial_softplus"
      ),
      id = "id",
      tree = tree
    ),
    paste0(
      "Variables following the 'negative_binomial_softplus' response ",
      "distribution must be integers greater than or equal to zero in ",
      "the data. Try using the as.integer() function to convert variables ",
      "to integers."
    ),
    fixed = TRUE
  )
  expect_error(
    coev_make_stancode(
      data = d,
      variables = list(
        v = "negative_binomial_softplus",
        z = "negative_binomial_softplus" # sd squared <= mean
      ),
      id = "id",
      tree = tree
    ),
    paste0(
      "No overdispersion or potentially underdispersion for ",
      "variable 'z' (sd^2 <= mean), do not use the ",
      "'negative_binomial_softplus' response distribution ",
      "for this variable."
    ),
    fixed = TRUE
  )
  expect_error(
    coev_fit(
      data = d,
      variables = list(
        w = "normal", # not numeric
        y = "normal"
      ),
      id = "id",
      tree = tree
    ),
    paste0(
      "Variables following the 'normal' response distribution ",
      "must be numeric in the data."
    ),
    fixed = TRUE
  )
  expect_error(
    coev_fit(
      data = d,
      variables = list(
        v = "gamma_log",
        w = "gamma_log" # not positive real
      ),
      id = "id",
      tree = tree
    ),
    paste0(
      "Variables following the 'gamma_log' response distribution must ",
      "be positive reals in the data."
    ),
    fixed = TRUE
  )
  expect_error(
    coev_fit(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = c(1, 2), # not character of length one
      tree = tree
    ),
    "Argument 'id' must be a character of length one.",
    fixed = TRUE
  )
  expect_error(
    coev_fit(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "testing", # incorrect id
      tree = tree
    ),
    "Argument 'id' is not a valid column name in the data.",
    fixed = TRUE
  )
  expect_error(
    coev_fit(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = "testing" # not of class phylo
    ),
    paste0(
      "Argument 'tree' must be a phylogenetic tree object of class phylo or ",
      "multiPhylo."
    ),
    fixed = TRUE
  )
  expect_error(
    coev_fit(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = ape::rtree(n, br = NULL) # no branch lengths
    ),
    "All trees in 'tree' argument must include branch lengths.",
    fixed = TRUE
  )
  expect_error(
    coev_fit(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = ape::unroot(tree) # unrooted
    ),
    "All trees in 'tree' argument must be rooted.",
    fixed = TRUE
  )
  expect_error(
    {
      tree2 <- tree
      tree2$edge.length[25] <- 0
      coev_fit(
        data = d,
        variables = list(
          x = "bernoulli_logit",
          y = "ordered_logistic"
        ),
        id = "id",
        tree = tree2
      )
    },
    "All trees in 'tree' argument must have positive non-zero branch lengths.",
    fixed = TRUE
  )
  expect_error(
    {
      d2 <- d
      d2$id <- rep("test", n) # not correct tip labels
      coev_fit(
        data = d2,
        variables = list(
          x = "bernoulli_logit",
          y = "ordered_logistic"
        ),
        id = "id",
        tree = tree
      )
    },
    "The id variable in the data does not match tree tip labels exactly.",
    fixed = TRUE
  )
  expect_error(
    {
      tree2 <- c(tree, ape::di2multi(tree, tol = 0.01)) # collapse internal node
      coev_fit(
        data = d,
        variables = list(
          x = "bernoulli_logit",
          y = "ordered_logistic"
        ),
        id = "id",
        tree = tree2
      )
    },
    paste0(
      "All trees in 'tree' argument must have the same number of ",
      "internal nodes and branches."
    ),
    fixed = TRUE
  )
  expect_error(
    {
      d2 <- d
      d2$id[1] <- NA
      tree2 <- tree
      tree2$tip.label[1] <- NA
      coev_fit(
        data = d2,
        variables = list(
          x = "bernoulli_logit",
          y = "ordered_logistic"
        ),
        id = "id",
        tree = tree2
      )
    },
    "The id variable in the data must not contain NAs.",
    fixed = TRUE
  )
  expect_error(
    coev_fit(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree,
      effects_mat = "testing" # not of class matrix
    ),
    "Argument 'effects_mat' must be a matrix.",
    fixed = TRUE
  )
  expect_error(
    coev_fit(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree,
      effects_mat = matrix(1) # not boolean matrix
    ),
    "Argument 'effects_mat' must be a boolean matrix.",
    fixed = TRUE
  )
  expect_error(
    coev_fit(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree,
      effects_mat = matrix(TRUE) # no row/col names
    ),
    "Argument 'effects_mat' does not have valid row or column names.",
    fixed = TRUE
  )
  expect_error(
    coev_fit(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree,
      effects_mat = matrix(TRUE, dimnames = list("fail", "fail")) # invalid name
    ),
    paste0(
      "Row and column names for argument 'effects_mat' do not match ",
      "variable names exactly."
    ),
    fixed = TRUE
  )
  expect_error(
    coev_fit(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree,
      # autoregressive effect = FALSE
      effects_mat = matrix(c(TRUE, TRUE, TRUE, FALSE),
                           ncol = 2, nrow = 2, byrow = TRUE,
                           dimnames = list(c("x", "y"), c("x", "y")))
    ),
    "Argument 'effects_mat' must specify TRUE for all autoregressive effects.",
    fixed = TRUE
  )
  expect_error(
    coev_fit(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree,
      complete_cases = "testing"
    ),
    "Argument 'complete_cases' must be a logical of length one.",
    fixed = TRUE
  )
  expect_error(
    coev_fit(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree,
      dist_mat = "testing" # not of class matrix
    ),
    "Argument 'dist_mat' must be a matrix.",
    fixed = TRUE
  )
  expect_error(
    coev_fit(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree,
      dist_mat = matrix(letters) # matrix not numeric
    ),
    "Argument 'dist_mat' must be a numeric matrix.",
    fixed = TRUE
  )
  expect_error(
    coev_fit(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree,
      dist_mat = matrix(1:100, nrow = 10) # matrix not symmetric
    ),
    "Argument 'dist_mat' must be a symmetric matrix.",
    fixed = TRUE
  )
  expect_error(
    coev_fit(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree,
      # matrix symmetric but diagonal not zero
      dist_mat = matrix(rep(1, 100), nrow = 10)
    ),
    "Argument 'dist_mat' must have zeroes on the diagonal of the matrix.",
    fixed = TRUE
  )
  expect_error(
    coev_fit(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree,
      dist_mat = matrix(0) # no row/col names
    ),
    "Argument 'dist_mat' does not have valid row or column names.",
    fixed = TRUE
  )
  expect_error(
    coev_fit(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree,
      dist_mat = matrix(rep(0, 100), nrow = 10,
                        # matrix row/col names do not match tips
                        dimnames = list(letters[1:10], letters[1:10]))
    ),
    paste0(
      "Row and column names for argument 'dist_mat' do not ",
      "match tree tip labels exactly."
    ),
    fixed = TRUE
  )
  expect_error(
    coev_fit(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree,
      dist_cov = FALSE
    ),
    "Argument 'dist_cov' is not a character string.",
    fixed = TRUE
  )
  expect_error(
    coev_fit(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree,
      dist_cov = c("fail", "fail")
    ),
    "Argument 'dist_cov' is not of length 1.",
    fixed = TRUE
  )
  expect_error(
    coev_fit(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree,
      dist_cov = "fail"
    ),
    paste0(
      "Argument 'dist_cov' currently only supports 'exp_quad', ",
      "'exponential', and 'matern32'."
    ),
    fixed = TRUE
  )
  expect_error(
    coev_fit(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree,
      prior = "testing" # not a list
    ),
    "Argument 'prior' is not a list.",
    fixed = TRUE
  )
  expect_error(
    coev_fit(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree,
      prior = list("testing") # not a named list
    ),
    "Argument 'prior' is not a named list.",
    fixed = TRUE
  )
  expect_error(
    coev_fit(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree,
      prior = list(testing = "testing") # incorrect name
    ),
    paste0(
      "Argument 'prior' list contains names that are not allowed. Please ",
      "use only the following names: 'b', 'eta_anc', 'A_offdiag', 'A_diag', ",
      "'L_R', 'Q_sigma', 'c', 'phi', 'shape', 'sigma_dist', 'rho_dist', ",
      "'sigma_residual', and 'L_residual'"
    ),
    fixed = TRUE
  )
  expect_error(
    coev_fit(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree,
      # duplicate names
      prior = list(A_diag = "normal(0,2)", A_diag = "normal(0,2)")
    ),
    "Argument 'prior' contains duplicate names.",
    fixed = TRUE
  )
  expect_error(
    coev_fit(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree,
      scale = "testing"
    ),
    "Argument 'scale' must be a logical of length one.",
    fixed = TRUE
  )
  expect_error(
    coev_fit(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree,
      estimate_correlated_drift = "testing"
    ),
    "Argument 'estimate_correlated_drift' must be a logical of length one.",
    fixed = TRUE
  )
  expect_error(
    coev_fit(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree,
      estimate_residual = "testing"
    ),
    "Argument 'estimate_residual' must be a logical of length one.",
    fixed = TRUE
  )
  expect_error(
    coev_fit(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree,
      log_lik = "testing"
    ),
    "Argument 'log_lik' must be a logical of length one.",
    fixed = TRUE
  )
  expect_error(
    coev_fit(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree,
      prior_only = "testing"
    ),
    "Argument 'prior_only' must be a logical of length one.",
    fixed = TRUE
  )
  expect_error(
    coev_fit(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree,
      backend = 0
    ),
    "Argument 'backend' must be a character string of length one.",
    fixed = TRUE
  )
  expect_error(
    coev_fit(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree,
      backend = "testing"
    ),
    "Argument 'backend' must be either 'cmdstanr' or 'nutpie'.",
    fixed = TRUE
  )
})

test_that("coev_fit() fits simple model without error", {
  # simulate data
  withr::with_seed(1, {
    n <- 3
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = tree$tip.label,
      x = rnorm(n),
      y = rnorm(n)
    )
  })
  # fit model
  fit <-
    coev_fit(
      data = d,
      variables = list(
        x = "normal",
        y = "normal"
      ),
      id = "id",
      tree = tree,
      chains = 1,
      seed = 1,
      refresh = 0
    )
  # expect no errors for model fitting or summaries
  sw <- suppressWarnings
  expect_no_error(sw(fit))
  expect_no_error(sw(summary(fit)))
  expect_no_error(sw(print(fit)))
  expect_no_error(sw(print(summary(fit))))
  # expect no error for stancode and standata methods
  expect_no_error(sw(stancode(fit)))
  expect_no_error(sw(standata(fit)))
  expect_output(sw(stancode(fit)))
  expect_true(sw(methods::is(standata(fit), "list")))
})

test_that("coev_fit() fits test fixtures", {
  # load models
  m1 <- readRDS(test_path("fixtures", "coevfit_example_01.rds"))
  m2 <- readRDS(test_path("fixtures", "coevfit_example_02.rds"))
  m3 <- readRDS(test_path("fixtures", "coevfit_example_03.rds"))
  m4 <- readRDS(test_path("fixtures", "coevfit_example_04.rds"))
  m1 <- reload_fit(m1, filename = "coevfit_example_01-1.csv")
  m2 <- reload_fit(m2, filename = "coevfit_example_02-1.csv")
  m3 <- reload_fit(m3, filename = "coevfit_example_03-1.csv")
  m4 <- reload_fit(m4, filename = "coevfit_example_04-1.csv")
  # suppress warnings
  sw <- suppressWarnings
  # expect no errors for model fitting or summaries
  expect_no_error(sw(m1))
  expect_no_error(sw(m2))
  expect_no_error(sw(m3))
  expect_no_error(sw(m4))
  expect_no_error(sw(summary(m1)))
  expect_no_error(sw(summary(m2)))
  expect_no_error(sw(summary(m3)))
  expect_no_error(sw(summary(m4)))
  expect_output(sw(print(m1)))
  expect_output(sw(print(m2)))
  expect_output(sw(print(m3)))
  expect_output(sw(print(m4)))
  expect_output(sw(print(summary(m1))))
  expect_output(sw(print(summary(m2))))
  expect_output(sw(print(summary(m3))))
  expect_output(sw(print(summary(m4))))
  # expect error if prob for summary is outside of range 0 - 1
  expect_error(sw(summary(m1, prob = -0.01)))
  expect_error(sw(summary(m1, prob =  1.01)))
  # expect no errors for extract_samples method
  expect_no_error(sw(extract_samples(m1)))
  expect_no_error(sw(extract_samples(m2)))
  expect_no_error(sw(extract_samples(m3)))
  expect_no_error(sw(extract_samples(m4)))
  expect_true(sw(methods::is(extract_samples(m1), "list")))
  expect_true(sw(methods::is(extract_samples(m2), "list")))
  expect_true(sw(methods::is(extract_samples(m3), "list")))
  expect_true(sw(methods::is(extract_samples(m4), "list")))
  # expect no error for stancode and standata methods
  expect_no_error(sw(stancode(m1)))
  expect_no_error(sw(stancode(m2)))
  expect_no_error(sw(stancode(m3)))
  expect_no_error(sw(stancode(m4)))
  expect_no_error(sw(standata(m1)))
  expect_no_error(sw(standata(m2)))
  expect_no_error(sw(standata(m3)))
  expect_no_error(sw(standata(m4)))
  expect_output(sw(stancode(m1)))
  expect_output(sw(stancode(m2)))
  expect_output(sw(stancode(m3)))
  expect_output(sw(stancode(m4)))
  expect_true(sw(methods::is(standata(m1), "list")))
  expect_true(sw(methods::is(standata(m2), "list")))
  expect_true(sw(methods::is(standata(m3), "list")))
  expect_true(sw(methods::is(standata(m4), "list")))
  # expect following errors for plot method
  expect_error(
    plot(m4, parameters = "test"),
    "Argument 'parameters' contains invalid parameter names.",
    fixed = TRUE
  )
  expect_error(
    plot(m4, npars = "test"),
    "Argument 'npars' is not a positive number.",
    fixed = TRUE
  )
  expect_error(
    plot(m4, plot = "test"),
    "Argument 'plot' is not logical.",
    fixed = TRUE
  )
  # plot method works as expected
  plot(m4)
  plot(m4, parameters = "A[1,1]")
  plot(m4, combo = c("hist", "trace"))
  plot(m4, npars = 3)
  plot(m4, plot = FALSE)
})

test_that("effects_mat argument to coev_fit() works as expected", {
  # load model
  m <- readRDS(test_path("fixtures", "coevfit_example_05.rds"))
  m <- reload_fit(m, filename = "coevfit_example_05-1.csv")
  # suppress warnings
  sw <- suppressWarnings
  # expect no errors for model fitting or summaries
  expect_no_error(sw(m))
  expect_no_error(sw(summary(m)))
  expect_output(sw(print(m)))
  expect_output(sw(print(summary(m))))
  # expect no errors for extract_samples method
  expect_no_error(sw(extract_samples(m)))
  expect_true(sw(methods::is(extract_samples(m), "list")))
  # expect no errors for stancode or standata methods
  expect_no_error(sw(stancode(m)))
  expect_no_error(sw(standata(m)))
  expect_output(sw(stancode(m)))
  expect_true(sw(methods::is(standata(m), "list")))
  # expect effects_mat correct in model output
  effects_mat <- matrix(
    c(TRUE, TRUE,
      FALSE, TRUE),
    byrow = TRUE,
    nrow = 2,
    ncol = 2,
    dimnames = list(c("w", "x"), c("w", "x"))
  )
  expect_true(
    identical(m$effects_mat, +effects_mat)
  )
  # correct parameter estimated to be zero (A[2,1])
  expect_true(all(as.vector(m$fit$draws()[, , "A[2,1]"]) == 0))
  # other parameters estimated as normal
  expect_true(!all(as.vector(m$fit$draws()[, , "A[1,1]"]) == 0))
  expect_true(!all(as.vector(m$fit$draws()[, , "A[1,2]"]) == 0))
  expect_true(!all(as.vector(m$fit$draws()[, , "A[2,2]"]) == 0))
})

test_that("coev_fit() works with missing data", {
  # load model
  m <- readRDS(test_path("fixtures", "coevfit_example_06.rds"))
  m <- reload_fit(m, filename = "coevfit_example_06-1.csv")
  # suppress warnings
  sw <- suppressWarnings
  # fitted without error
  expect_no_error(sw(m))
  expect_no_error(sw(summary(m)))
  expect_output(sw(print(m)))
  expect_output(sw(print(summary(m))))
  # expect no errors for extract_samples method
  expect_no_error(sw(extract_samples(m)))
  expect_true(sw(methods::is(extract_samples(m), "list")))
  # expect no errors for stancode or standata methods
  expect_no_error(sw(stancode(m)))
  expect_no_error(sw(standata(m)))
  expect_output(sw(stancode(m)))
  expect_true(sw(methods::is(standata(m), "list")))
})

test_that("coev_fit() works with repeated observations", {
  # load model
  m <- readRDS(test_path("fixtures", "coevfit_example_07.rds"))
  m <- reload_fit(m, filename = "coevfit_example_07-1.csv")
  # suppress warnings
  sw <- suppressWarnings
  # fitted without error
  expect_no_error(sw(m))
  expect_no_error(sw(summary(m)))
  expect_output(sw(print(m)))
  expect_output(sw(print(summary(m))))
  # expect no errors for stancode or standata methods
  expect_no_error(sw(stancode(m)))
  expect_no_error(sw(standata(m)))
  expect_output(sw(stancode(m)))
  expect_true(sw(methods::is(standata(m), "list")))
})

test_that("coev_fit() works with multiPhylo object", {
  # load model
  m <- readRDS(test_path("fixtures", "coevfit_example_08.rds"))
  m <- reload_fit(m, filename = "coevfit_example_08-1.csv")
  # suppress warnings
  sw <- suppressWarnings
  # fitted without error
  expect_no_error(sw(m))
  expect_no_error(sw(summary(m)))
  expect_output(sw(print(m)))
  expect_output(sw(print(summary(m))))
  # expect no errors for stancode or standata methods
  expect_no_error(sw(stancode(m)))
  expect_no_error(sw(standata(m)))
  expect_output(sw(stancode(m)))
  expect_true(sw(methods::is(standata(m), "list")))
})

test_that("coev_fit() works when Q off diagonals == 0", {
  # load model
  m <- readRDS(test_path("fixtures", "coevfit_example_09.rds"))
  m <- reload_fit(m, filename = "coevfit_example_09-1.csv")
  # suppress warnings
  sw <- suppressWarnings
  # fitted without error
  expect_no_error(sw(m))
  expect_no_error(sw(summary(m)))
  expect_output(sw(print(m)))
  expect_output(sw(print(summary(m))))
  # expect no errors for stancode or standata methods
  expect_no_error(sw(stancode(m)))
  expect_no_error(sw(standata(m)))
  expect_output(sw(stancode(m)))
  expect_true(sw(methods::is(standata(m), "list")))
})

test_that("coev_fit() works with measurement error", {
  # load model
  m <- readRDS(test_path("fixtures", "coevfit_example_10.rds"))
  m <- reload_fit(m, filename = "coevfit_example_10-1.csv")
  # suppress warnings
  sw <- suppressWarnings
  # fitted without error
  expect_no_error(sw(m))
  expect_no_error(sw(summary(m)))
  expect_output(sw(print(m)))
  expect_output(sw(print(summary(m))))
  # expect no errors for stancode or standata methods
  expect_no_error(sw(stancode(m)))
  expect_no_error(sw(standata(m)))
  expect_output(sw(stancode(m)))
  expect_true(sw(methods::is(standata(m), "list")))
})
