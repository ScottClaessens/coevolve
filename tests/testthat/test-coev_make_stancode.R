test_that("coev_make_stancode() produces expected errors", {
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
  expect_error(
    coev_make_stancode(
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
  expect_error(
    coev_make_stancode(
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
    coev_make_stancode(
      data = d,
      variables = "testing", # not a named list
      id = "id",
      tree = tree
    ),
    "Argument 'variables' is not a named list.",
    fixed = TRUE
  )
  expect_error(
    coev_make_stancode(
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
    coev_make_stancode(
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
    coev_make_stancode(
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
    coev_make_stancode(
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
    coev_make_stancode(
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
    coev_make_stancode(
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
    coev_make_stancode(
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
        z = "negative_binomial_softplus" # sd^2 <= mean
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
    coev_make_stancode(
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
    coev_make_stancode(
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
    coev_make_stancode(
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
    coev_make_stancode(
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
    coev_make_stancode(
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
    coev_make_stancode(
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
    coev_make_stancode(
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
      coev_make_stancode(
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
      coev_make_stancode(
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
      coev_make_stancode(
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
      d2 <- d; d2$id[1] <- NA
      tree2 <- tree; tree2$tip.label[1] <- NA
      coev_make_stancode(
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
    coev_make_stancode(
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
    coev_make_stancode(
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
    coev_make_stancode(
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
    coev_make_stancode(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree,
      effects_mat = matrix(TRUE, dimnames = list("fail","fail")) # invalid names
    ),
    paste0(
      "Row and column names for argument 'effects_mat' do not match ",
      "variable names exactly."
    ),
    fixed = TRUE
  )
  expect_error(
    coev_make_stancode(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree,
      # autoregressive effect = FALSE
      effects_mat = matrix(c(T,T,T,F), ncol = 2, nrow = 2, byrow = TRUE,
                           dimnames = list(c("x","y"),c("x","y")))
    ),
    "Argument 'effects_mat' must specify TRUE for all autoregressive effects.",
    fixed = TRUE
  )
  expect_error(
    coev_make_stancode(
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
    coev_make_stancode(
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
    coev_make_stancode(
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
    coev_make_stancode(
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
    coev_make_stancode(
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
    coev_make_stancode(
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
    coev_make_stancode(
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
    coev_make_stancode(
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
    coev_make_stancode(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree,
      dist_cov = c("fail","fail")
    ),
    "Argument 'dist_cov' is not of length 1.",
    fixed = TRUE
  )
  expect_error(
    coev_make_stancode(
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
    coev_make_stancode(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        w = "normal"
      ),
      id = "id",
      tree = tree,
      measurement_error = "fail"
    ),
    "Argument 'measurement_error' is not a named list.",
    fixed = TRUE
  )
  expect_error(
    coev_make_stancode(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        w = "normal"
      ),
      id = "id",
      tree = tree,
      measurement_error = list(x = "x_se")
    ),
    paste0(
      "Argument 'measurement_error' contains variables that were not ",
      "declared as normally-distributed variables in the model."
    ),
    fixed = TRUE
  )
  expect_error(
    coev_make_stancode(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        w = "normal"
      ),
      id = "id",
      tree = tree,
      measurement_error = list(w = "w_se")
    ),
    paste0(
      "Argument 'measurement_error' refers to measurement error columns ",
      "that are not valid column names in the data."
    ),
    fixed = TRUE
  )
  d$w_se <- rep(-1, 20)
  expect_error(
    coev_make_stancode(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        w = "normal"
      ),
      id = "id",
      tree = tree,
      measurement_error = list(w = "w_se")
    ),
    paste0(
      "Standard errors in measurement error columns must be zero or ",
      "positive reals."
    ),
    fixed = TRUE
  )
  d$w_se <- rexp(20)
  d$w_se[1] <- NA
  expect_error(
    coev_make_stancode(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        w = "normal"
      ),
      id = "id",
      tree = tree,
      measurement_error = list(w = "w_se")
    ),
    paste0(
      "Standard errors in measurement error columns must not be NA ",
      "in rows where there is observed data for the focal variable."
    ),
    fixed = TRUE
  )
  expect_error(
    coev_make_stancode(
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
    coev_make_stancode(
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
    coev_make_stancode(
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
      "'sigma_group', and 'L_group'"
    ),
    fixed = TRUE
  )
  expect_error(
    coev_make_stancode(
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
    coev_make_stancode(
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
    coev_make_stancode(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree,
      estimate_Q_offdiag = "testing"
    ),
    "Argument 'estimate_Q_offdiag' must be a logical of length one.",
    fixed = TRUE
  )
  expect_error(
    coev_make_stancode(
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
})

test_that("coev_make_stancode() returns a character of length one", {
  # simulate data
  withr::with_seed(1, {
    n <- 20
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = tree$tip.label,
      x = rbinom(n, size = 1, prob = 0.5),
      y = ordered(sample(1:4, size = n, replace = TRUE))
    )
  })
  # make stan code
  sc <-
    coev_make_stancode(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree
    )
  # expect string of length one
  expect_no_error(sc)
  expect_type(sc, "character")
  expect_length(sc, 1)
})

test_that("coev_make_stancode() creates Stan code with correct syntax", {
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
  # make stan code without distance matrix
  sc1 <-
    coev_make_stancode(
      data = d,
      variables = list(
          v = "negative_binomial_softplus",
          w = "normal",
          x = "bernoulli_logit",
          y = "ordered_logistic",
          z = "poisson_softplus"
      ),
      id = "id",
      tree = tree
    )
  # make stan code with distance matrix
  dist_mat <- as.matrix(dist(rnorm(n)))
  rownames(dist_mat) <- colnames(dist_mat) <- d$id
  sc2 <-
    coev_make_stancode(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic",
        z = "poisson_softplus"
      ),
      id = "id",
      tree = tree,
      dist_mat = dist_mat
    )
  # check stan code is syntactically correct
  expect_no_error(sc1)
  expect_no_error(sc2)
  expect_true(
    cmdstanr::cmdstan_model(
      stan_file = cmdstanr::write_stan_file(sc1),
      compile = FALSE
    )$check_syntax(quiet = TRUE)
  )
  expect_true(
    cmdstanr::cmdstan_model(
      stan_file = cmdstanr::write_stan_file(sc2),
      compile = FALSE
    )$check_syntax(quiet = TRUE)
  )
  # expect message when distance matrix included
  expect_message(
    coev_make_stancode(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic",
        z = "poisson_softplus"
      ),
      id = "id",
      tree = tree,
      dist_mat = dist_mat
    ),
    paste0(
      "Note: Distance matrix detected. Gaussian processes over spatial ",
      "distances have been included for each variable in the model."
    )
  )
})

test_that("Setting manual priors in coev_make_stancode() works as expected", {
  # simulate data
  withr::with_seed(1, {
    n <- 20
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = tree$tip.label,
      x = rbinom(n, size = 1, prob = 0.5),
      y = rbinom(n, size = 1, prob = 0.5)
    )
  })
  # invalid priors should throw an error
  suppressMessages({
    expect_error(
      coev_make_stancode(
        data = d,
        variables = list(
          x = "bernoulli_logit",
          y = "bernoulli_logit"
        ),
        id = "id",
        tree = tree,
        prior = list(A_diag = "testing")
      )
    )
  })
  # valid priors should throw no errors
  expect_no_error(
    coev_make_stancode(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "bernoulli_logit"
      ),
      id = "id",
      tree = tree,
      prior = list(
        b           = "normal(0, 2)",
        eta_anc     = "normal(0, 2)",
        A_offdiag   = "normal(0, 2)",
        A_diag      = "normal(0, 2)",
        L_R         = "lkj_corr_cholesky(3)",
        Q_sigma     = "normal(0, 2)",
        c           = "normal(0, 3)",
        phi         = "normal(1, 1)",
        shape       = "gamma(0.02, 0.02)",
        sigma_dist  = "exponential(2)",
        rho_dist    = "exponential(6)",
        sigma_group = "exponential(2)",
        L_group     = "lkj_corr_cholesky(3)"
      )
    )
  )
})

test_that("coev_make_stancode() produces neg binomial priors as expected", {
  # simulate data
  withr::with_seed(1, {
    n <- 20
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = tree$tip.label,
      x = as.integer(rnbinom(n, mu = 2, size = 1)),
      y = as.integer(rnbinom(n, mu = 2, size = 1))
    )
  })
  # default prior for phi should work
  expect_no_error(
    coev_make_stancode(
      data = d,
      variables = list(
        x = "negative_binomial_softplus",
        y = "negative_binomial_softplus"
      ),
      id = "id",
      tree = tree
    )
  )
  # setting a prior for phi should also work
  expect_no_error(
    coev_make_stancode(
      data = d,
      variables = list(
        x = "negative_binomial_softplus",
        y = "negative_binomial_softplus"
      ),
      id = "id",
      tree = tree,
      prior = list(phi = "std_normal()")
    )
  )
})

test_that("coev_make_stancode() works with missing data", {
  # simulate data
  withr::with_seed(1, {
    n <- 20
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = tree$tip.label,
      x = rbinom(n, size = 1, prob = 0.5),
      y = rbinom(n, size = 1, prob = 0.5)
    )
  })
  # some data are missing
  d$x[c(1,2)] <- NA
  d$y[c(1,3)] <- NA
  # get stan code
  sc <- coev_make_stancode(
    data = d,
    variables = list(
      x = "bernoulli_logit",
      y = "bernoulli_logit"
    ),
    id = "id",
    tree = tree
  )
  # runs without error
  expect_no_error(sc)
  # syntactically correct
  expect_true(
    cmdstanr::cmdstan_model(
      stan_file = cmdstanr::write_stan_file(sc),
      compile = FALSE
    )$check_syntax(quiet = TRUE)
  )
})

test_that("coev_make_stancode() works with repeated observations", {
  # simulate data
  withr::with_seed(1, {
    n <- 10
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = rep(tree$tip.label, each = 10),
      x = rbinom(n*10, size = 1, prob = 0.5),
      y = rbinom(n*10, size = 1, prob = 0.5)
    )
  })
  # get stan code
  sc <- coev_make_stancode(
    data = d,
    variables = list(
      x = "bernoulli_logit",
      y = "bernoulli_logit"
    ),
    id = "id",
    tree = tree
  )
  # runs without error
  expect_no_error(sc)
  # syntactically correct
  expect_true(
    cmdstanr::cmdstan_model(
      stan_file = cmdstanr::write_stan_file(sc),
      compile = FALSE
    )$check_syntax(quiet = TRUE)
  )
  # expect message with repeated observations
  expect_message(
    coev_make_stancode(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "bernoulli_logit"
      ),
      id = "id",
      tree = tree
    ),
    paste0(
      "Note: Repeated observations detected. Group-level varying effects ",
      "have been included for each variable in the model."
    )
  )
})

test_that("coev_make_stancode() works with tibbles", {
  # simulate data
  withr::with_seed(1, {
    n <- 20
    tree <- ape::rcoal(n)
    d <- tibble::tibble(
      id = tree$tip.label,
      x = rbinom(n, size = 1, prob = 0.5),
      y = ordered(sample(1:4, size = n, replace = TRUE))
    )
  })
  # make stan code
  sc <-
    coev_make_stancode(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree
    )
  # expect string of length one
  expect_no_error(sc)
  expect_type(sc, "character")
  expect_length(sc, 1)
})

test_that("coev_make_stancode() works with multiPhylo object", {
  # simulate data
  withr::with_seed(1, {
    n <- 10
    tree <- c(ape::rcoal(n), ape::rcoal(n))
    d <- data.frame(
      id = tree[[1]]$tip.label,
      x = rbinom(n, size = 1, prob = 0.5),
      y = rbinom(n, size = 1, prob = 0.5)
    )
  })
  # get stan code
  sc <- coev_make_stancode(
    data = d,
    variables = list(
      x = "bernoulli_logit",
      y = "bernoulli_logit"
    ),
    id = "id",
    tree = tree
  )
  # runs without error
  expect_no_error(sc)
  # syntactically correct
  expect_true(
    cmdstanr::cmdstan_model(
      stan_file = cmdstanr::write_stan_file(sc),
      compile = FALSE
    )$check_syntax(quiet = TRUE)
  )
})

test_that("coev_make_stancode() works when constraining Q offdiag to zero", {
  # simulate data
  withr::with_seed(1, {
    n <- 10
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = tree$tip.label,
      x = rbinom(n, size = 1, prob = 0.5),
      y = rbinom(n, size = 1, prob = 0.5)
    )
  })
  # get stan code
  sc <- coev_make_stancode(
    data = d,
    variables = list(
      x = "bernoulli_logit",
      y = "bernoulli_logit"
    ),
    id = "id",
    tree = tree,
    estimate_Q_offdiag = FALSE # Q off diagonal constrained to zero
  )
  # runs without error
  expect_no_error(sc)
  # syntactically correct
  expect_true(
    cmdstanr::cmdstan_model(
      stan_file = cmdstanr::write_stan_file(sc),
      compile = FALSE
    )$check_syntax(quiet = TRUE)
  )
})

test_that("GP covariance kernels produce syntactically correct Stan code", {
  # simulate data
  withr::with_seed(1, {
    n <- 10
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = tree$tip.label,
      x = rbinom(n, size = 1, prob = 0.5),
      y = rbinom(n, size = 1, prob = 0.5)
    )
    dist_mat <- as.matrix(dist(rnorm(n)))
    rownames(dist_mat) <- colnames(dist_mat) <- tree$tip.label
  })
  # get stan code
  sc1 <- coev_make_stancode(
    data = d,
    variables = list(
      x = "bernoulli_logit",
      y = "bernoulli_logit"
    ),
    id = "id",
    tree = tree,
    dist_mat = dist_mat,
    dist_cov = "exp_quad"
  )
  sc2 <- coev_make_stancode(
    data = d,
    variables = list(
      x = "bernoulli_logit",
      y = "bernoulli_logit"
    ),
    id = "id",
    tree = tree,
    dist_mat = dist_mat,
    dist_cov = "exponential"
  )
  sc3 <- coev_make_stancode(
    data = d,
    variables = list(
      x = "bernoulli_logit",
      y = "bernoulli_logit"
    ),
    id = "id",
    tree = tree,
    dist_mat = dist_mat,
    dist_cov = "matern32"
  )
  # runs without error
  expect_no_error(sc1)
  expect_no_error(sc2)
  expect_no_error(sc3)
  # syntactically correct
  expect_true(
    cmdstanr::cmdstan_model(
      stan_file = cmdstanr::write_stan_file(sc1),
      compile = FALSE
    )$check_syntax(quiet = TRUE)
  )
  expect_true(
    cmdstanr::cmdstan_model(
      stan_file = cmdstanr::write_stan_file(sc2),
      compile = FALSE
    )$check_syntax(quiet = TRUE)
  )
  expect_true(
    cmdstanr::cmdstan_model(
      stan_file = cmdstanr::write_stan_file(sc3),
      compile = FALSE
    )$check_syntax(quiet = TRUE)
  )
})

test_that("coev_make_stancode() works with gamma_log distribution", {
  # simulate data
  withr::with_seed(1, {
    n <- 20
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = tree$tip.label,
      x = rgamma(n, shape = 1, rate = 1),
      y = rgamma(n, shape = 1, rate = 1)
    )
  })
  # default priors should work
  expect_no_error(
    coev_make_stancode(
      data = d,
      variables = list(
        x = "gamma_log",
        y = "gamma_log"
      ),
      id = "id",
      tree = tree
    )
  )
  # setting a prior for shape parameter should also work
  expect_no_error(
    coev_make_stancode(
      data = d,
      variables = list(
        x = "gamma_log",
        y = "gamma_log"
      ),
      id = "id",
      tree = tree,
      prior = list(shape = "gamma(0.05, 0.05)")
    )
  )
})
