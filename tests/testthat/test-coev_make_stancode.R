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
    "Argument 'data' must be coercible to a data.frame."
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
    "Argument 'data' does not contain observations."
  )
  expect_error(
    coev_make_stancode(
      data = d,
      variables = "testing", # not a named list
      id = "id",
      tree = tree
    ),
    "Argument 'variables' is not a named list."
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
    "Some variable names are not valid column names in the data."
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
      "'ordered_logistic', 'poisson_softplus', 'normal', 'student_t', ",
      "'lognormal', and 'negative_binomial_softplus' are not yet supported."
    )
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
    "Must be at least two coevolving variables."
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
      )
  )
  expect_error(
    coev_make_stancode(
      data = d,
      variables = list(
        w = "student_t", # not numeric
        y = "student_t"
      ),
      id = "id",
      tree = tree
    ),
    paste0(
      "Variables following the 'student_t' response distribution ",
      "must be numeric in the data."
    )
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
    "Argument 'id' must be a character of length one."
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
    "Argument 'id' is not a valid column name in the data."
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
    "Argument 'tree' must be an phylogenetic tree object of class phylo."
  )
  expect_error(
    coev_make_stancode(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = ape::rtree(n, br = NULL)
    ),
    "Argument 'tree' does not include branch lengths."
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
    "The id variable in the data does not match tree tip labels exactly."
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
    "The id variable in the data must not contain NAs."
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
    "Argument 'effects_mat' must be a matrix."
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
    "Argument 'effects_mat' must be a boolean matrix."
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
    "Argument 'effects_mat' does not have valid row or column names."
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
    )
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
    "Argument 'effects_mat' must specify TRUE for all autoregressive effects."
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
    "Argument 'dist_mat' must be a matrix."
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
    "Argument 'dist_mat' must be a numeric matrix."
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
    "Argument 'dist_mat' must be a symmetric matrix."
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
    "Argument 'dist_mat' must have zeroes on the diagonal of the matrix."
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
    "Argument 'dist_mat' does not have valid row or column names."
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
      )
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
    "Argument 'prior' is not a list."
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
    "Argument 'prior' is not a named list."
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
      "'Q_diag', 'c', 'phi', 'nu', 'sigma_dist', 'rho_dist', 'sigma_group', ",
      "and 'L_group'"
    )
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
    "Argument 'prior' contains duplicate names."
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
    "Argument 'scale' is not logical."
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
    "Argument 'prior_only' is not logical."
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
      u = rnorm(n),
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
          u = "student_t",
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
        b         = "normal(0, 2)",
        eta_anc   = "normal(0, 2)",
        A_offdiag = "normal(0, 2)",
        A_diag    = "normal(0, 2)",
        Q_diag    = "normal(0, 2)"
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
