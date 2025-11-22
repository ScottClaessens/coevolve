test_that("coev_make_standata() produces expected errors", {
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
    coev_make_standata(
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
    coev_make_standata(
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
    coev_make_standata(
      data = d,
      variables = "testing", # not a named list
      id = "id",
      tree = tree
    ),
    "Argument 'variables' is not a named list.",
    fixed = TRUE
  )
  expect_error(
    coev_make_standata(
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
    coev_make_standata(
      data = d,
      variables = list(
        x = "test", # incorrect response distribution
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
    coev_make_standata(
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
    coev_make_standata(
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
    coev_make_standata(
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
    coev_make_standata(
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
    coev_make_standata(
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
    coev_make_standata(
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
    coev_make_standata(
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
    coev_make_standata(
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
    coev_make_standata(
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
    coev_make_standata(
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
    coev_make_standata(
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
    coev_make_standata(
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
    coev_make_standata(
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
    coev_make_standata(
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
    coev_make_standata(
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
      coev_make_standata(
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
      coev_make_standata(
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
      coev_make_standata(
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
      coev_make_standata(
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
    coev_make_standata(
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
    coev_make_standata(
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
    coev_make_standata(
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
    coev_make_standata(
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
    coev_make_standata(
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
    coev_make_standata(
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
    coev_make_standata(
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
    coev_make_standata(
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
    coev_make_standata(
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
    coev_make_standata(
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
    coev_make_standata(
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
    coev_make_standata(
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
    coev_make_standata(
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
    coev_make_standata(
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
    coev_make_standata(
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
    coev_make_standata(
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
    coev_make_standata(
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
    coev_make_standata(
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
    coev_make_standata(
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
    coev_make_standata(
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
    coev_make_standata(
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
    coev_make_standata(
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
    coev_make_standata(
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

test_that("coev_make_standata() returns a list with correct names for Stan", {
  # simulate data
  withr::with_seed(1, {
    n <- 20
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = tree$tip.label,
      u = rgamma(n, shape = 1, rate = 1),
      v = as.integer(rnbinom(n, mu = 4, size = 1)),
      w = rnorm(n),
      x = rbinom(n, size = 1, prob = 0.5),
      y = ordered(sample(1:4, size = n, replace = TRUE)),
      z = rpois(n, 3)
    )
  })
  # make stan data
  sd1 <-
    coev_make_standata(
      data = d,
      variables = list(
        u = "gamma_log",
        v = "negative_binomial_softplus",
        w = "normal",
        x = "bernoulli_logit",
        y = "ordered_logistic",
        z = "poisson_softplus"
      ),
      id = "id",
      tree = tree
    )
  # expect list with correct names and prior_only = 0
  expect_no_error(sd1)
  expect_type(sd1, "list")
  expect_equal(
    names(sd1),
    c("N_tips", "N_tree", "N_obs", "J", "N_seg", "node_seq", "parent", "ts",
      "tip", "effects_mat", "num_effects", "y", "miss", "tip_id",
      "N_unique_lengths", "unique_lengths", "length_index", "tip_to_seg",
      "prior_only")
  )
  expect_equal(sd1$prior_only, 0)
  # include distance matrix
  withr::with_seed(1, {
    dist_mat <- as.matrix(dist(rnorm(n)))
    rownames(dist_mat) <- colnames(dist_mat) <- tree$tip.label
  })
  sd2 <-
    coev_make_standata(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree,
      dist_mat = dist_mat
    )
  # expect list with correct names and prior_only = 0
  expect_no_error(sd2)
  expect_type(sd2, "list")
  expect_equal(
    names(sd2),
    c("N_tips", "N_tree", "N_obs", "J", "N_seg", "node_seq", "parent", "ts",
      "tip", "effects_mat", "num_effects", "y", "miss", "tip_id",
      "N_unique_lengths", "unique_lengths", "length_index", "tip_to_seg",
      "dist_mat", "prior_only")
  )
  expect_equal(sd2$prior_only, 0)
  # set prior only
  sd3 <-
    coev_make_standata(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree,
      prior_only = TRUE
    )
  expect_no_error(sd3)
  expect_type(sd3, "list")
  expect_equal(
    names(sd3),
    c("N_tips", "N_tree", "N_obs", "J", "N_seg", "node_seq", "parent", "ts",
      "tip", "effects_mat", "num_effects", "y", "miss", "tip_id",
      "N_unique_lengths", "unique_lengths", "length_index", "tip_to_seg",
      "prior_only")
  )
  expect_equal(sd3$prior_only, 1)
})

#' @srrstats {BS2.1a} Testing that input data is dimensionally commensurate
test_that("coev_make_standata() produces dimensionally commensurate data", {
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
  # make stan data
  sd <-
    coev_make_standata(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree
    )
  # test for dimensionally commensurate data list
  expect_true(is.matrix(sd$y))
  expect_equal(length(sd$y[, 1]), length(sd$y[, 2]))
  expect_equal(nrow(sd$y), n)
  expect_equal(ncol(sd$y), 2)
})

test_that("coev_make_standata() works with missing data", {
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
  # some data are missing
  # row 1 missing both x and y
  # row 2 missing x only
  # row 3 missing y only
  d$x[c(1, 2)] <- NA
  d$y[c(1, 3)] <- NA
  # make stan data with missing values to be imputed
  # should return expected message with missing data
  expect_message(
    {sd1 <- coev_make_standata(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree,
      complete_cases = FALSE
    )},
    paste0(
      "Note: Missing values (NAs) detected. These values have been ",
      "included in the Stan data list and will be imputed by the model. Set ",
      "complete_cases = TRUE to exclude taxa with missing values."
    ),
    fixed = TRUE
  )
  # make stan data with listwise deletion
  sd2 <-
    coev_make_standata(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree,
      complete_cases = TRUE
    )
  # runs without error
  expect_no_error(sd1)
  expect_no_error(sd2)
  # N_tips correct
  expect_equal(sd1$N_tips, 20)
  expect_equal(sd2$N_tips, 17)
  expect_equal(nrow(sd1$y), 20)
  expect_equal(nrow(sd2$y), 17)
  expect_equal(nrow(sd1$miss), 20)
  expect_equal(nrow(sd2$miss), 17)
  # missing matrix correct
  miss1 <- as.matrix(ifelse(is.na(d[, c("x", "y")]), 1, 0))
  miss2 <-
    matrix(0, nrow = sum(apply(d, 1, function(x) all(!is.na(x)))), ncol = 2)
  rownames(miss1) <- NULL
  colnames(miss2) <- c("x", "y")
  expect_identical(sd1$miss, miss1)
  expect_identical(sd2$miss, miss2)
  # dataset correct
  d1 <- d[, c("x", "y")]
  d2 <- d[, c("x", "y")]
  d1$x <- ifelse(is.na(d1$x), -9999, d1$x)
  d1$y <- ifelse(is.na(d1$y), -9999, d1$y)
  d2 <- d2[apply(d2, 1, function(x) all(!is.na(x))), ]
  d2$y <- as.numeric(d2$y)
  d1 <- as.matrix(d1, rownames.force = FALSE)
  d2 <- as.matrix(d2, rownames.force = FALSE)
  expect_identical(sd1$y, d1)
  expect_identical(sd2$y, d2)
})

test_that("coev_make_standata() works with repeated observations", {
  # simulate data with repeated observations
  withr::with_seed(1, {
    n <- 10
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = rep(tree$tip.label, each = 10),
      x = rnorm(n * 10),
      y = rnorm(n * 10)
    )
  })
  # make stan data
  sd <-
    coev_make_standata(
      data = d,
      variables = list(
        x = "normal",
        y = "normal"
      ),
      id = "id",
      tree = tree
    )
  # run without error
  expect_no_error(sd)
  # expect correct output
  expect_true(sd$N_tips == 10)
  expect_true(sd$N_obs == 100)
  expect_true(sum(sd$tip) == 10)
  expect_true(nrow(sd$y) == 100)
  expect_true(nrow(sd$miss) == 100)
  expect_true(identical(sd$tip_id, match(d$id, tree$tip.label)))
})

test_that("coev_make_stancode() scales data correctly", {
  # simulate data
  withr::with_seed(1, {
    n <- 20
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = tree$tip.label,
      w = rnorm(n, 10, 0.5),
      x = as.integer(rnbinom(n, mu = 4, size = 1)),
      y = rnorm(n, 5, 2),
      z = rlnorm(n, 1, 1.5)
    )
  })
  # check stan data is the same as real data when scale is FALSE
  expect_warning(
    {sd1 <-
      coev_make_standata(
        data = d,
        variables = list(
          w = "normal",
          x = "negative_binomial_softplus"
        ),
        id = "id",
        tree = tree,
        scale = FALSE
      )},
    paste0(
      "When scale is FALSE, continuous variables are left unstandardised in ",
      "the Stan data list. Users should take care to set sensible priors ",
      "for model fitting, rather than use default priors."
    )
  )
  expect_identical(sd1$y[, "w"], as.numeric(d$w))
  expect_identical(sd1$y[, "x"], as.numeric(d$x))
  # check stan data is correctly standardised when scale = TRUE
  expect_no_warning(
    {
      sd2 <-
        coev_make_standata(
          data = d,
          variables = list(
            w = "normal",
            x = "negative_binomial_softplus"
          ),
          id = "id",
          tree = tree,
          scale = TRUE
        )
    }
  )
  expect_identical(sd2$y[, "w"], as.numeric(scale(d$w)))
  expect_identical(sd2$y[, "x"], as.numeric(d$x))
})

test_that("coev_make_standata() works with tibbles", {
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
  # make stan data
  sd <-
    coev_make_standata(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree
    )
  # expect list with correct names and prior_only = 0
  expect_no_error(sd)
  expect_type(sd, "list")
  expect_equal(
    names(sd),
    c("N_tips", "N_tree", "N_obs", "J", "N_seg", "node_seq", "parent", "ts",
      "tip", "effects_mat", "num_effects", "y", "miss", "tip_id",
      "N_unique_lengths", "unique_lengths", "length_index", "tip_to_seg",
      "prior_only")
  )
  expect_equal(sd$prior_only, 0)
})

test_that("coev_make_standata() works with multiPhylo object", {
  # simulate data with repeated observations
  withr::with_seed(1, {
    n <- 10
    tree <- c(ape::rcoal(n), ape::rcoal(n))
    d <- data.frame(
      id = tree[[1]]$tip.label,
      x = rnorm(n),
      y = rnorm(n)
    )
  })
  # make stan data
  sd <-
    coev_make_standata(
      data = d,
      variables = list(
        x = "normal",
        y = "normal"
      ),
      id = "id",
      tree = tree
    )
  # run without error
  expect_no_error(sd)
  # expect correct output
  expect_true(sd$N_tree == 2)
  expect_true(nrow(sd$node_seq) == 2)
  expect_true(nrow(sd$parent) == 2)
  expect_true(nrow(sd$ts) == 2)
  expect_true(nrow(sd$tip) == 2)
})

test_that("coev_make_standata() works with measurement error", {
  # simulate data
  withr::with_seed(1, {
    n <- 20
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = tree$tip.label,
      x = rnorm(n),
      y = rnorm(n),
      z = rbinom(n, 1, 0.5),
      x_se = 0, # no measurement error
      y_se = rexp(20, 5)
    )
  })
  # make stan data with measurement error when scale = TRUE
  sd1 <-
    coev_make_standata(
      data = d,
      variables = list(
        x = "normal",
        y = "normal",
        z = "bernoulli_logit"
      ),
      id = "id",
      tree = tree,
      measurement_error = list(
        x = "x_se",
        y = "y_se"
      )
    )
  # make stan data with measurement error when scale is FALSE
  sd2 <-
    suppressWarnings(
      coev_make_standata(
        data = d,
        variables = list(
          x = "normal",
          y = "normal",
          z = "bernoulli_logit"
        ),
        id = "id",
        tree = tree,
        measurement_error = list(
          x = "x_se",
          y = "y_se"
        ),
        scale = FALSE
      )
    )
  # check squared standard errors are scaled correctly in stan data
  expect_identical(sd1$se[, 1], (d$x_se / sd(d$x))^2)
  expect_identical(sd1$se[, 2], (d$y_se / sd(d$y))^2)
  expect_identical(sd2$se[, 1], d$x_se^2)
  expect_identical(sd2$se[, 2], d$y_se^2)
  # se set to zero for non-normal variables
  expect_true(all(sd1$se[, 3] == 0))
  expect_true(all(sd2$se[, 3] == 0))
})

test_that("coev_make_standata() does correct caching for Stan", {
  # simulate data
  withr::with_seed(1, {
    n <- 20
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = tree$tip.label,
      x = rnorm(n),
      y = rnorm(n)
    )
  })
  # make stan data list
  sd <-
    coev_make_standata(
      data = d,
      variables = list(
        x = "normal",
        y = "normal"
      ),
      id = "id",
      tree = tree
    )
  # test caching computations
  lengths <- as.vector(sd$ts)[-1]
  unique_lengths <- sort(unique(lengths))
  expect_equal(sd$N_unique_lengths, length(unique_lengths))
  expect_equal(sd$unique_lengths, unique_lengths)
  expect_equal(as.vector(sd$length_index), c(0, match(lengths, unique_lengths)))
  stan_tip_to_seg <- matrix(0L, 1L, length(tree$tip.label))
  for (seg in 1:sd$N_seg) {
    if (sd$tip[1, seg] == 1) {
      tip_node <- sd$node_seq[1, seg]
      if (tip_node >= 1 && tip_node <= length(tree$tip.label)) {
        stan_tip_to_seg[1, tip_node] <- seg
      }
    }
  }
  expect_equal(sd$tip_to_seg, stan_tip_to_seg)
})
