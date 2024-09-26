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
    "Argument 'data' must be coercible to a data.frame."
  )
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
    "Argument 'data' does not contain observations."
  )
  expect_error(
    coev_make_standata(
      data = d,
      variables = "testing", # not a named list
      id = "id",
      tree = tree
    ),
    "Argument 'variables' is not a named list."
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
    "Some variable names are not valid column names in the data."
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
      "'ordered_logistic', 'poisson_softplus', 'normal', 'student_t', ",
      "'lognormal', and 'negative_binomial_softplus' are not yet supported."
    )
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
    "Must be at least two coevolving variables."
  )
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
      )
  )
  expect_error(
    coev_make_standata(
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
    coev_make_standata(
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
    coev_make_standata(
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
    "All trees in 'tree' argument must include branch lengths."
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
    "All trees in 'tree' argument must be rooted."
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
    "The id variable in the data does not match tree tip labels exactly."
  )
  expect_error(
    {
      d2 <- d; d2$id[1] <- NA
      tree2 <- tree; tree2$tip.label[1] <- NA
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
    "The id variable in the data must not contain NAs."
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
    "Argument 'effects_mat' must be a matrix."
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
    "Argument 'effects_mat' must be a boolean matrix."
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
    "Argument 'effects_mat' does not have valid row or column names."
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
      effects_mat = matrix(TRUE, dimnames = list("fail","fail")) # invalid names
    ),
    paste0(
      "Row and column names for argument 'effects_mat' do not match ",
      "variable names exactly."
    )
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
      effects_mat = matrix(c(T,T,T,F), ncol = 2, nrow = 2, byrow = TRUE,
                           dimnames = list(c("x","y"),c("x","y")))
    ),
    "Argument 'effects_mat' must specify TRUE for all autoregressive effects."
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
    "Argument 'dist_mat' must be a matrix."
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
    "Argument 'dist_mat' must be a numeric matrix."
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
    "Argument 'dist_mat' must be a symmetric matrix."
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
    "Argument 'dist_mat' must have zeroes on the diagonal of the matrix."
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
    "Argument 'dist_mat' does not have valid row or column names."
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
      )
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
    "Argument 'prior' is not a list."
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
    "Argument 'prior' is not a named list."
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
      "'Q_diag', 'c', 'phi', 'nu', 'sigma_dist', 'rho_dist', 'sigma_group', ",
      "and 'L_group'"
    )
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
    "Argument 'prior_only' is not logical."
  )
})

test_that("coev_make_standata() returns a list with correct names for Stan", {
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
  # make stan data
  sd1 <-
    coev_make_standata(
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
  # expect list with correct names and prior_only = 0
  expect_no_error(sd1)
  expect_type(sd1, "list")
  expect_equal(
    names(sd1),
    c("N_tips", "N_tree", "N_obs", "J", "N_seg", "node_seq", "parent", "ts",
      "tip", "effects_mat", "num_effects", "y", "miss", "tip_id", "prior_only")
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
      "tip", "effects_mat", "num_effects", "y", "miss", "tip_id", "dist_mat",
      "prior_only")
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
      "tip", "effects_mat", "num_effects", "y", "miss", "tip_id", "prior_only")
  )
  expect_equal(sd3$prior_only, 1)
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
  d$x[c(1,2)] <- NA
  d$y[c(1,3)] <- NA
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
  # runs without error
  expect_no_error(sd)
  # N_tips = 19 because one row removed
  expect_equal(sd$N_tips, 19)
  expect_equal(nrow(sd$y), 19)
  expect_equal(nrow(sd$miss), 19)
  # missing matrix correct
  d <- d[-1,c("x","y")]
  miss <- as.matrix(ifelse(is.na(d), 1, 0))
  rownames(miss) <- NULL
  expect_identical(sd$miss, miss)
  # dataset correct
  d$x <- ifelse(is.na(d$x), -9999, d$x)
  d$y <- ifelse(is.na(d$y), -9999, d$y)
  d <- as.matrix(d, rownames.force = FALSE)
  expect_identical(sd$y, d)
})

test_that("coev_make_standata() works with repeated observations", {
  # simulate data with repeated observations
  withr::with_seed(1, {
    n <- 10
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = rep(tree$tip.label, each = 10),
      x = rnorm(n*10),
      y = rnorm(n*10)
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
  # check stan data is the same as real data when scale = FALSE
  expect_warning(
    {sd1 <-
      coev_make_standata(
        data = d,
        variables = list(
          w = "normal",
          x = "negative_binomial_softplus",
          y = "student_t",
          z = "lognormal"
        ),
        id = "id",
        tree = tree,
        scale = FALSE
      )},
    paste0(
      "When scale = FALSE, continuous and positive real variables are left ",
      "unstandardised in the Stan data list. Users should take care to set ",
      "sensible priors for model fitting, rather than use default priors."
      )
    )
  expect_identical(sd1$y[,"w"], as.numeric(d$w))
  expect_identical(sd1$y[,"x"], as.numeric(d$x))
  expect_identical(sd1$y[,"y"], as.numeric(d$y))
  expect_identical(sd1$y[,"z"], as.numeric(d$z))
  # check stan data is correctly standardised when scale = TRUE
  expect_no_warning(
    {sd2 <-
      coev_make_standata(
        data = d,
        variables = list(
          w = "normal",
          x = "negative_binomial_softplus",
          y = "student_t",
          z = "lognormal"
        ),
        id = "id",
        tree = tree,
        scale = TRUE
      )}
    )
  expect_identical(sd2$y[,"w"], as.numeric(scale(d$w)))
  expect_identical(sd2$y[,"x"], as.numeric(d$x))
  expect_identical(sd2$y[,"y"], as.numeric(scale(d$y)))
  expect_identical(sd2$y[,"z"], as.numeric(exp(scale(log(d$z)))))
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
      "tip", "effects_mat", "num_effects", "y", "miss", "tip_id", "prior_only")
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
