test_that("coev_make_standata() produces expected errors", {
  # simulate data
  withr::with_seed(1, {
    n <- 20
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = tree$tip.label,
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
      "Response distributions other than 'bernoulli_logit', 'ordered_logistic'",
      ", 'poisson_softmax', 'normal', and 'lognormal' are not yet supported."
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
    "Variables following the 'bernoulli_logit' response distribution must be integers with values of 0/1 in the data."
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
    "Variables following the 'ordered_logistic' response distribution must be ordered factors in the data."
  )
  expect_error(
    coev_make_standata(
      data = d,
      variables = list(
        y = "poisson_softmax", # not integer >= 0
        z = "poisson_softmax"
      ),
      id = "id",
      tree = tree
    ),
    "Variables following the 'poisson_softmax' response distribution must be integers greater than or equal to zero in the data."
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
    "Variables following the 'normal' response distribution must be numeric in the data."
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
    "Argument 'tree' must be an phylogenetic tree object of class phylo."
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
    {
      d2 <- d; d2$y[1] <- NA
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
    "Coevolving variables in the data must not contain NAs."
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
      effects_mat = matrix(TRUE, dimnames = list("fail","fail")) # invalid row/col names
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
      effects_mat = matrix(c(T,T,T,F), ncol = 2, nrow = 2, byrow = TRUE,
                           dimnames = list(c("x","y"),c("x","y"))) # autoregressive effect = FALSE
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
      dist_mat = matrix(rep(1, 100), nrow = 10) # matrix symmetric but diagonal not zero
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
    "Row and column names for argument 'dist_mat' do not match tree tip labels exactly."
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
      "use only the following names: 'alpha', 'b', 'sigma', 'eta_anc', ",
      "'c', 'sigma_dist', and 'rho_dist'"
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
      prior = list(alpha = "normal(0,2)", alpha = "normal(0,2)") # duplicate names
    ),
    "Argument 'prior' contains duplicate names."
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
      x = rbinom(n, size = 1, prob = 0.5),
      y = ordered(sample(1:4, size = n, replace = TRUE))
    )
  })
  # make stan data
  sd1 <-
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
  expect_no_error(sd1)
  expect_type(sd1, "list")
  expect_equal(names(sd1), c("N", "J", "N_seg", "node_seq", "parent", "ts",
                            "tip", "effects_mat", "num_effects", "y1", "y2",
                            "prior_only"))
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
  expect_equal(names(sd2), c("N", "J", "N_seg", "node_seq", "parent", "ts",
                             "tip", "effects_mat", "num_effects", "y1", "y2",
                             "dist_mat", "prior_only"))
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
  expect_equal(names(sd3), c("N", "J", "N_seg", "node_seq", "parent", "ts",
                             "tip", "effects_mat", "num_effects",
                             "y1", "y2", "prior_only"))
  expect_equal(sd3$prior_only, 1)
})
