test_that("coev_fit() produces expected errors", {
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
    coev_fit(
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
    coev_fit(
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
    coev_fit(
      data = d,
      variables = "testing", # not a named list
      id = "id",
      tree = tree
    ),
    "Argument 'variables' is not a named list."
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
    "Some variable names are not valid column names in the data."
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
    "Response distributions other than 'bernoulli_logit', 'ordered_logistic', 'poisson_log', and 'normal' are not yet supported."
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
    "Must be at least two coevolving variables."
  )
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
    "Variables following the 'bernoulli_logit' response distribution must be integers with values of 0/1 in the data."
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
    "Variables following the 'ordered_logistic' response distribution must be ordered factors in the data."
  )
  expect_error(
    coev_fit(
      data = d,
      variables = list(
        y = "poisson_log", # not integer >= 0
        z = "poisson_log"
      ),
      id = "id",
      tree = tree
    ),
    "Variables following the 'poisson_log' response distribution must be integers greater than or equal to zero in the data."
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
    "Variables following the 'normal' response distribution must be numeric in the data."
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
    "Argument 'id' must be a character of length one."
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
    "Argument 'id' is not a valid column name in the data."
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
    "Argument 'id' must be an phylogenetic tree object of class phylo."
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
    "The id variable in the data does not match tree tip labels exactly."
  )
  expect_error(
    {
      d2 <- d; d2$id[1] <- NA
      tree2 <- tree; tree2$tip.label[1] <- NA
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
    "The id variable in the data must not contain NAs."
  )
  expect_error(
    {
      d2 <- d; d2$y[1] <- NA
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
    "Coevolving variables in the data must not contain NAs."
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
    "Argument 'dist_mat' must be a matrix."
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
    "Argument 'dist_mat' must be a numeric matrix."
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
    "Argument 'dist_mat' must be a symmetric matrix."
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
      dist_mat = matrix(rep(1, 100), nrow = 10) # matrix symmetric but diagonal not zero
    ),
    "Argument 'dist_mat' must have zeroes on the diagonal of the matrix."
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
      dist_mat = matrix(rep(0, 100), nrow = 10) # matrix row/col names do not match tips
    ),
    "Row and column names for argument 'dist_mat' do not match tree tip labels exactly."
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
    "Argument 'prior' is not a list."
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
    "Argument 'prior' is not a named list."
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
      "use only the following names: 'alpha', 'b', 'sigma', 'eta_anc', ",
      "'c', 'sigma_dist', and 'rho_dist'"
    )
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
      prior = list(alpha = "normal(0,2)", alpha = "normal(0,2)") # duplicate names
    ),
    "Argument 'prior' contains duplicate names."
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
    "Argument 'prior_only' is not logical."
  )
})

test_that("coev_fit() fits the model without error", {
  # simulate data
  withr::with_seed(1, {
    n <- 5
    tree <- ape::rcoal(n)
    d <- data.frame(
      id = tree$tip.label,
      w = rnorm(n),
      x = rbinom(n, size = 1, prob = 0.5),
      y = ordered(sample(1:4, size = n, replace = TRUE)),
      z = rpois(n, 3)
    )
  })
  # model without distance matrix
  m1 <- coev_fit(
    data = d,
    variables = list(
      w = "normal",
      x = "bernoulli_logit",
      y = "ordered_logistic",
      z = "poisson_log"
    ),
    id = "id",
    tree = tree,
    parallel_chains = 4,
    iter_warmup = 100,
    iter_sampling = 100,
    seed = 1
  )
  # model with distance matrix
  dist_mat <- as.matrix(dist(rnorm(n)))
  rownames(dist_mat) <- colnames(dist_mat) <- d$id
  m2 <- coev_fit(
    data = d,
    variables = list(
      x = "bernoulli_logit",
      y = "ordered_logistic"
    ),
    id = "id",
    tree = tree,
    dist_mat = dist_mat,
    parallel_chains = 4,
    iter_warmup = 100,
    iter_sampling = 100,
    seed = 1
  )
  # prior only model
  m3 <- coev_fit(
    data = d,
    variables = list(
      x = "bernoulli_logit",
      y = "ordered_logistic"
    ),
    id = "id",
    tree = tree,
    parallel_chains = 4,
    iter_warmup = 100,
    iter_sampling = 100,
    seed = 1
  )
  # expect no errors for model fitting or summaries
  expect_no_error(m1)
  expect_no_error(m2)
  expect_no_error(m3)
  expect_no_error(summary(m1))
  expect_no_error(summary(m2))
  expect_no_error(summary(m3))
  expect_output(print(m1))
  expect_output(print(m2))
  expect_output(print(m3))
  expect_output(print(summary(m1)))
  expect_output(print(summary(m2)))
  expect_output(print(summary(m3)))
  # expect error if prob for summary is outside of range 0 - 1
  expect_error(summary(m1, prob = -0.01))
  expect_error(summary(m2, prob = -0.01))
  expect_error(summary(m3, prob = -0.01))
  expect_error(summary(m1, prob =  1.01))
  expect_error(summary(m2, prob =  1.01))
  expect_error(summary(m3, prob =  1.01))
})
