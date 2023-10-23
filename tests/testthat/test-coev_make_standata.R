test_that("coev_make_standata() produces expected errors", {
  # simulate data
  withr::with_seed(1, {
    n <- 20
    tree <- ape::rtree(n)
    d <- data.frame(
      id = tree$tip.label,
      x = rbinom(n, size = 1, prob = 0.5),
      y = ordered(sample(1:4, size = n, replace = TRUE))
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
        x = "poisson", # incorrect response distribution
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree
    ),
    "Response distributions other than 'bernoulli_logit' and 'ordered_logistic' are not yet supported."
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
    "Argument 'id' must be an phylogenetic tree object of class phylo."
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
})


test_that("coev_make_standata() returns a list with correct names for stan", {
  # simulate data
  withr::with_seed(1, {
    n <- 20
    tree <- ape::rtree(n)
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
  # expect list with correct names
  expect_no_error(sd)
  expect_type(sd, "list")
  expect_equal(names(sd), c("N", "J", "N_seg", "node_seq", "parent", "ts", "y"))
})
