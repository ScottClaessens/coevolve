test_that("coev_make_stancode() throws error if passed data that is not coercible to data.frame", {
  # simulate data
  set.seed(1)
  n <- 20
  tree <- ape::rtree(n)
  d <- log # e.g. should not accept functions
  # expect error
  expect_error(
    coev_make_stancode(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree
    )
  )
})

test_that("coev_make_stancode() throws error if passed data that does not contain observations", {
  # simulate data
  set.seed(1)
  n <- 20
  tree <- ape::rtree(n)
  d <- data.frame()
  # expect error
  expect_error(
    coev_make_stancode(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree
    )
  )
})

test_that("coev_make_stancode() throws error if variables argument is not a named list", {
  # simulate data
  set.seed(1)
  n <- 20
  tree <- ape::rtree(n)
  d <- data.frame()
  # expect error
  expect_error(
    coev_make_stancode(
      data = d,
      variables = "testing", # not a named list
      id = "id",
      tree = tree
    )
  )
})

test_that("coev_make_stancode() throws error if variables do not refer to columns in data", {
  # simulate data
  set.seed(1)
  n <- 20
  tree <- ape::rtree(n)
  d <- data.frame(
    id = tree$tip.label,
    x = rbinom(n, size = 1, prob = 0.5),
    y = ordered(sample(1:4, size = n, replace = TRUE))
  )
  # expect error
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
    )
  )
})

test_that("coev_make_stancode() throws error if response distributions are not valid", {
  # simulate data
  set.seed(1)
  n <- 20
  tree <- ape::rtree(n)
  d <- data.frame(
    id = tree$tip.label,
    x = rbinom(n, size = 1, prob = 0.5),
    y = ordered(sample(1:4, size = n, replace = TRUE))
  )
  # expect error
  expect_error(
    coev_make_stancode(
      data = d,
      variables = list(
        x = "poisson", # incorrect response distribution
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree
    )
  )
})

test_that("coev_make_stancode() throws error if only one variable declared", {
  # simulate data
  set.seed(1)
  n <- 20
  tree <- ape::rtree(n)
  d <- data.frame(
    id = tree$tip.label,
    x = rbinom(n, size = 1, prob = 0.5),
    y = ordered(sample(1:4, size = n, replace = TRUE))
  )
  # expect error
  expect_error(
    coev_make_stancode(
      data = d,
      variables = list(
        # only x declared
        x = "bernoulli_logit"
      ),
      id = "id",
      tree = tree
    )
  )
})

test_that("coev_make_stancode() throws error if any bernoulli_logit variables are not binary integers", {
  # simulate data
  set.seed(1)
  n <- 20
  tree <- ape::rtree(n)
  d <- data.frame(
    id = tree$tip.label,
    # x is not a binary integer, but y is
    x = as.integer(sample(1:4, size = n, replace = TRUE)),
    y = as.integer(sample(0:1, size = n, replace = TRUE))
  )
  # expect error
  expect_error(
    coev_make_stancode(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "bernoulli_logit"
      ),
      id = "id",
      tree = tree
    )
  )
})

test_that("coev_make_stancode() throws error if any ordered_logistic variables are not ordered factors", {
  # simulate data
  set.seed(1)
  n <- 20
  tree <- ape::rtree(n)
  d <- data.frame(
    id = tree$tip.label,
    # x is not an ordered factor, but y is
    x = as.integer(sample(1:4, size = n, replace = TRUE)),
    y = ordered(sample(1:3, size = n, replace = TRUE))
  )
  # expect error
  expect_error(
    coev_make_stancode(
      data = d,
      variables = list(
        x = "ordered_logistic",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree
    )
  )
})

test_that("coev_make_stancode() throws error if id is not a character of length one", {
  # simulate data
  set.seed(1)
  n <- 20
  tree <- ape::rtree(n)
  d <- data.frame(
    id = tree$tip.label,
    x = rbinom(n, size = 1, prob = 0.5),
    y = ordered(sample(1:4, size = n, replace = TRUE))
  )
  # expect error
  expect_error(
    coev_make_stancode(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = c(1, 2), # not character of length one
      tree = tree
    )
  )
})

test_that("coev_make_stancode() throws error if id is not a valid column name in data", {
  # simulate data
  set.seed(1)
  n <- 20
  tree <- ape::rtree(n)
  d <- data.frame(
    id = tree$tip.label,
    x = rbinom(n, size = 1, prob = 0.5),
    y = ordered(sample(1:4, size = n, replace = TRUE))
  )
  # expect error
  expect_error(
    coev_make_stancode(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "testing", # incorrect id
      tree = tree
    )
  )
})

test_that("coev_make_stancode() throws error if id in data does not match tip labels", {
  # simulate data
  set.seed(1)
  n <- 20
  tree <- ape::rtree(n)
  d <- data.frame(
    # id variable does not match tip labels exactly
    id = sample(letters, size = n, replace = TRUE),
    x = rbinom(n, size = 1, prob = 0.5),
    y = ordered(sample(1:4, size = n, replace = TRUE))
  )
  # expect error
  expect_error(
    coev_make_stancode(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = tree
    )
  )
})

test_that("coev_make_stancode() throws error if tree is not of class phylo", {
  # simulate data
  set.seed(1)
  n <- 20
  tree <- ape::rtree(n)
  d <- data.frame(
    id = tree$tip.label,
    x = rbinom(n, size = 1, prob = 0.5),
    y = ordered(sample(1:4, size = n, replace = TRUE))
  )
  # expect error
  expect_error(
    coev_make_stancode(
      data = d,
      variables = list(
        x = "bernoulli_logit",
        y = "ordered_logistic"
      ),
      id = "id",
      tree = "testing" # not of class phylo
    )
  )
})

test_that("coev_make_stancode() returns a character of length one", {
  # simulate data
  set.seed(1)
  n <- 20
  tree <- ape::rtree(n)
  d <- data.frame(
    id = tree$tip.label,
    x = rbinom(n, size = 1, prob = 0.5),
    y = ordered(sample(1:4, size = n, replace = TRUE))
  )
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

test_that("coev_make_stancode() creates Stan code that is syntactically correct", {
  # simulate data
  set.seed(1)
  n <- 20
  tree <- ape::rtree(n)
  d <- data.frame(
    id = tree$tip.label,
    x = rbinom(n, size = 1, prob = 0.5),
    y = ordered(sample(1:4, size = n, replace = TRUE))
  )
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
  # check stan code is syntactically correct
  expect_no_error(sc)
  expect_true(
    cmdstanr::cmdstan_model(
      stan_file = cmdstanr::write_stan_file(sc),
      compile = FALSE
      )$check_syntax(quiet = TRUE)
    )
})
