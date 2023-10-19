test_that("coev_make_stancode() creates Stan code that is syntactically correct", {
  # simulate data
  set.seed(1)
  n <- 100
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
  # check stan syntax is correct
  expect_true(
    cmdstanr::cmdstan_model(
      stan_file = cmdstanr::write_stan_file(sc),
      compile = FALSE
      )$check_syntax(quiet = TRUE)
    )
})
