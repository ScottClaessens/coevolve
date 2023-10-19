test_that("coev_make_stancode() creates Stan code that is syntactically correct", {
  # simulate data
  set.seed(1)
  n <- 100
  d <- data.frame(
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
      )
    )
  # check stan syntax is correct
  expect_true(
    cmdstanr::cmdstan_model(
      stan_file = cmdstanr::write_stan_file(sc),
      compile = FALSE
      )$check_syntax(quiet = TRUE)
    )
})
