test_that("coev_get_delta_theta() produces expected errors and output", {
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
  # fit model
  m <- coev_fit(
    data = d,
    variables = list(
      x = "bernoulli_logit",
      y = "bernoulli_logit"
    ),
    id = "id",
    tree = tree,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    seed = 1
  )
  # expect the following errors
  expect_error(
    coev_get_delta_theta(object = "fail"),
    "Argument 'object' must be a fitted coevolutionary model of class coevfit."
  )
  expect_error(
    coev_get_delta_theta(object = m, response = c(NA, "fail")),
    "Argument 'response' must be a character string of length one."
  )
  expect_error(
    coev_get_delta_theta(object = m, response = "fail"),
    "Argument 'response' must be a variable included in the fitted model."
  )
  expect_error(
    coev_get_delta_theta(object = m, response = "x", predictor = c(NA, "fail")),
    "Argument 'predictor' must be a character string of length one."
  )
  expect_error(
    coev_get_delta_theta(object = m, response = "x", predictor = "fail"),
    "Argument 'predictor' must be a variable included in the fitted model."
  )
  expect_error(
    coev_get_delta_theta(object = m, response = "x", predictor = "x"),
    "Argument 'response' and 'predictor' must refer to different variables."
  )
  # should run without error and produce draws_array object
  expect_no_error(
    coev_get_delta_theta(m, response = "x", predictor = "y")
  )
  expect_true(
    methods::is(
      coev_get_delta_theta(m, response = "x", predictor = "y"),
      "draws_array"
    )
  )
})
