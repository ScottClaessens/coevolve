test_that("coev_plot_delta_theta() produces expected errors and output", {
  # simulate data
  withr::with_seed(1, {
    n <- 5
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
    iter_warmup = 500,
    iter_sampling = 500,
    seed = 1
  )
  # expect the following errors
  expect_error(
    coev_plot_delta_theta(object = "fail"),
    "Argument 'object' must be a fitted coevolutionary model of class coevfit."
  )
  expect_error(
    coev_plot_delta_theta(object = m, variables = NA),
    "Argument 'variables' must be a character vector."
  )
  expect_error(
    coev_plot_delta_theta(object = m, variables = "fail"),
    "Argument 'variables' must be of length > 1."
  )
  expect_error(
    coev_plot_delta_theta(object = m, variables = c("x", "y", "fail")),
    "Some variables in 'variables' are not included in the fitted model."
  )
  expect_error(
    coev_plot_delta_theta(object = m, variables = c("x", "y", "y")),
    "Argument 'variables' contains duplicates."
  )
  # should run without error and produce ggplot object
  expect_no_error(coev_plot_delta_theta(m))
  expect_true(methods::is(coev_plot_delta_theta(m), "ggplot"))
})
