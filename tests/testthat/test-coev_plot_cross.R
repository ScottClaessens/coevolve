test_that("coev_plot_cross() produces expected errors", {
  # expect the following errors
  expect_error(
    coev_plot_cross(object = "fail"),
    "Argument 'object' must be a fitted coevolutionary model of class coevfit."
  )
})

test_that("coev_plot_cross() produces ggplot object", {
  # simulate data
  withr::with_seed(1, {
    n <- 5
    tree <- ape::rtree(n)
    d <- data.frame(
      id = tree$tip.label,
      x = rbinom(n, size = 1, prob = 0.5),
      y = rbinom(n, size = 1, prob = 0.5)
    )
  })
  m <- coev_fit(
    data = d,
    variables = list(
      x = "bernoulli_logit",
      y = "bernoulli_logit"
    ),
    id = "id",
    tree = tree,
    parallel_chains = 4,
    iter_warmup = 100,
    iter_sampling = 100,
    seed = 1
  )
  # should run without error and produce ggplot object
  expect_no_error(coev_plot_cross(m))
  expect_true(methods::is(coev_plot_cross(m), "ggplot"))
})
