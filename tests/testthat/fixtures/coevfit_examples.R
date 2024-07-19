set.seed(1234)
library(coevolve)

# simulate data
n <- 5
tree <- ape::rcoal(n)
d <- data.frame(
  id = tree$tip.label,
  v = rnorm(n),
  w = rbinom(n, size = 1, prob = 0.5),
  x = ordered(sample(1:4, size = n, replace = TRUE)),
  y = as.integer(rnbinom(n, mu = 3, size = 1)),
  z = as.integer(rnbinom(n, mu = 3, size = 1))
)

# stan arguments
warmup <- 50
iter <- 50
chains <- 1

# fit model without distance matrix
coevfit_example1 <-
  coev_fit(
    data = d,
    variables = list(
      v = "normal",
      w = "bernoulli_logit",
      x = "ordered_logistic",
      y = "poisson_softplus"
    ),
    id = "id",
    tree = tree,
    chains = chains,
    iter_warmup = warmup,
    iter_sampling = iter,
    seed = 1
    )

# fit model with distance matrix
dist_mat <- as.matrix(dist(rnorm(n)))
rownames(dist_mat) <- colnames(dist_mat) <- d$id
coevfit_example2 <-
  coev_fit(
    data = d,
    variables = list(
      w = "bernoulli_logit",
      x = "ordered_logistic"
    ),
    id = "id",
    tree = tree,
    dist_mat = dist_mat,
    chains = chains,
    iter_warmup = warmup,
    iter_sampling = iter,
    seed = 1
    )

# fit prior only model
coevfit_example3 <-
  coev_fit(
    data = d,
    variables = list(
      w = "bernoulli_logit",
      x = "ordered_logistic"
    ),
    id = "id",
    tree = tree,
    chains = chains,
    iter_warmup = warmup,
    iter_sampling = iter,
    seed = 1,
    prior_only = TRUE
    )

# fit negative binomial model
coevfit_example4 <-
  coev_fit(
    data = d,
    variables = list(
      y = "negative_binomial_softplus",
      z = "poisson_softplus"
    ),
    id = "id",
    tree = tree,
    chains = chains,
    iter_warmup = warmup,
    iter_sampling = iter,
    seed = 1
    )

# fit model with effects matrix
effects_mat <- matrix(
  c(TRUE, TRUE,
    FALSE,TRUE),
  byrow = TRUE,
  nrow = 2,
  ncol = 2,
  dimnames = list(c("w","x"),c("w","x"))
)
coevfit_example5 <-
  coev_fit(
    data = d,
    variables = list(
      w = "bernoulli_logit",
      x = "ordered_logistic"
    ),
    id = "id",
    tree = tree,
    effects_mat = effects_mat,
    chains = chains,
    iter_warmup = warmup,
    iter_sampling = iter,
    seed = 1
    )

# fit model with missing data
d$w[c(1,2)] <- NA
d$x[c(1,3)] <- NA
coevfit_example6 <-
  coev_fit(
    data = d,
    variables = list(
      w = "bernoulli_logit",
      x = "ordered_logistic"
    ),
    id = "id",
    tree = tree,
    chains = chains,
    iter_warmup = warmup,
    iter_sampling = iter,
    seed = 1
    )

# update cmdstanr file locations
update_file_location <- function(coevfit) {
  coevfit$fit$save_output_files(
    dir = "./tests/testthat/fixtures",
    basename = deparse(substitute(coevfit)),
    timestamp = FALSE,
    random = FALSE
  )
}
suppressMessages({
  update_file_location(coevfit_example1)
  update_file_location(coevfit_example2)
  update_file_location(coevfit_example3)
  update_file_location(coevfit_example4)
  update_file_location(coevfit_example5)
  update_file_location(coevfit_example6)
})

# save coevfit objects as rds files
save_coevfit_rds <- function(coevfit) {
  saveRDS(
    coevfit,
    file = paste0(
      "./tests/testthat/fixtures/",
      deparse(substitute(coevfit)),
      ".rds"
      )
    )
}
save_coevfit_rds(coevfit_example1)
save_coevfit_rds(coevfit_example2)
save_coevfit_rds(coevfit_example3)
save_coevfit_rds(coevfit_example4)
save_coevfit_rds(coevfit_example5)
save_coevfit_rds(coevfit_example6)
