#' Fit Bayesian dynamic coevolutionary model in Stan
#'
#' @param data An object of class \code{data.frame} (or one that can be coerced
#'   to that class) containing data of all variables used in the model.
#' @param variables A named list identifying variables that should coevolve in
#'   the model and their associated response distributions (e.g.
#'   \code{list(var1 = bernoulli_logit, var2 = ordered_logistic)}). Must identify
#'   at least two variables. Variable names must refer to valid column names in data.
#'   Currently, the only supported response distributions are \code{bernoulli_logit}
#'   and \code{ordered_logistic}.
#' @param id A character of length one identifying the variable in the data that links rows to tips
#'   on the phylogeny. Must refer to a valid column name in the data. The id column
#'   must exactly match the tip labels in the phylogeny.
#' @param tree A phylogenetic tree object of class \code{phylo}.
#' @param prior A list of priors for the model.
#' @param ... Additional arguments for \pkg{cmdstanr::sampling()}.
#'
#' @return Fitted model of class \code{CmdStanModel}.
#' @export
#'
#' @examples
#' \dontrun{
#' # simulate data
#' set.seed(1)
#' n <- 20
#' tree <- ape::rtree(n)
#' d <- data.frame(
#'   id = tree$tip.label,
#'   x = rbinom(n, size = 1, prob = 0.5),
#'   y = ordered(sample(1:4, size = n, replace = TRUE))
#' )
#' # fit dynamic coevolutionary model
#' coev_fit(
#'   data = d,
#'   variables = list(
#'     x = "bernoulli_logit",
#'     y = "ordered_logistic"
#'   ),
#'   id = "id",
#'   tree = tree,
#'   # additional arguments for cmdstanr::sample()
#'   chains = 4,
#'   parallel_chains = 4,
#'   iter_warmup = 500
#' )
#' }
coev_fit <- function(data, variables, id, tree, prior = NULL, ...) {
  # check arguments
  run_checks(data, variables, id, tree)
  # write stan code for model
  sc <- coev_make_stancode(data, variables, id, tree, prior)
  # get data list for stan
  sd <- coev_make_standata(data, variables, id, tree, prior)
  # fit model
  model <-
    cmdstanr::cmdstan_model(
      stan_file = cmdstanr::write_stan_file(sc),
      compile = TRUE
    )$sample(
      data = sd,
      ...
    )
  # return object of class 'coevfit'
  out <-
    list(
      fit = model,
      data = data,
      data_name = deparse(substitute(data)),
      variables = variables,
      id = id,
      tree = tree
    )
  class(out) <- "coevfit"
  return(out)
}
