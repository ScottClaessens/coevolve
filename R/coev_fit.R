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
#' @param tree A phylogenetic tree object of class \code{phylo}.
#' @param prior A list of priors for the model.
#' @param ... Additional arguments for \pkg{Stan}.
#'
#' @return The model.
#' @export
#'
#' @examples #coev_fit()
coev_fit <- function(data, variables, tree, prior = NULL, ...) {
  ## 0. Check arguments
  ## 1. Write Stan code for model
  ## 2. Get data list for Stan
  ## 3. Fit model
  ## 4. Return fitted model
}
