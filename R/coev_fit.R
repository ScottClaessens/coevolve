#' Fit Bayesian dynamic coevolutionary model in Stan
#'
#' @param data An object of class \code{data.frame} (or one that can be coerced
#'   to that class) containing data of all variables used in the model.
#' @param variables A named list identifying variables that should coevolve in
#'   the model and their associated response distributions as character strings (e.g.
#'   \code{list(var1 = "bernoulli_logit", var2 = "ordered_logistic")}). Must identify
#'   at least two variables. Variable names must refer to valid column names in data.
#'   Currently, the only supported response distributions are \code{bernoulli_logit},
#'   \code{ordered_logistic}, \code{poisson_softmax}, \code{normal}, and \code{lognormal}.
#' @param id A character of length one identifying the variable in the data that
#'   links rows to tips on the phylogeny. Must refer to a valid column name in
#'   the data. The id column must exactly match the tip labels in the phylogeny.
#' @param tree A phylogenetic tree object of class \code{phylo}.
#' @param dist_mat (optional) A distance matrix with row and column names exactly
#'   matching the tip labels in the phylogeny. If specified, the model will
#'   additionally control for spatial location by including a separate Gaussian
#'   Process over locations for every coevolving variable in the model.
#' @param prior (optional) A named list of priors for the model. If not specified,
#'   the model uses default priors (see Stan code). Alternatively, the user can
#'   specify a named list of priors. The list must contain non-duplicated entries
#'   for any of the following variables: the autoregressive and cross-effects
#'   (\code{alpha}), the drift scale parameters (\code{sigma}), the continuous
#'   time intercepts (\code{b}), the ancestral states for the traits (\code{eta_anc}),
#'   the cutpoints for ordinal variables (\code{c}), the sigma parameter(s) for
#'   Gaussian Processes over locations (\code{sigma_dist}), and the rho parameter(s)
#'   for Gaussian Processes over locations (\code{rho_dist}). These must be
#'   entered with valid prior strings, e.g. \code{list(alpha = "normal(0, 2)")}.
#'   Invalid prior strings will throw an error when the function internally checks
#'   the syntax of resulting Stan code.
#' @param prior_only Logical. If \code{FALSE} (default), the model is fitted to
#'   the data and returns a posterior distribution. If \code{TRUE}, the model
#'   samples from the prior only, ignoring the likelihood.
#' @param ... Additional arguments for \pkg{cmdstanr::sampling()}.
#'
#' @return Fitted model of class \code{coevfit}.
#' @export
#'
#' @examples
#' \dontrun{
#' # simulate data
#' n <- 20
#' tree <- ape::rcoal(n)
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
#'   iter_warmup = 500,
#'   seed = 1
#' )
#' }
coev_fit <- function(data, variables, id, tree,
                     dist_mat = NULL, prior = NULL,
                     prior_only = FALSE, ...) {
  # check arguments
  run_checks(data, variables, id, tree, dist_mat, prior, prior_only)
  # write stan code for model
  sc <- coev_make_stancode(data, variables, id, tree, dist_mat, prior, prior_only)
  # get data list for stan
  sd <- coev_make_standata(data, variables, id, tree, dist_mat, prior, prior_only)
  # fit model
  model <-
    cmdstanr::cmdstan_model(
      stan_file = cmdstanr::write_stan_file(sc),
      compile = TRUE
    )$sample(
      data = sd,
      show_exceptions = FALSE,
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
      tree = tree,
      dist_mat = sd$dist_mat,
      prior_only = prior_only
    )
  class(out) <- "coevfit"
  return(out)
}
