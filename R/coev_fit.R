#' Fit Bayesian dynamic coevolutionary model in Stan
#'
#' @param data An object of class \code{data.frame} (or one that can be coerced
#'   to that class) containing data of all variables used in the model.
#' @param variables A named list identifying variables that should coevolve in
#'   the model and their associated response distributions as character strings
#'   (e.g. \code{list(var1 = "bernoulli_logit", var2 = "ordered_logistic")}).
#'   Must identify at least two variables. Variable names must refer to valid
#'   column names in data. Currently, the only supported response distributions
#'   are \code{bernoulli_logit}, \code{ordered_logistic},
#'   \code{poisson_softplus}, \code{normal}, \code{student_t}, \code{lognormal},
#'   and \code{negative_binomial_softplus}.
#' @param id A character of length one identifying the variable in the data that
#'   links rows to tips on the phylogeny. Must refer to a valid column name in
#'   the data. The id column must exactly match the tip labels in the phylogeny.
#' @param tree A phylogenetic tree object of class \code{phylo}.
#' @param effects_mat (optional) A boolean matrix with row and column names
#'   exactly matching the variables declared for the model. If not specified,
#'   all cross-lagged effects will be estimated in the model. If specified, the
#'   model will only estimate cross-lagged effects where cells in the matrix are
#'   TRUE and will ignore cross-lagged effects where cells in the matrix are
#'   FALSE. In the matrix, columns represent predictor variables and rows
#'   represent outcome variables. All autoregressive effects (e.g., X -> X) must
#'   be TRUE in the matrix.
#' @param dist_mat (optional) A distance matrix with row and column names
#'   exactly matching the tip labels in the phylogeny. If specified, the model
#'   will additionally control for spatial location by including a separate
#'   Gaussian Process over locations for every coevolving variable in the model.
#' @param prior (optional) A named list of priors for the model. If not
#'   specified, the model uses default priors (see Stan code). Alternatively,
#'   the user can specify a named list of priors. The list must contain
#'   non-duplicated entries for any of the following variables: the
#'   autoregressive effects (\code{A_diag}), the cross effects
#'   (\code{A_offdiag}), the drift scale parameters (\code{Q_diag}), the
#'   continuous time intercepts (\code{b}), the ancestral states for the traits
#'   (\code{eta_anc}), the cutpoints for ordinal variables (\code{c}), the
#'   overdispersion parameters for negative binomial variables (\code{phi}),
#'   the degrees of freedom parameters for Student t variables (\code{nu}),
#'   the sigma parameters for Gaussian Processes over locations
#'   (\code{sigma_dist}), the rho parameters for Gaussian Processes over
#'   locations (\code{rho_dist}), the standard deviation parameters for
#'   non-phylogenetic group-level varying effects (\code{sigma_group}), and the
#'   Cholesky factor for the non-phylogenetic group-level correlation matrix
#'   (\code{L_group}). These must be entered with valid prior strings, e.g.
#'   \code{list(A_offdiag = "normal(0, 2)")}. Invalid prior strings will throw
#'   an error when the function internally checks the syntax of resulting Stan
#'   code.
#' @param scale Logical. If \code{TRUE} (default), continuous and positive real
#'   variables following the \code{normal}, \code{student_t}, and
#'   \code{lognormal} response distributions are standardised before fitting the
#'   model. This approach is recommended when using default priors to improve
#'   efficiency and ensure accurate inferences. If \code{FALSE}, variables are
#'   left unstandardised for model fitting. In this case, users should take care
#'   to set sensible priors on variables.
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
#' # fit dynamic coevolutionary model
#' fit <- coev_fit(
#'   data = authority$data,
#'   variables = list(
#'     political_authority = "ordered_logistic",
#'     religious_authority = "ordered_logistic"
#'   ),
#'   id = "language",
#'   tree = authority$phylogeny,
#'   # additional arguments for cmdstanr::sample()
#'   chains = 4,
#'   parallel_chains = 4,
#'   seed = 1
#'   )
#'
#' # print model summary
#' summary(fit)
#' }
coev_fit <- function(data, variables, id, tree,
                     effects_mat = NULL, dist_mat = NULL,
                     prior = NULL, scale = TRUE, prior_only = FALSE, ...) {
  # check arguments
  run_checks(data, variables, id, tree, effects_mat,
             dist_mat, prior, scale, prior_only)
  # write stan code for model
  sc <- coev_make_stancode(data, variables, id, tree, effects_mat,
                           dist_mat, prior, scale, prior_only)
  # get data list for stan
  sd <- coev_make_standata(data, variables, id, tree, effects_mat,
                           dist_mat, prior, scale, prior_only)
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
      stan_code = sc,
      stan_data = sd,
      effects_mat = sd$effects_mat,
      dist_mat = sd$dist_mat,
      scale = scale,
      prior_only = prior_only
    )
  class(out) <- "coevfit"
  return(out)
}
