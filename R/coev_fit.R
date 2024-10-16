#' Fit Bayesian dynamic coevolutionary model in Stan
#'
#' Fit Bayesian dynamic coevolutionary models in \pkg{Stan} for full Bayesian
#' inference. The model allows users to assess causal directionality
#' (X -> Y vs. Y -> X) and contingencies (X, then Y) in evolution. Several data
#' types are supported, including binary, ordinal, count, continuous, and
#' positive real variables. The model can additionally account for missing data,
#' repeated observations, and controls for spatial proximity. Model fit can be
#' assessed and compared with posterior predictive checks and cross-validation.
#'
#' The \code{coev_fit} function generates the \pkg{Stan} code for the model,
#' generates the data list, and then compiles and fits the model using the
#' \pkg{cmdstanr} package.
#'
#' @param data An object of class \code{data.frame} (or one that can be coerced
#'   to that class) containing data of all variables used in the model.
#' @param variables A named list identifying variables that should coevolve in
#'   the model and their associated response distributions as character strings
#'   (e.g. \code{list(var1 = "bernoulli_logit", var2 = "ordered_logistic")}).
#'   Must identify at least two variables. Variable names must refer to valid
#'   column names in data. Currently, the only supported response distributions
#'   are \code{bernoulli_logit}, \code{ordered_logistic},
#'   \code{poisson_softplus}, \code{negative_binomial_softplus}, \code{normal},
#'   and \code{gamma_log}.
#' @param id A character of length one identifying the variable in the data that
#'   links rows to tips on the phylogeny. Must refer to a valid column name in
#'   the data. The id column must exactly match the tip labels in the phylogeny.
#' @param tree A phylogenetic tree object of class \code{phylo} or
#'   \code{multiPhylo}. The tree(s) must be rooted and must include positive
#'   non-zero branch lengths. All trees in \code{multiPhylo} objects must have
#'   the same number of internal nodes and branches.
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
#' @param dist_cov A string specifying the covariance kernel used for Gaussian
#'   Processes over locations. Currently supported are \code{"exp_quad"}
#'   (exponentiated-quadratic kernel; default), \code{"exponential"}
#'   (exponential kernel), and \code{"matern32"} (Matern 3/2 kernel).
#' @param prior (optional) A named list of priors for the model. If not
#'   specified, the model uses default priors (see \code{help(coev_fit)}).
#'   Alternatively, the user can specify a named list of priors. The list must
#'   contain non-duplicated entries for any of the following parameters: the
#'   autoregressive effects (\code{A_diag}), the cross effects
#'   (\code{A_offdiag}), the Cholesky factor for the drift matrix (\code{L_R}),
#'   the drift std. dev. parameters (\code{Q_sigma}), the continuous time
#'   intercepts (\code{b}), the ancestral states for the traits
#'   (\code{eta_anc}), the cutpoints for ordinal variables (\code{c}), the
#'   overdispersion parameters for negative binomial variables (\code{phi}),
#'   the shape parameters for gamma variables (\code{shape}), the sigma
#'   parameters for Gaussian Processes over locations (\code{sigma_dist}), the
#'   rho parameters for Gaussian Processes over locations (\code{rho_dist}), the
#'   standard deviation parameters for non-phylogenetic group-level varying
#'   effects (\code{sigma_group}), and the Cholesky factor for the
#'   non-phylogenetic group-level correlation matrix (\code{L_group}). These
#'   must be entered with valid prior strings, e.g.
#'   \code{list(A_offdiag = "normal(0, 2)")}. Invalid prior strings will throw
#'   an error when the function internally checks the syntax of resulting Stan
#'   code.
#' @param scale Logical. If \code{TRUE} (default), variables following the
#'   \code{normal} and \code{gamma_log} response distributions are scaled before
#'   fitting the model. Continuous variables following the \code{normal}
#'   distribution are standardised (e.g., mean centered and divided by their
#'   standard deviation) and positive real variables following the
#'   \code{gamma_log} distribution are divided by the mean value without
#'   centering. This approach is recommended when using default priors to
#'   improve efficiency and ensure accurate inferences. If \code{FALSE},
#'   variables are left unscaled for model fitting. In this case, users should
#'   take care to set sensible priors on variables.
#' @param estimate_Q_offdiag Logical. If \code{TRUE} (default), the model
#'   estimates the off-diagonals for the \eqn{Q} drift matrix (i.e., correlated
#'   drift). If \code{FALSE}, the off-diagonals for the \eqn{Q} drift matrix
#'   are set to zero.
#' @param prior_only Logical. If \code{FALSE} (default), the model is fitted to
#'   the data and returns a posterior distribution. If \code{TRUE}, the model
#'   samples from the prior only, ignoring the likelihood.
#' @param ... Additional arguments for \pkg{cmdstanr::sampling()}.
#'
#' @return Fitted model of class \code{coevfit}.
#'
#' @author Scott Claessens \email{scott.claessens@@gmail.com}, Erik Ringen
#'   \email{erikjacob.ringen@@uzh.ch}
#'
#' @details Fit a Bayesian dynamic coevolutionary model in Stan. A general
#'   overview is provided in the vignette \code{vignette("coevolve")}
#'
#'   \bold{Available response distributions}
#'
#'   Variables that should coevolve in the model are declared in the
#'   \code{variables} argument, using a named list with associated response
#'   distributions. For example: \code{list(x = "bernoulli_logit", y =
#'   "ordered_logistic")}. Currently, the only supported response distributions
#'   are \code{bernoulli_logit}, \code{ordered_logistic},
#'   \code{poisson_softplus}, \code{negative_binomial_softplus}, \code{normal},
#'   and \code{gamma_log}.
#'
#'   \bold{Default prior distributions}
#'
#'   If priors are not explicitly declared by the user, default priors are used.
#'   Default priors were chosen to be weakly regularising, which improves model
#'   fitting and conservatism in parameter estimates. Default priors for Stan
#'   parameters are as follows:
#'
#'   - \code{A_diag} (autoregressive effects) = \code{std_normal()}
#'   - \code{A_offdiag} (cross effects) = \code{std_normal()}
#'   - \code{L_R} (Cholesky factor for drift matrix) =
#'   \code{lkj_corr_cholesky(4)}
#'   - \code{Q_sigma} (drift std. dev. parameters) = \code{std_normal()}
#'   - \code{b} (continuous time intercepts) = \code{std_normal()}
#'   - \code{eta_anc} (trait ancestral states) = \code{std_normal()}
#'   - \code{c} (ordinal cutpoints) = \code{normal(0, 2)}
#'   - \code{shape} (shape parameters for gamma distributions) =
#'   \code{gamma(0.01, 0.01)}
#'   - \code{sigma_dist} (sigma for Gaussian process over locations) =
#'   \code{exponential(1)}
#'   - \code{rho_dist} (rho for Gaussian process over locations) =
#'   \code{exponential(5)}
#'   - \code{sigma_group} (standard deviation for group-level varying effects) =
#'   \code{exponential(1)}
#'   - \code{L_group} (Cholesky factor for group-level varying effects) =
#'   \code{lkj_corr_cholesky(2)}
#'
#'   The default prior for \code{phi} (the overdispersion parameter for the
#'   negative-binomial distribution) is scaled automatically based on the
#'   variance of the data.
#'
#'   We recommend that users assess the suitability of these default priors by
#'   fitting the model with \code{prior_only = TRUE} and then plotting a prior
#'   predictive check for all variables using the
#'   \code{coev_plot_predictive_check()} function.
#'
#'   \bold{Handling missing data}
#'
#'   In order to retain the most information, the \code{coev_fit()} function
#'   removes cases only when data are missing for all coevolving variables. For
#'   cases where data are missing for some coevolving variables but not others,
#'   the model retains these cases and averages over the missingness in the
#'   model.
#'
#'   \bold{Dealing with repeated observations}
#'
#'   If taxa appear in the dataset multiple times (i.e., there are repeated
#'   observations), the model will automatically include non-phylogenetic
#'   group-level varying effects and correlations to account for this
#'   clustering. These parameters represent the between-taxa variation and
#'   correlation that remains after accounting for the coevolutionary process.
#'
#'   \bold{Controlling for spatial location}
#'
#'   If users declare a distance matrix with the \code{dist_mat} argument,
#'   the model will include several Gaussian processes (one per coevolving
#'   variable) over spatial locations.
#'
#' @references
#' Ringen, E., Martin, J. S., & Jaeggi, A. (2021). Novel phylogenetic methods
#' reveal that resource-use intensification drives the evolution of "complex"
#' societies. \emph{EcoEvoRXiv}. \code{doi:10.32942/osf.io/wfp95}
#'
#' Sheehan, O., Watts, J., Gray, R. D., Bulbulia, J., Claessens, S., Ringen,
#' E. J., & Atkinson, Q. D. (2023). Coevolution of religious and political
#' authority in Austronesian societies. \emph{Nature Human Behaviour},
#' \emph{7}(1), 38-45. \code{10.1038/s41562-022-01471-y}
#'
#' @seealso \code{\link{coev_make_stancode}}, \code{\link{coev_make_standata}}
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
#'
#' @export
coev_fit <- function(data, variables, id, tree,
                     effects_mat = NULL, dist_mat = NULL,
                     dist_cov = "exp_quad",
                     prior = NULL, scale = TRUE,
                     estimate_Q_offdiag = TRUE,
                     prior_only = FALSE, ...) {
  # check arguments
  run_checks(data, variables, id, tree, effects_mat, dist_mat,
             dist_cov, prior, scale, estimate_Q_offdiag, prior_only)
  # write stan code for model
  sc <- coev_make_stancode(data, variables, id, tree, effects_mat,
                           dist_mat, dist_cov, prior, scale,
                           estimate_Q_offdiag, prior_only)
  # get data list for stan
  sd <- coev_make_standata(data, variables, id, tree, effects_mat,
                           dist_mat, dist_cov, prior, scale,
                           estimate_Q_offdiag, prior_only)
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
      tree_name = deparse(substitute(tree)),
      stan_code = sc,
      stan_data = sd,
      effects_mat = sd$effects_mat,
      dist_mat = sd$dist_mat,
      dist_cov = dist_cov,
      scale = scale,
      estimate_Q_offdiag = estimate_Q_offdiag,
      prior_only = prior_only
    )
  class(out) <- "coevfit"
  return(out)
}
