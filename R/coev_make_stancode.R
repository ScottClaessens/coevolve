#' Make Stan code for dynamic coevolutionary model
#'
#' Make the \pkg{Stan} code for the Bayesian dynamic coevolutionary model.
#' \pkg{Stan} code is generated, checked for syntactical errors, and then
#' returned as a character string.
#'
#' @srrstats {G1.3, G1.4, G2.1a} Function documentation begins here, with
#'   expected data types and definitions of statistical terminology and inputs
#' @srrstats {G2.0a} Secondary documentation on expected argument length (see
#'   "id" and "dist_cov")
#' @srrstats {G2.3, G2.3b} Documenting that character parameters are
#'   strictly case-sensitive (see "id" and "dist_cov")
#' @srrstats {G2.5} Secondary documentation of ordered factors (see "variables")
#' @srrstats {G2.14, BS3.0} Option for missing data handling (see
#'   "complete_cases") and documentation of missing data handling
#' @srrstats {G3.1, G3.1a} Users can choose the covariance function underlying
#'   the spatial Gaussian Process (see "dist_cov")
#' @srrstats {BS1.1} Data entry is described in the secondary documentation for
#'   the "data" parameter and in code examples
#' @srrstats {BS1.2, BS1.2c} Specification of prior distributions is described
#'   in secondary documentation for the "prior" parameter
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
#'   and \code{gamma_log}. Bernoulli variables must be 0/1 integers, ordered
#'   variables must be ordered factors, Poisson and negative binomial variables
#'   must be positive integers, normal variables must be continuous numeric,
#'   and gamma variables must be positive numeric.
#' @param id A character of length one identifying the variable in the data that
#'   links rows to tips on the phylogeny (strictly case-sensitive). Must refer
#'   to a valid column name in the data. The id column must exactly match the
#'   tip labels in the phylogeny.
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
#' @param complete_cases (optional) Logical. If \code{FALSE} (default), all
#'   missing values are imputed by the model. If \code{TRUE}, taxa with missing
#'   data are excluded.
#' @param dist_mat (optional) A distance matrix with row and column names
#'   exactly matching the tip labels in the phylogeny. If specified, the model
#'   will additionally control for spatial location by including a separate
#'   Gaussian Process over locations for every coevolving variable in the model.
#' @param dist_cov A string of length one specifying the covariance kernel used
#'   for Gaussian Processes over locations (strictly case-sensitive). Currently
#'   supported are \code{"exp_quad"} (exponentiated-quadratic kernel; default),
#'   \code{"exponential"} (exponential kernel), and \code{"matern32"}
#'   (Matern 3/2 kernel).
#' @param measurement_error (optional) A named list of coevolving variables and
#'   their associated columns in the dataset containing standard errors. Only
#'   valid for normally-distributed variables. For example, if we declare
#'   \code{variables = list(x = "normal", y = "normal")}, then we could set
#'   \code{measurement_error = list(x = "x_std_err")} to tell the function to
#'   include measurement error on \code{x} using standard errors from the
#'   \code{x_std_err} column of the dataset.
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
#'   residual standard deviations when there are repeated observations
#'   (\code{sigma_residual}), and the Cholesky factor for the residual
#'   correlations when there are repeated observations (\code{L_residual}).
#'   These must be entered with valid prior strings, e.g.
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
#' @param estimate_correlated_drift Logical. If \code{TRUE} (default), the model
#'   estimates the off-diagonals for the \deqn{Q} drift matrix (i.e., correlated
#'   drift). If \code{FALSE}, the off-diagonals for the \deqn{Q} drift matrix
#'   are set to zero.
#' @param estimate_residual Logical. If \code{TRUE} (default), the model
#'   estimates residual standard deviations and residual correlations when there
#'   are repeated observations for taxa. If \code{FALSE}, residual standard
#'   deviations and residual correlations are not estimated. The latter may be
#'   preferable in cases where repeated observations are sparse (e.g., only some
#'   taxa have only few repeated observations). This argument only applies when
#'   repeated observations are present in the data.
#' @param log_lik Logical. Set to \code{FALSE} by default. If \code{TRUE}, the
#'   model returns the pointwise log likelihood, which can be used to calculate
#'   WAIC and LOO.
#' @param prior_only Logical. If \code{FALSE} (default), the model is fitted to
#'   the data and returns a posterior distribution. If \code{TRUE}, the model
#'   samples from the prior only, ignoring the likelihood.
#'
#' @returns A character string containing the \pkg{Stan} code to fit the dynamic
#'   coevolutionary model
#'
#' @author Scott Claessens \email{scott.claessens@@gmail.com}, Erik Ringen
#'   \email{erikjacob.ringen@@uzh.ch}
#'
#' @details For further details, see \code{help(coev_fit)}
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
#' @seealso \code{\link{coev_make_standata}}, \code{\link{coev_fit}}
#'
#' @examples
#' # make stan code
#' stan_code <- coev_make_stancode(
#'   data = authority$data,
#'   variables = list(
#'     political_authority = "ordered_logistic",
#'     religious_authority = "ordered_logistic"
#'   ),
#'   id = "language",
#'   tree = authority$phylogeny
#' )
#'
#' # include effects matrix
#' effects_mat <-
#'   matrix(
#'     c(TRUE, TRUE,
#'       FALSE, TRUE),
#'     nrow = 2,
#'     dimnames = list(
#'       c("political_authority", "religious_authority"),
#'       c("political_authority", "religious_authority")
#'     )
#'   )
#' stan_code <- coev_make_stancode(
#'   data = authority$data,
#'   variables = list(
#'     political_authority = "ordered_logistic",
#'     religious_authority = "ordered_logistic"
#'   ),
#'   id = "language",
#'   tree = authority$phylogeny,
#'   effects_mat = effects_mat
#' )
#'
#' # include distance matrix
#' stan_code <- coev_make_stancode(
#'   data = authority$data,
#'   variables = list(
#'     political_authority = "ordered_logistic",
#'     religious_authority = "ordered_logistic"
#'   ),
#'   id = "language",
#'   tree = authority$phylogeny,
#'   dist_mat = authority$distance_matrix
#' )
#'
#' # include measurement error
#' d <- authority$data
#' d$x <- rnorm(nrow(d))
#' d$y <- rnorm(nrow(d))
#' d$x_std_err <- rexp(nrow(d))
#' d$y_std_err <- rexp(nrow(d))
#' stan_code <- coev_make_stancode(
#'   data = d,
#'   variables = list(
#'     x = "normal",
#'     y = "normal"
#'   ),
#'   id = "language",
#'   tree = authority$phylogeny,
#'   measurement_error = list(
#'     x = "x_std_err",
#'     y = "y_std_err"
#'   )
#' )
#'
#' # set manual priors
#' stan_code <- coev_make_stancode(
#'   data = authority$data,
#'   variables = list(
#'     political_authority = "ordered_logistic",
#'     religious_authority = "ordered_logistic"
#'   ),
#'   id = "language",
#'   tree = authority$phylogeny,
#'   prior = list(A_offdiag = "normal(0, 2)")
#' )
#'
#' @export
coev_make_stancode <- function(data, variables, id, tree,
                               effects_mat = NULL, complete_cases = FALSE,
                               dist_mat = NULL, dist_cov = "exp_quad",
                               measurement_error = NULL,
                               prior = NULL, scale = TRUE,
                               estimate_correlated_drift = TRUE,
                               estimate_residual = TRUE,
                               log_lik = FALSE,
                               prior_only = FALSE) {
  #' @srrstats {BS2.1} Pre-processing routines in this function ensure that all
  #'   input data is dimensionally commensurate
  # check arguments
  run_checks(data, variables, id, tree, effects_mat, complete_cases, dist_mat,
             dist_cov, measurement_error, prior, scale,
             estimate_correlated_drift, estimate_residual, log_lik, prior_only)
  # coerce data argument to data frame
  #' @srrstats {G2.7, G2.10} Accepts multiple tabular forms, ensures data frame
  data <- as.data.frame(data)
  # extract distributions and variable names from named list
  #' @srrstats {G2.4, G2.4c} Convert to character
  distributions <- as.character(variables)
  variables <- names(variables)
  # get default priors
  priors <-
    list(
      b              = "std_normal()",
      eta_anc        = "std_normal()",
      A_offdiag      = "std_normal()",
      A_diag         = "std_normal()",
      L_R            = "lkj_corr_cholesky(4)",
      Q_sigma        = "std_normal()",
      c              = "normal(0, 2)",
      shape          = "gamma(0.01, 0.01)",
      sigma_dist     = "exponential(1)",
      rho_dist       = "exponential(5)",
      sigma_residual = "exponential(1)",
      L_residual     = "lkj_corr_cholesky(4)"
    )
  # note: default prior for phi (overdispersion) set within the model code
  # replace priors if user has explicitly set them
  if (!is.null(prior)) {
    for (i in names(prior)) {
      priors[[i]] <- prior[[i]]
    }
  }
  # put stan code together
  sc <- paste0(
    "// Generated with coevolve ",
    utils::packageVersion("coevolve"),
    "\n\n",
    write_functions_block(),
    "\n\n",
    write_data_block(measurement_error, dist_mat),
    "\n\n",
    write_transformed_data_block(distributions, priors),
    "\n\n",
    write_parameters_block(data, variables, distributions, id, dist_mat,
                           estimate_correlated_drift, estimate_residual),
    "\n\n",
    write_transformed_pars_block(data, distributions, id, dist_mat,
                                 dist_cov, estimate_correlated_drift,
                                 estimate_residual),
    "\n\n",
    write_model_block(data, distributions, id, dist_mat, priors,
                      measurement_error, estimate_correlated_drift,
                      estimate_residual),
    "\n\n",
    write_gen_quantities_block(data, distributions, id, dist_mat,
                               measurement_error, estimate_correlated_drift,
                               estimate_residual, log_lik)
  )
  #' @srrstats {BS2.2, BS2.3, BS2.4, BS2.5} Checking distributional parameters
  #'   by confirming that the Stan code is syntactically correct
  # check that stan code is syntactically correct
  # if not (likely due to invalid prior string) return error
  cmdstanr::cmdstan_model(
    stan_file = cmdstanr::write_stan_file(sc),
    compile = FALSE
  )$check_syntax(quiet = TRUE)
  # produce warnings for gaussian processes and/or random effects
  if (!is.null(dist_mat)) {
    message(
      paste0(
        "Note: Distance matrix detected. Gaussian processes over spatial ",
        "distances have been included for each variable in the model ",
        "using the '", dist_cov, "' covariance kernel."
      )
    )
  }
  if (any(duplicated(data[, id])) && estimate_residual) {
    message(
      paste0(
        "Note: Repeated observations detected. Residual standard deviations ",
        "and correlations have been included in the model. To turn off this ",
        "behaviour, set estimate_residual = FALSE."
      )
    )
  }
  # produce warning that repeated models with mix of gaussian and non-gaussian
  # is experimental at this stage
  if (any(duplicated(data[, id])) && estimate_residual &&
        "normal" %in% distributions && !all(distributions == "normal")) {
    message(
      paste0(
        "Note: Repeated observations models with a mixture of ",
        "normally-distributed and non-normally-distributed variables are ",
        "currently experimental. Be sure to check models for convergence."
      )
    )
  }
  # return stan code
  sc
}

#' Internal function for rendering Stan code from whisker template
#'
#' @srrstats {G1.4a} Non-exported function documented here
#'
#' @description Renders Stan code from whisker template.
#'
#' @returns Character string
#'
#' @noRd
render_stan_template <- function(filepath, data = parent.frame()) {
  whisker::whisker.render(
    template = readLines(
      system.file(
        filepath,
        package = "coevolve"
      )
    ),
    data = data
  )
}

#' Internal function for writing the Stan functions block
#'
#' @srrstats {G1.4a} Non-exported function documented here
#'
#' @description Writes the Stan functions block for
#'   \code{\link{coev_make_stancode}}.
#'
#' @returns Character string
#'
#' @noRd
write_functions_block <- function() {
  render_stan_template(
    filepath = "stan/templates/01-functions.stan"
  )
}

#' Internal function for writing the Stan data block
#'
#' @srrstats {G1.4a} Non-exported function documented here
#'
#' @description Writes the Stan data block for \code{\link{coev_make_stancode}}.
#'
#' @returns Character string
#'
#' @noRd
write_data_block <- function(measurement_error, dist_mat) {
  render_stan_template(
    filepath = "stan/templates/02-data.stan",
    data = list(
      measurement_error = !is.null(measurement_error),
      dist_mat = !is.null(dist_mat)
    )
  )
}

#' Internal function for writing the Stan transformed data block
#'
#' @srrstats {G1.4a} Non-exported function documented here
#'
#' @description Writes the Stan transformed data block for
#'   \code{\link{coev_make_stancode}}.
#'
#' @returns Character string
#'
#' @noRd
write_transformed_data_block <- function(distributions, priors) {
  # sequence of variables for template
  variable_seq <- lapply(seq_along(distributions), function(j) list(j = j))
  # negative binomial variables for template
  if ("negative_binomial_softplus" %in% distributions && is.null(priors$phi)) {
    neg_binomial_seq <- lapply(
      which(distributions == "negative_binomial_softplus"),
      function(j) list(j = j)
    )
  } else {
    neg_binomial_seq <- FALSE
  }
  # render template
  render_stan_template(
    filepath = "stan/templates/03-transformed-data.stan",
    data = list(
      variable_seq = variable_seq,
      neg_binomial_seq = neg_binomial_seq
    )
  )
}

#' Internal function for writing the Stan parameters block
#'
#' @srrstats {G1.4a} Non-exported function documented here
#'
#' @description Writes the Stan parameters block for
#'   \code{\link{coev_make_stancode}}.
#'
#' @returns Character string
#'
#' @noRd
write_parameters_block <- function(data, variables, distributions, id, dist_mat,
                                   estimate_correlated_drift,
                                   estimate_residual) {
  # ordered variables for template
  if ("ordered_logistic" %in% distributions) {
    ordered_seq <- lapply(
      which(distributions == "ordered_logistic"),
      function(j) {
        list(
          j = j,
          # calculate number of cut points (number of levels - 1)
          #' @srrstats {G2.4, G2.4b} Convert to continuous calculate cutpoints
          #' @srrstats {G2.15} Software does not assume non-missingness (na.rm)
          num_cuts = max(as.numeric(data[, variables[j]]), na.rm = TRUE) - 1
        )
      }
    )
  } else {
    ordered_seq <- FALSE
  }
  # negative binomial variables for template
  if ("negative_binomial_softplus" %in% distributions) {
    neg_binomial_seq <- lapply(
      which(distributions == "negative_binomial_softplus"),
      function(j) list(j = j)
    )
  } else {
    neg_binomial_seq <- FALSE
  }
  # gamma variables for template
  if ("gamma_log" %in% distributions) {
    gamma_seq <- lapply(
      which(distributions == "gamma_log"),
      function(j) list(j = j)
    )
  } else {
    gamma_seq <- FALSE
  }
  # render template
  render_stan_template(
    filepath = "stan/templates/04-parameters.stan",
    data = list(
      estimate_correlated_drift = estimate_correlated_drift,
      ordered_seq = ordered_seq,
      neg_binomial_seq = neg_binomial_seq,
      gamma_seq = gamma_seq,
      dist_mat = !is.null(dist_mat),
      repeated_measures = any(duplicated(data[, id])) && estimate_residual
    )
  )
}

#' Internal function for writing the Stan transformed parameters block
#'
#' @srrstats {G1.4a} Non-exported function documented here
#'
#' @description Writes the Stan transformed parameters block for
#'   \code{\link{coev_make_stancode}}.
#'
#' @returns Character string
#'
#' @noRd
write_transformed_pars_block <- function(data, distributions, id, dist_mat,
                                         dist_cov, estimate_correlated_drift,
                                         estimate_residual) {
  # calculate tdrift?
  tdrift <-
    (any(duplicated(data[, id])) && estimate_residual) ||
    !("normal" %in% distributions)
  # calculate residual?
  residual <-
    any(duplicated(data[, id])) && estimate_residual &&
    !("normal" %in% distributions)
  # get gaussian process kernel code for template
  dist_cov_code <- NULL
  if (dist_cov == "exp_quad") {
    # exponentiated quadratic kernel
    dist_cov_code <- paste0(
      "sigma_dist[j] * exp(-(square(dist_mat[i,m]) / rho_dist[j]))"
    )
  } else if (dist_cov == "exponential") {
    # exponential kernel
    dist_cov_code <-
      "sigma_dist[j] * exp(-(dist_mat[i,m] / rho_dist[j]))"
  } else if (dist_cov == "matern32") {
    # matern 3/2 kernel
    dist_cov_code <- paste0(
      "sigma_dist[j] * (1 + ((sqrt(3.0) * dist_mat[i,m]) / rho_dist[j])) * ",
      "exp(-(sqrt(3.0) * dist_mat[i,m]) / rho_dist[j])"
    )
  }
  # render template
  render_stan_template(
    filepath = "stan/templates/05-transformed-parameters.stan",
    data = list(
      estimate_correlated_drift = estimate_correlated_drift,
      no_correlated_drift = !estimate_correlated_drift,
      dist_mat = !is.null(dist_mat),
      tdrift = tdrift,
      residual = residual,
      dist_cov_code = dist_cov_code
    )
  )
}

#' Internal function for writing the Stan model block
#'
#' @srrstats {G1.4a} Non-exported function documented here
#'
#' @description Writes the Stan model block for
#'   \code{\link{coev_make_stancode}}.
#'
#' @returns Character string
#'
#' @noRd
write_model_block <- function(data, distributions, id, dist_mat, priors,
                              measurement_error, estimate_correlated_drift,
                              estimate_residual) {
  # check for repeated
  repeated <- any(duplicated(data[, id])) && estimate_residual
  # add priors for terminal_drift when:
  # 1. there are no repeated measures and no gaussian variables  OR
  # 2. there are repeated measures and estimate_residual = TRUE
  add_terminal_drift_prior <-
    (!any(duplicated(data[, id])) && !("normal" %in% distributions)) ||
    (repeated)
  # which variables are ordinal?
  ordered_seq <- FALSE
  if ("ordered_logistic" %in% distributions) {
    ordered_seq <- lapply(
      which(distributions == "ordered_logistic"),
      function(j) list(j = j)
    )
  }
  # loop over phi priors for negative binomial variables
  prior_phi_default <- FALSE
  prior_phi_manual <- FALSE
  if ("negative_binomial_softplus" %in% distributions) {
    if (is.null(priors$phi)) {
      # use default prior
      prior_phi_default <- lapply(
        which(distributions == "negative_binomial_softplus"),
        function(j) list(j = j)
      )
    } else {
      # use manually specified prior
      prior_phi_manual <- lapply(
        which(distributions == "negative_binomial_softplus"),
        function(j) list(j = j, prior_phi = priors$phi)
      )
    }
  }
  # which variables are gamma?
  gamma_seq <- FALSE
  if ("gamma_log" %in% distributions) {
    gamma_seq <- lapply(
      which(distributions == "gamma_log"),
      function(j) list(j = j)
    )
  }
  # add priors for residual sds and cors
  add_priors_residual_sds_cors <- FALSE
  if (repeated) {
    add_priors_residual_sds_cors <- list(
      normal_present = FALSE,
      normal_absent = FALSE,
      prior_sigma_residual = priors$sigma_residual,
      prior_L_residual = priors$L_residual
    )
    if ("normal" %in% distributions) {
      add_priors_residual_sds_cors$normal_present <-
        list(
          is_normal = lapply(
            which(distributions == "normal"),
            function(j) list(j = j)
          ),
          is_not_normal = lapply(
            which(distributions != "normal"),
            function(j) list(j = j)
          )
        )
    } else {
      add_priors_residual_sds_cors$normal_absent <- TRUE
    }
  }
  # function to get linear model
  lmod <- function(j) {
    paste0(
      "eta[t,tip_id[i]][", j, "]",
      ifelse(
        !is.null(dist_mat),
        paste0(" + dist_v[tip_id[i],", j, "]"),
        ""
      )
    )
  }
  # loop over variables to detect normal variables
  set <- list(
    is_normal = lapply(
      which(distributions == "normal"),
      function(j) list(j = j, lmod = lmod(j))
    ),
    is_not_normal = lapply(
      which(distributions != "normal"),
      function(j) list(j = j)
    ),
    measurement_error = !is.null(measurement_error)
  )
  # set residuals when repeated observations and normal variables present
  init_residuals <- FALSE
  set_residuals <- FALSE
  if (repeated && "normal" %in% distributions) {
    init_residuals <- TRUE
    set_residuals <- set
    set_residuals$cov_matrix <- ifelse(
      !is.null(measurement_error),
      "cholesky_decompose(residual_cov)",
      "diag_pre_multiply(sigma_residual, L_residual)"
    )
  }
  # set tdrift when no repeated observations and normal variables present
  init_tdrift <- FALSE
  set_tdrift <- FALSE
  if (!(repeated) && "normal" %in% distributions) {
    init_tdrift <- TRUE
    set_tdrift <- set
    set_tdrift$cov_matrix <- ifelse(
      !is.null(measurement_error),
      "cholesky_decompose(add_diag(VCV_tips[t, tip_id[i]], se[i,]))",
      "L_VCV_tips[t, tip_id[i]]"
    )
  }
  # function to write likelihoods for non-continuous variables
  write_likelihood <- function(distribution, j) {
    # linear model
    lmod <- lmod(j)
    # append tdrift
    lmod <- ifelse(
      !(repeated) && "normal" %in% distributions,
      paste0(lmod, " + tdrift[", j, "]"),
      paste0(lmod, " + tdrift[t,tip_id[i]][", j, "]")
    )
    # append residual_v
    lmod <- ifelse(
      repeated && !("normal" %in% distributions),
      paste0(lmod, " + residual_v[i,", j, "]"),
      lmod
    )
    # append residuals
    lmod <- ifelse(
      repeated && "normal" %in% distributions,
      paste0(lmod, " + residuals[", j, "]"),
      lmod
    )
    # put together
    if (distribution == "bernoulli_logit") {
      paste0("bernoulli_logit_lpmf(to_int(y[i,", j, "]) | ", lmod, ")")
    } else if (distribution == "ordered_logistic") {
      paste0("ordered_logistic_lpmf(to_int(y[i,", j, "]) | ", lmod,
             ", c", j, ")")
    } else if (distribution == "poisson_softplus") {
      paste0("poisson_lpmf(to_int(y[i,", j, "]) | ",
             "mean(obs", j, ") * log1p_exp(", lmod, "))")
    } else if (distribution == "negative_binomial_softplus") {
      paste0("neg_binomial_2_lpmf(to_int(y[i,", j, "]) | ",
             "mean(obs", j, ") * log1p_exp(", lmod, "), phi", j, ")")
    } else if (distribution == "gamma_log") {
      paste0("gamma_lpdf(y[i,", j, "] | shape", j, ", shape", j, " / exp(",
             lmod, "))")
    }
  }
  # loop over likelihoods for non-continuous variables
  likelihoods <- FALSE
  if (!all(distributions == "normal")) {
    likelihoods <- lapply(
      which(distributions != "normal"),
      function(j) {
        list(j = j, likelihood = write_likelihood(distributions[j], j))
      }
    )
  }
  # render template
  render_stan_template(
    filepath = "stan/templates/06-model.stan",
    data = list(
      prior_b = priors$b,
      prior_eta_anc = priors$eta_anc,
      prior_A_offdiag = priors$A_offdiag,
      prior_A_diag = priors$A_diag,
      prior_L_R = priors$L_R,
      prior_Q_sigma = priors$Q_sigma,
      prior_c = priors$c,
      prior_phi_default = prior_phi_default,
      prior_phi_manual = prior_phi_manual,
      prior_shape = priors$shape,
      prior_sigma_dist = priors$sigma_dist,
      prior_rho_dist = priors$rho_dist,
      add_terminal_drift_prior = add_terminal_drift_prior,
      estimate_correlated_drift = estimate_correlated_drift,
      ordered_seq = ordered_seq,
      gamma_seq = gamma_seq,
      dist_mat = !is.null(dist_mat),
      add_priors_residual_sds_cors = add_priors_residual_sds_cors,
      init_residuals = init_residuals,
      init_tdrift = init_tdrift,
      set_residuals = set_residuals,
      set_tdrift = set_tdrift,
      likelihoods = likelihoods
    )
  )
}

#' Internal function for writing the Stan generated quantities block
#'
#' @srrstats {G1.4a} Non-exported function documented here
#'
#' @description Writes the Stan generated quantities block for
#'   \code{\link{coev_make_stancode}}.
#'
#' @returns Character string
#'
#' @noRd
write_gen_quantities_block <- function(data, distributions, id, dist_mat,
                                       measurement_error,
                                       estimate_correlated_drift,
                                       estimate_residual, log_lik) {
  # check if repeated
  repeated <- any(duplicated(data[, id])) && estimate_residual
  # terminal drift cov matrix
  terminal_drift_cov_matrix <- ifelse(
    !is.null(measurement_error) && !(repeated),
    "add_diag(VCV_tips[t, i], se[i,])",
    "VCV_tips[t, i]"
  )
  # normal present and log_lik
  normal_present_and_log_lik <- "normal" %in% distributions && log_lik
  # function to get linear model
  lmod <- function(j) {
    paste0(
      "eta[t,tip_id[i]][", j, "]",
      ifelse(
        !is.null(dist_mat),
        paste0(" + dist_v[tip_id[i],", j, "]"),
        ""
      )
    )
  }
  # loop over variables to detect normal variables
  set <- list(
    is_normal = FALSE,
    is_not_normal_and_normal_present = FALSE,
    is_not_normal_and_normal_present = FALSE
  )
  if ("normal" %in% distributions) {
    set$is_normal <- lapply(
      which(distributions == "normal"),
      function(j) list(j = j, lmod = lmod(j))
    )
    set$is_not_normal_and_normal_present <- lapply(
      which(distributions != "normal"),
      function(j) list(j = j)
    )
  } else {
    set$is_not_normal_and_normal_absent <- lapply(
      which(distributions != "normal"),
      function(j) list(j = j)
    )
  }
  # set residuals and tdrifts
  init_residuals <- FALSE
  init_tdrifts <- FALSE
  set_residuals <- FALSE
  set_tdrifts <- FALSE
  if (log_lik) {
    if (repeated) {
      init_residuals <- TRUE
      set_residuals <- set
    } else {
      init_tdrifts <- TRUE
      set_tdrifts <- set
    }
  }
  # set mu_cond and sigma_cond
  set_mu_cond_and_sigma_cond <- FALSE
  measurement_error_list <- list(
    measurement_error = !is.null(measurement_error),
    no_measurement_error = is.null(measurement_error)
  )
  if ("normal" %in% distributions && log_lik) {
    if (repeated) {
      set_mu_cond_and_sigma_cond <- list(
        repeated = measurement_error_list,
        no_repeated = FALSE
      )
    } else {
      set_mu_cond_and_sigma_cond <- list(
        repeated = FALSE,
        no_repeated = measurement_error_list
      )
    }
  }
  # function to write log likelihood statement
  write_log_lik_statement <- function(distribution, j) {
    # linear model
    lmod <- lmod(j)
    # append residuals or tdrifts
    if (repeated) {
      lmod <- paste0(lmod, " + residuals[", j, "]")
    } else {
      lmod <- paste0(lmod, " + tdrifts[", j, "]")
    }
    # put together
    if (distribution == "bernoulli_logit") {
      paste0("bernoulli_logit_lpmf(to_int(y[i,", j, "]) | ", lmod, ")")
    } else if (distribution == "ordered_logistic") {
      paste0("ordered_logistic_lpmf(to_int(y[i,", j, "]) | ",
             lmod, ", c", j, ")")
    } else if (distribution == "poisson_softplus") {
      paste0("poisson_lpmf(to_int(y[i,", j, "]) | mean(obs", j,
             ") * log1p_exp(", lmod, "))")
    } else if (distribution == "negative_binomial_softplus") {
      paste0("neg_binomial_2_lpmf(to_int(y[i,", j, "]) | mean(obs", j,
             ") * log1p_exp(", lmod, "), phi", j, ")")
    } else if (distribution == "gamma_log") {
      paste0("gamma_lpdf(y[i,", j, "] | shape", j, ", shape", j, " / exp(",
             lmod, "))")
    } else if (distribution == "normal") {
      paste0(
        "normal_lpdf(",
        ifelse(
          repeated,
          paste0("residuals[", j, "]"),
          paste0("tdrifts[", j, "]")
        ),
        " | mu_cond[", j, "], sigma_cond[", j, "])"
      )
    }
  }
  # function to write posterior prediction statement
  write_yrep_statement <- function(distribution, j) {
    # linear model
    lmod <- paste0(lmod(j), " + terminal_drift_rep[t,tip_id[i]][", j, "]")
    # append residuals_rep
    if (repeated) {
      lmod <- paste0(lmod, " + residuals_rep[", j, "]")
    }
    # put together
    if (distribution == "bernoulli_logit") {
      paste0("bernoulli_logit_rng(", lmod, ")")
    } else if (distribution == "ordered_logistic") {
      paste0("ordered_logistic_rng(", lmod, ", c", j, ")")
    } else if (distribution == "poisson_softplus") {
      paste0("poisson_rng(mean(obs", j, ") * log1p_exp(", lmod, "))")
    } else if (distribution == "negative_binomial_softplus") {
      paste0("neg_binomial_2_rng(mean(obs", j, ") * log1p_exp(",
             lmod, "), phi", j, ")")
    } else if (distribution == "gamma_log") {
      paste0("gamma_rng(shape", j, ", shape", j, " / exp(", lmod, "))")
    } else if (distribution == "normal") {
      lmod
    }
  }
  # loop over variables to get posterior predictions (and log likelihoods)
  posterior_preds_and_log_liks <-
    lapply(
      seq_along(distributions),
      function(j) {
        list(
          log_lik = if (log_lik) {
            list(
              j = j,
              log_lik_statement = write_log_lik_statement(distributions[j], j)
            )
          } else {
            FALSE
          },
          yrep = list(
            j = j,
            yrep_statement = write_yrep_statement(distributions[j], j)
          )
        )
      }
    )
  # render template
  render_stan_template(
    filepath = "stan/templates/07-generated-quantities.stan",
    data = list(
      log_lik = log_lik,
      estimate_correlated_drift = estimate_correlated_drift,
      repeated = repeated,
      terminal_drift_cov_matrix = terminal_drift_cov_matrix,
      normal_present_and_log_lik = normal_present_and_log_lik,
      init_residuals = init_residuals,
      init_tdrifts = init_tdrifts,
      set_residuals = set_residuals,
      set_tdrifts = set_tdrifts,
      set_mu_cond_and_sigma_cond = set_mu_cond_and_sigma_cond,
      posterior_preds_and_log_liks =
        posterior_preds_and_log_liks
    )
  )
}
