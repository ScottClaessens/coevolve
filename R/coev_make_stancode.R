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
    "\n",
    write_functions_block(),
    "\n",
    write_data_block(measurement_error, dist_mat),
    "\n",
    write_transformed_data_block(distributions, priors),
    "\n",
    write_parameters_block(data, variables, distributions, id, dist_mat,
                           estimate_correlated_drift, estimate_residual),
    "\n",
    write_transformed_pars_block(data, distributions, id, dist_mat,
                                 dist_cov, estimate_correlated_drift,
                                 estimate_residual, measurement_error),
    "\n",
    write_model_block(data, distributions, id, dist_mat, priors,
                      measurement_error, estimate_correlated_drift,
                      estimate_residual),
    "\n",
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
  return(sc)
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
  paste0(
    "functions {\n",
    "  // Charles Driver's solver for the asymptotic Q matrix\n",
    "  matrix ksolve (matrix A, matrix Q) {\n",
    "    int d = rows(A);\n",
    "    int d2 = (d * d - d) %/% 2;\n",
    "    matrix [d + d2, d + d2] O;\n",
    "    vector [d + d2] triQ;\n",
    "    matrix[d,d] AQ;\n",
    "    int z = 0;         // z is row of output\n",
    "    for (j in 1:d) {   // for column reference of solution vector\n",
    "      for (i in 1:j) { // and row reference...\n",
    "        if (j >= i) {  // if i and j denote a covariance parameter\n",
    "          int y = 0;   // start new output row\n",
    "          z += 1;      // shift current output row down\n",
    "          for (ci in 1:d) {   // for columns and\n",
    "            for (ri in 1:d) { // rows of solution\n",
    "              if (ci >= ri) { // when in upper tri (inc diag)\n",
    "                y += 1;       // move to next column of output\n",
    "                if (i == j) { // if output row is a diag element\n",
    "                  if (ri == i) O[z, y] = 2 * A[ri, ci];\n",
    "                  if (ci == i) O[z, y] = 2 * A[ci, ri];\n",
    "                }\n",
    "                if (i != j) { // if output row is not a diag element\n",
    "                  //if column matches row, sum both A diags\n",
    "                  if (y == z) O[z, y] = A[ri, ri] + A[ci, ci];\n",
    "                  if (y != z) { // otherwise...\n",
    "                    // if solution element is related to output row...\n",
    "                    if (ci == ri) { // if solution element is variance\n",
    "                      // if variance of solution corresponds to row\n",
    "                      if (ci == i) O[z, y] = A[j, ci];\n",
    "                      // if variance of solution corresponds to col\n",
    "                      if (ci == j) O[z, y] = A[i, ci];\n",
    "                    }\n",
    "                    //if solution element is a related covariance\n",
    "                    if (ci != ri && (ri == i || ri == j ",
    "|| ci == i || ci == j )) {\n",
    "                      // for row 1,2 / 2,1 of output,\n",
    "                      // if solution row ri 1 (match)\n",
    "                      // and column ci 3, we need A[2,3]\n",
    "                      if (ri == i) O[z, y] = A[j, ci];\n",
    "                      if (ri == j) O[z, y] = A[i, ci];\n",
    "                      if (ci == i) O[z, y] = A[j, ri];\n",
    "                      if (ci == j) O[z, y] = A[i, ri];\n",
    "                    }\n",
    "                  }\n",
    "                }\n",
    "                if (is_nan(O[z, y])) O[z, y] = 0;\n",
    "              }\n",
    "            }\n",
    "          }\n",
    "        }\n",
    "      }\n",
    "    }\n",
    "    z = 0; // get upper tri of Q\n",
    "    for (j in 1:d) {\n",
    "      for (i in 1:j) {\n",
    "        z += 1;\n",
    "        triQ[z] = Q[i, j];\n",
    "      }\n",
    "    }\n",
    "    triQ = -O \\ triQ; // get upper tri of asymQ\n",
    "    z = 0; // put upper tri of asymQ into matrix\n",
    "    for (j in 1:d) {\n",
    "      for (i in 1:j) {\n",
    "        z += 1;\n",
    "        AQ[i, j] = triQ[z];\n",
    "        if (i != j) AQ[j, i] = triQ[z];\n",
    "      }\n",
    "    }\n",
    "    return AQ;\n",
    "  }\n",
    "  \n",
    "  // return number of matches of y in vector x\n",
    "  int num_matches(vector x, real y) {\n",
    "    int n = 0;\n",
    "    for (i in 1:rows(x))\n",
    "      if (x[i] == y)\n",
    "        n += 1;\n",
    "    return n;\n",
    "  }\n",
    "  \n",
    "  // return indices in vector x where x == y\n",
    "  array[] int which_equal(vector x, real y) {\n",
    "    array [num_matches(x, y)] int match_positions;\n",
    "    int pos = 1;\n",
    "    for (i in 1:rows(x)) {\n",
    "      if (x[i] == y) {\n",
    "        match_positions[pos] = i;\n",
    "        pos += 1;\n",
    "      }\n",
    "    }\n",
    "    return match_positions;\n",
    "  }\n",
    "}"
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
  sc_data <- paste0(
    "data{\n",
    "  int<lower=1> N_tips; // number of tips\n",
    "  int<lower=1> N_tree; // number of trees\n",
    "  int<lower=1> N_obs; // number of observations\n",
    "  int<lower=2> J; // number of response traits\n",
    "  int<lower=1> N_seg; // total number of segments in the trees\n",
    "  array[N_tree, N_seg] int<lower=1> node_seq; // index of tree nodes\n",
    "  array[N_tree, N_seg] int<lower=0> parent; // index of parent nodes\n",
    "  array[N_tree, N_seg] real ts; // time since parent\n",
    "  array[N_tree, N_seg] int<lower=0,upper=1> tip; // segment ends in tip\n",
    "  array[J,J] int<lower=0,upper=1> effects_mat; // effects matrix\n",
    "  int<lower=2> num_effects; // number of effects being estimated\n",
    "  matrix[N_obs,J] y; // observed data\n",
    "  matrix[N_obs,J] miss; // are data points missing?\n",
    ifelse(
      !is.null(measurement_error),
      "  matrix[N_obs,J] se; // squared standard errors\n", ""
    ),
    "  array[N_obs] int<lower=1> tip_id; // group index between 1 and N_tips\n",
    "  int<lower=1> N_unique_lengths; // number of unique branch lengths\n",
    "  array[N_unique_lengths] real unique_lengths; // unique branch lengths for caching\n",
    "  array[N_tree, N_seg] int<lower=0> length_index; // mapping from segments to unique lengths\n",
    "  array[N_tree, N_tips] int<lower=0> tip_to_seg; // mapping from tips to segments\n"
  )
  # add distance matrix if user has defined one
  if (!is.null(dist_mat)) {
    sc_data <-
      paste0(
        sc_data,
        "  matrix[N_tips,N_tips] dist_mat; // distance matrix\n"
      )
  }
  # add prior_only data variable
  paste0(
    sc_data,
    "  int<lower=0,upper=1> prior_only; // should likelihood be ignored?\n}"
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
  sc_transformed_data <- "transformed data{\n"
  for (j in seq_along(distributions)) {
    sc_transformed_data <- paste0(
      sc_transformed_data,
      "  vector[to_int(N_obs - sum(col(miss, ", j,
      ")))] obs", j, "; // observed data for variable ", j, "\n"
    )
  }
  for (j in seq_along(distributions)) {
    if (distributions[j] == "negative_binomial_softplus" &&
          is.null(priors$phi)) {
      sc_transformed_data <-
        paste0(
          sc_transformed_data,
          "  real inv_overdisp", j, "; // best guess for phi", j, "\n"
        )
    }
  }
  for (j in seq_along(distributions)) {
    sc_transformed_data <- paste0(
      sc_transformed_data,
      "  obs", j, " = col(y, ", j, ")[which_equal(col(miss, ", j, "), 0)];\n"
    )
  }
  for (j in seq_along(distributions)) {
    if (distributions[j] == "negative_binomial_softplus" &&
          is.null(priors$phi)) {
      sc_transformed_data <-
        paste0(
          sc_transformed_data,
          "  inv_overdisp", j, " = (mean(obs", j, ")^2) / ",
          "(sd(obs", j, ")^2 - mean(obs", j, "));\n"
        )
    }
  }
  paste0(sc_transformed_data, "}")
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
  sc_parameters <- paste0(
    "parameters{\n",
    "  vector<upper=0>[J] A_diag; // autoregressive terms of A\n",
    "  vector[num_effects - J] A_offdiag; // cross-lagged terms of A\n"
  )
  # add cholesky factor for Q matrix if estimating Q off diagonals
  if (estimate_correlated_drift) {
    sc_parameters <- paste0(
      sc_parameters,
      "  cholesky_factor_corr[J] L_R; // lower-tri choleksy decomp corr mat\n"
    )
  }
  sc_parameters <- paste0(
    sc_parameters,
    "  vector<lower=0>[J] Q_sigma; // std deviation parameters of the Q mat\n",
    "  vector[J] b; // SDE intercepts\n",
    "  array[N_tree] vector[J] eta_anc; // ancestral states\n",
    "  array[N_tree, N_seg - 1] vector[J] z_drift; // stochastic drift\n",
    "  array[N_tree] matrix[N_tips, J] terminal_drift; // drift for the tips\n"
  )
  for (i in seq_along(distributions)) {
    # add cut points for ordinal_logistic distributions
    if (distributions[i] == "ordered_logistic") {
      # calculate number of cut points (number of levels - 1)
      #' @srrstats {G2.4, G2.4b} Convert to continuous to calculate cutpoints
      #' @srrstats {G2.15} Software does not assume non-missingness (na.rm)
      num_cuts <- max(as.numeric(data[, variables[i]]), na.rm = TRUE) - 1
      sc_parameters <-
        paste0(
          sc_parameters,
          "  ordered[", num_cuts, "] c", i, "; ",
          "// cut points for variable ", i, "\n"
        )
    }
    # add overdispersion parameters for negative_binomial_softplus distributions
    if (distributions[i] == "negative_binomial_softplus") {
      sc_parameters <-
        paste0(
          sc_parameters,
          "  real<lower=0> phi", i, "; ",
          "// neg binom inverse overdispersion parameter for variable ", i, "\n"
        )
    }
    # add shape parameters for gamma_log distributions
    if (distributions[i] == "gamma_log") {
      sc_parameters <-
        paste0(
          sc_parameters,
          "  real<lower=0> shape", i, "; ",
          "// gamma shape parameter for variable ", i, "\n"
        )
    }
  }
  # add gaussian process parameters if distance matrix specified by user
  if (!is.null(dist_mat)) {
    sc_parameters <- paste0(
      sc_parameters,
      "  matrix[N_tips,J] dist_z; // spatial covariance random effects\n",
      "  vector<lower=0>[J] rho_dist; // covariance declining with distance\n",
      "  vector<lower=0>[J] sigma_dist; // maximum covariance\n"
    )
  }
  # add residual sds and cors if there are repeated measures
  if (any(duplicated(data[, id])) && estimate_residual) {
    sc_parameters <- paste0(
      sc_parameters,
      "  matrix[J,N_obs] residual_z;\n",
      "  vector<lower=0>[J] sigma_residual;\n",
      "  cholesky_factor_corr[J] L_residual;\n"
    )
  }
  paste0(sc_parameters, "}")
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
                                         estimate_residual, measurement_error) {
  sc_transformed_parameters <- paste0(
    "transformed parameters{\n",
    "  array[N_tree, N_seg] vector[J] eta;\n",
    "  matrix[J,J] A = diag_matrix(A_diag); // selection matrix\n",
    ifelse(
      estimate_correlated_drift,
      paste0(
        "  matrix[J,J] Q = diag_matrix(Q_sigma) * (L_R * L_R')",
        " * diag_matrix(Q_sigma); // drift matrix\n"
      ),
      "  matrix[J,J] Q = diag_matrix(Q_sigma^2); // drift matrix\n"
    ),
    "  matrix[J,J] Q_inf; // asymptotic covariance matrix\n",
    "  array[N_tree, N_seg] matrix[J,J] VCV_tips; // vcov matrix for drift\n"
  )
  # add distance random effects if distance matrix specified by user
  if (!is.null(dist_mat)) {
    sc_transformed_parameters <- paste0(
      sc_transformed_parameters,
      "  matrix[N_tips,J] dist_v; // distance covariance random effects\n"
    )
  }
  # add tdrift if repeated observations or no gaussian traits
  if ((any(duplicated(data[, id])) && estimate_residual) ||
        !("normal" %in% distributions)) {
    sc_transformed_parameters <- paste0(
      sc_transformed_parameters,
      "  array[N_tree,N_tips] vector[J] tdrift; // terminal drift\n"
    )
  }
  # add residual sds and cors if there are repeated measures and using the
  # non-centered parameterisation (there are no gaussian variables)
  if (any(duplicated(data[, id])) && estimate_residual &&
        !("normal" %in% distributions)) {
    sc_transformed_parameters <- paste0(
      sc_transformed_parameters,
      "  matrix[N_obs,J] residual_v; // residual pars\n",
      "  // scale and correlate residual pars\n",
      "  residual_v = (diag_pre_multiply(sigma_residual, L_residual)",
      " * residual_z)';\n"
    )
  }
  sc_transformed_parameters <- paste0(
    sc_transformed_parameters,
    "  // fill off diagonal of A matrix\n",
    "  {\n",
    "    int ticker = 1;\n",
    "    for (i in 1:J) {\n",
    "      for (j in 1:J) {\n",
    "        if (i != j) {\n",
    "          if (effects_mat[i,j] == 1) {\n",
    "            A[i,j] = A_offdiag[ticker];\n",
    "            ticker += 1;\n",
    "          } else if (effects_mat[i,j] == 0) {\n",
    "            A[i,j] = 0;\n",
    "          }\n",
    "        }\n",
    "      }\n",
    "    }\n",
    "  }\n",
    "  // calculate asymptotic covariance\n",
    "  Q_inf = ksolve(A, Q);\n",
    ifelse(
      estimate_residual && any(duplicated(data[, id])),
      paste0(
        "  // cache residual covariance components (computed once per iteration)\n",
        ifelse(
          !is.null(measurement_error),
          "  matrix[J,J] residual_cov_base = quad_form_diag(L_residual * L_residual', sigma_residual);\n",
          ""
        ),
        "  matrix[J,J] L_residual_scaled = diag_pre_multiply(sigma_residual, L_residual);\n"
      ),
      ""
    ),
    "  array[N_unique_lengths] matrix[J,J] L_VCV_tips_cache;\n",
    "  {\n",
    "    array[N_unique_lengths] matrix[J,J] A_delta_cache;\n",
    "    array[N_unique_lengths] matrix[J,J] VCV_cache;\n",
    "    array[N_unique_lengths] matrix[J,J] A_solve_cache;\n",
    "    array[N_unique_lengths] vector[J] A_solve_b_cache;\n",
    "    for (u in 1:N_unique_lengths) {\n",
    "      A_delta_cache[u] = matrix_exp(A * unique_lengths[u]);\n",
    "      VCV_cache[u] = Q_inf - quad_form_sym(Q_inf, A_delta_cache[u]');\n",
    "      L_VCV_tips_cache[u] = cholesky_decompose(VCV_cache[u]);\n",
      "      A_solve_cache[u] = A \\ add_diag(A_delta_cache[u], -1);\n",
    "      for (i in 1:J) {\n",
    "        for (j in 1:i) {\n",
    "          real val = 0.5 * (A_solve_cache[u][i, j] + A_solve_cache[u][j, i]);\n",
    "          A_solve_cache[u][i, j] = val;\n",
    "          A_solve_cache[u][j, i] = val;\n",
    "        }\n",
    "      }\n",
      "      A_solve_b_cache[u] = A_solve_cache[u] * b;\n",
    "    }\n",
    "    for (t in 1:N_tree) {\n",
    "    // setting ancestral states and placeholders\n",
    "    eta[t, node_seq[t, 1]] = eta_anc[t];\n",
    "    VCV_tips[t, node_seq[t, 1]] = diag_matrix(rep_vector(-99, J));\n",
    "    for (i in 2:N_seg) {\n",
    "      matrix[J,J] A_delta;\n",
    "      matrix[J,J] VCV;\n",
    "      vector[J] drift_seg;\n",
    "      matrix[J,J] L_VCV;\n",
    "      vector[J] A_solve_b;\n",
      "      if (length_index[t, i] > 0) {\n",
      "        A_delta = A_delta_cache[length_index[t, i]];\n",
      "        VCV = VCV_cache[length_index[t, i]];\n",
      "        L_VCV = L_VCV_tips_cache[length_index[t, i]];\n",
      "        A_solve_b = A_solve_b_cache[length_index[t, i]];\n",
    "      } else {\n",
    "        A_delta = matrix_exp(A * ts[t, i]);\n",
    "        VCV = Q_inf - quad_form_sym(Q_inf, A_delta');\n",
    "        L_VCV = cholesky_decompose(VCV);\n",
    "        A_solve_b = (A \\ add_diag(A_delta, -1)) * b;\n",
    "      }\n",
    "      drift_seg = L_VCV * z_drift[t, i-1];\n",
    "      // if not a tip, add the drift parameter\n",
    "      if (tip[t, i] == 0) {\n",
    "        eta[t, node_seq[t, i]] = to_vector(\n",
    "          A_delta * eta[t, parent[t, i]] + A_solve_b + drift_seg\n",
    "        );\n",
    "        VCV_tips[t, node_seq[t, i]] = diag_matrix(rep_vector(-99, J));\n",
    "      }\n",
    "      // if is a tip, omit, we'll deal with it in the model block;\n",
    "      else {\n",
    "        eta[t, node_seq[t, i]] = to_vector(\n",
    "          A_delta * eta[t, parent[t, i]] + A_solve_b\n",
    "        );\n",
    "        VCV_tips[t, node_seq[t, i]] = VCV;\n",
    "      }\n",
      "    }\n",
      "    }\n",
      "  }\n"
  )
  # if repeated observations or no gaussian traits, calculate tdrift
  if ((any(duplicated(data[, id])) && estimate_residual) ||
        !("normal" %in% distributions)) {
    sc_transformed_parameters <- paste0(
      sc_transformed_parameters,
      "  for (t in 1:N_tree) {\n",
      "    for (i in 1:N_tips) {\n",
      "      if (tip_to_seg[t, i] > 0 && length_index[t, tip_to_seg[t, i]] > 0) {\n",
      "        tdrift[t,i] = L_VCV_tips_cache[length_index[t, tip_to_seg[t, i]]] * ",
      "to_vector(terminal_drift[t][i,]);\n",
      "      } else {\n",
      "        tdrift[t,i] = cholesky_decompose(VCV_tips[t,i]) * ",
      "to_vector(terminal_drift[t][i,]);\n",
      "      }\n",
      "    }\n",
      "  }\n"
    )
  }
  # get code for gaussian process kernel
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
  # add gaussian process functions if dist_mat specified by user
  if (!is.null(dist_mat)) {
    sc_transformed_parameters <- paste0(
      sc_transformed_parameters,
      "  // distance covariance functions\n",
      "  for (j in 1:J) {\n",
      "    matrix[N_tips,N_tips] dist_cov;\n",
      "    matrix[N_tips,N_tips] L_dist_cov;\n",
      "    for ( i in 1:(N_tips-1) )\n",
      "      for ( m in (i+1):N_tips ) {\n",
      "        dist_cov[i,m] = ", dist_cov_code, ";\n",
      "        dist_cov[m,i] = dist_cov[i,m];\n",
      "      }\n",
      "    for ( q in 1:N_tips )\n",
      "      dist_cov[q,q] = sigma_dist[j] + 0.01;\n",
      "    L_dist_cov = cholesky_decompose(dist_cov);\n",
      "    dist_v[,j] = L_dist_cov * dist_z[,j];\n",
      "  }\n"
    )
  }
  paste0(sc_transformed_parameters, "}")
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
  sc_model <- paste0(
    "model{\n",
    "  b ~ ", priors$b, ";\n",
    "  for (t in 1:N_tree) {\n",
    "    eta_anc[t] ~ ", priors$eta_anc, ";\n",
    "    for (i in 1:(N_seg - 1)) z_drift[t, i] ~ std_normal();\n",
    # add priors for terminal_drift when:
    # 1. there are no repeated measures and no gaussian variables  OR
    # 2. there are repeated measures and estimate_residual = TRUE
    ifelse(
      (!any(duplicated(data[, id])) && !("normal" %in% distributions)) ||
        (repeated),
      "    to_vector(terminal_drift[t]) ~ std_normal();\n",
      ""
    ),
    "  }\n",
    "  A_offdiag ~ ", priors$A_offdiag, ";\n",
    "  A_diag ~ ", priors$A_diag, ";\n",
    ifelse(
      estimate_correlated_drift,
      paste0("  L_R ~ ", priors$L_R, ";\n"),
      ""
    ),
    "  Q_sigma ~ ", priors$Q_sigma, ";\n"
  )
  # add priors for any cutpoint parameters
  for (j in seq_along(distributions)) {
    if (distributions[j] == "ordered_logistic") {
      sc_model <- paste0(sc_model, "  c", j, " ~ ", priors$c, ";\n")
    }
  }
  # add priors for any overdispersion parameters
  for (j in seq_along(distributions)) {
    if (distributions[j] == "negative_binomial_softplus") {
      sc_model <-
        paste0(
          sc_model,
          "  phi", j, " ~ ",
          ifelse(
            is.null(priors$phi),
            # if prior not set manually, use default prior
            paste0("normal(inv_overdisp", j, ", inv_overdisp", j, ")"),
            # else if prior set manually, use set prior
            priors$phi
          ),
          ";\n"
        )
    }
  }
  # add priors for any shape parameters
  for (j in seq_along(distributions)) {
    if (distributions[j] == "gamma_log") {
      sc_model <- paste0(sc_model, "  shape", j, " ~ ", priors$shape, ";\n")
    }
  }
  # add priors for any gaussian process parameters
  if (!is.null(dist_mat)) {
    sc_model <- paste0(
      sc_model,
      "  to_vector(dist_z) ~ std_normal();\n",
      "  sigma_dist ~ ", priors$sigma_dist, ";\n",
      "  rho_dist ~ ", priors$rho_dist, ";\n"
    )
  }
  # add priors for any residual sds and cors
  if (repeated) {
    sc_model <- paste0(sc_model, "  // priors for residual sds and cors\n")
    if ("normal" %in% distributions) {
      sc_model <- paste0(sc_model, "  for (i in 1:N_obs) {\n")
      for (j in seq_along(distributions)) {
        sc_model <- paste0(
          sc_model,
          ifelse(
            distributions[j] == "normal",
            paste0(
              "    if (miss[i,", j, "] == 0) residual_z[", j,
              ",i] ~ std_normal();\n"
            ),
            paste0(
              "    residual_z[", j, ",i] ~ std_normal();\n"
            )
          )
        )
      }
      sc_model <- paste0(sc_model, "  }\n")
    } else {
      sc_model <- paste0(sc_model, "  to_vector(residual_z) ~ std_normal();\n")
    }
    sc_model <- paste0(
      sc_model,
      "  sigma_residual ~ ", priors$sigma_residual, ";\n",
      "  L_residual ~ ", priors$L_residual, ";\n"
    )
  }
  # add likelihood
  sc_model <- paste0(
    sc_model,
    "  if (!prior_only) {\n",
    "    for (i in 1:N_obs) {\n",
    "      vector[N_tree] lp = rep_vector(0.0, N_tree);\n",
    "      for (t in 1:N_tree) {\n",
    # only initialise tdrift vector when no repeated observations and gaussian
    ifelse(
      !(repeated) &&
        "normal" %in% distributions,
      "        vector[J] tdrift;\n",
      ""
    ),
    # only initialise residuals vector when repeated observations and gaussian
    ifelse(
      repeated &&
        ("normal" %in% distributions),
      "        vector[J] residuals;\n",
      ""
    )
  )
  # set residuals when there are repeated observations:
  if (repeated) {
    if ("normal" %in% distributions) {
      for (j in seq_along(distributions)) {
        if (distributions[j] == "normal") {
          sc_model <- paste0(
            sc_model,
            "        if (miss[i,", j, "] == 0) {\n",
            "          residuals[", j, "] = y[i,", j, "] - (", lmod(j),
            " + tdrift[t,tip_id[i]][", j, "]);\n",
            "        } else {\n",
            "          residuals[", j, "] = residual_z[", j, ",i];\n",
            "        }\n"
          )
        } else {
          sc_model <- paste0(
            sc_model,
            "        residuals[", j, "] = residual_z[", j, ",i];\n"
          )
        }
      }
      # add residual_cov matrix if measurement error included
      if (!is.null(measurement_error)) {
        sc_model <- paste0(
          sc_model,
          "        matrix[J,J] residual_cov = diag_matrix(to_vector(se[i,]))",
          " + residual_cov_base;\n"
        )
      }
      # add multi-normal prior for residuals
      sc_model <- paste0(
        sc_model,
        "        lp[t] = multi_normal_cholesky_lpdf(residuals | ",
        "rep_vector(0.0, J), ",
        ifelse(
          !is.null(measurement_error),
          "cholesky_decompose(residual_cov)",
          "L_residual_scaled"
        ),
        ");\n"
      )
    }
  } else {
    # else, set tdrift when there are no repeated observations:
    if ("normal" %in% distributions) {
      for (j in seq_along(distributions)) {
        if (distributions[j] == "normal") {
          sc_model <- paste0(
            sc_model,
            "        if (miss[i,", j, "] == 0) {\n",
            "          tdrift[", j, "] = y[i,", j, "] - (", lmod(j), ");\n",
            "          terminal_drift[t][tip_id[i],", j, "] ~ std_normal();\n",
            "        } else {\n",
            "          tdrift[", j, "] = terminal_drift[t][tip_id[i],", j,
            "];\n",
            "        }\n"
          )
        } else {
          sc_model <- paste0(
            sc_model,
            "        tdrift[", j, "] = terminal_drift[t][tip_id[i],", j, "];\n"
          )
        }
      }
      # add multi-normal prior for tdrift
      sc_model <- paste0(
        sc_model,
        "        lp[t] = multi_normal_cholesky_lpdf(tdrift | rep_vector(0.0, ",
        "J), cholesky_decompose(",
        ifelse(
          !is.null(measurement_error),
          # add squared standard errors to diagonal of VCV_tips
          "add_diag(VCV_tips[t, tip_id[i]], se[i,])",
          "VCV_tips[t, tip_id[i]]"
        ),
        "));\n"
      )
    }
  }
  # linear models for non-continuous variables
  for (j in seq_along(distributions)) {
    if (distributions[j] == "bernoulli_logit") {
      sc_model <- paste0(
        sc_model,
        "        if (miss[i,", j, "] == 0) lp[t] += ",
        "bernoulli_logit_lpmf(to_int(y[i,", j, "]) | ", lmod(j),
        ifelse(
          !(repeated) &&
            "normal" %in% distributions,
          paste0(" + tdrift[", j, "]"),
          paste0(" + tdrift[t,tip_id[i]][", j, "]")
        ),
        ifelse(
          repeated &&
            !("normal" %in% distributions),
          paste0(" + residual_v[i,", j, "]"), ""
        ),
        ifelse(
          repeated &&
            "normal" %in% distributions,
          paste0(" + residuals[", j, "]"), ""
        ),
        ");\n"
      )
    } else if (distributions[j] == "ordered_logistic") {
      sc_model <- paste0(
        sc_model,
        "        if (miss[i,", j, "] == 0) lp[t] += ",
        "ordered_logistic_lpmf(to_int(y[i,", j, "]) | ", lmod(j),
        ifelse(
          !(repeated) &&
            "normal" %in% distributions,
          paste0(" + tdrift[", j, "]"),
          paste0(" + tdrift[t,tip_id[i]][", j, "]")
        ),
        ifelse(
          repeated &&
            !("normal" %in% distributions),
          paste0(" + residual_v[i,", j, "]"), ""
        ),
        ifelse(
          repeated &&
            "normal" %in% distributions,
          paste0(" + residuals[", j, "]"), ""
        ),
        ", c", j, ");\n"
      )
    } else if (distributions[j] == "poisson_softplus") {
      sc_model <- paste0(
        sc_model,
        "        if (miss[i,", j, "] == 0) lp[t] += ",
        "poisson_lpmf(to_int(y[i,", j, "]) | mean(obs", j, ") * log1p_exp(",
        lmod(j),
        ifelse(
          !(repeated) &&
            "normal" %in% distributions,
          paste0(" + tdrift[", j, "]"),
          paste0(" + tdrift[t,tip_id[i]][", j, "]")
        ),
        ifelse(
          repeated &&
            !("normal" %in% distributions),
          paste0(" + residual_v[i,", j, "]"), ""
        ),
        ifelse(
          repeated &&
            "normal" %in% distributions,
          paste0(" + residuals[", j, "]"), ""
        ),
        "));\n"
      )
    } else if (distributions[j] == "negative_binomial_softplus") {
      sc_model <- paste0(
        sc_model,
        "        if (miss[i,", j, "] == 0) lp[t] += ",
        "neg_binomial_2_lpmf(to_int(y[i,", j,
        "]) | mean(obs", j, ") * log1p_exp(", lmod(j),
        ifelse(
          !(repeated) &&
            "normal" %in% distributions,
          paste0(" + tdrift[", j, "]"),
          paste0(" + tdrift[t,tip_id[i]][", j, "]")
        ),
        ifelse(
          repeated &&
            !("normal" %in% distributions),
          paste0(" + residual_v[i,", j, "]"), ""
        ),
        ifelse(
          repeated &&
            "normal" %in% distributions,
          paste0(" + residuals[", j, "]"), ""
        ),
        "), phi", j, ");\n"
      )
    } else if (distributions[j] == "gamma_log") {
      sc_model <- paste0(
        sc_model,
        "        if (miss[i,", j, "] == 0) lp[t] += ",
        "gamma_lpdf(y[i,", j, "] | shape", j, ", shape", j, " / exp(", lmod(j),
        ifelse(
          !(repeated) &&
            "normal" %in% distributions,
          paste0(" + tdrift[", j, "]"),
          paste0(" + tdrift[t,tip_id[i]][", j, "]")
        ),
        ifelse(
          repeated &&
            !("normal" %in% distributions),
          paste0(" + residual_v[i,", j, "]"), ""
        ),
        ifelse(
          repeated &&
            "normal" %in% distributions,
          paste0(" + residuals[", j, "]"), ""
        ),
        "));\n"
      )
    }
  }
  paste0(
    sc_model,
    "      }\n",
    "      target += log_sum_exp(lp);\n",
    "    }\n",
    "  }\n",
    "}"
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
  sc_generated_quantities <-
    paste0(
      "generated quantities{\n",
      ifelse(log_lik, "  vector[N_obs*J] log_lik; // log-likelihood\n", ""),
      "  array[N_tree,N_obs,J] real yrep; // predictive checks\n"
    )
  if (estimate_correlated_drift) {
    sc_generated_quantities <-
      paste0(
        sc_generated_quantities,
        "  matrix[J,J] cor_R; // correlated drift\n",
        "  cor_R = multiply_lower_tri_self_transpose(L_R);\n"
      )
  }
  if (repeated) {
    sc_generated_quantities <-
      paste0(
        sc_generated_quantities,
        "  matrix[J,J] cor_residual; // residual correlations\n",
        "  cor_residual = multiply_lower_tri_self_transpose(L_residual);\n"
      )
  }
  sc_generated_quantities <-
    paste0(
      sc_generated_quantities,
      "  {\n",
      ifelse(
        log_lik,
        "    matrix[N_obs,J] log_lik_temp = rep_matrix(0.0, N_obs, J);\n",
        ""
      )
    )
  # calculate terminal drift for yrep
  sc_generated_quantities <- paste0(
    sc_generated_quantities,
    "    array[N_tree,N_tips] vector[J] terminal_drift_rep;\n",
    "    for (i in 1:N_tips) {\n",
    "      for (t in 1:N_tree) {\n",
    "        for (j in 1:J) terminal_drift_rep[t,i][j] = normal_rng(0, 1);\n",
    "        terminal_drift_rep[t,i] = cholesky_decompose(",
    ifelse(
      !is.null(measurement_error) && !(any(duplicated(data[, id])) &&
                                         estimate_residual),
      # add squared standard errors to diagonal of VCV_tips
      "add_diag(VCV_tips[t, i], se[i,])",
      "VCV_tips[t, i]"
    ),
    ") * terminal_drift_rep[t,i];\n",
    "      }\n",
    "    }\n",
    "    for (i in 1:N_obs) {\n",
    ifelse(
      log_lik,
      paste0(
        "      array[N_tree,N_obs,J] real lp = ",
        "rep_array(0.0, N_tree, N_obs, J);\n"
      ),
      ""
    ),
    "      for (t in 1:N_tree) {\n"
  )
  # only declare the following if there are gaussian distributions and log_lik
  if ("normal" %in% distributions && log_lik) {
    sc_generated_quantities <- paste0(
      sc_generated_quantities,
      "        vector[J] mu_cond;\n",
      "        vector[J] sigma_cond;\n"
    )
  }
  # only declare the following if log_lik = TRUE
  if (log_lik) {
    sc_generated_quantities <- paste0(
      sc_generated_quantities,
      ifelse(
        repeated,
        "        vector[J] residuals;\n",
        "        vector[J] tdrifts;\n"
      )
    )
  }
  # calculate residuals_rep for yrep if repeated observations
  if (repeated) {
    sc_generated_quantities <- paste0(
      sc_generated_quantities,
      "        vector[J] residuals_rep;\n",
      "        for (j in 1:J) residuals_rep[j] = normal_rng(0, 1);\n",
      "        residuals_rep = L_residual_scaled * residuals_rep;\n"
    )
  }
  # get tdrifts/residuals if log_lik = TRUE
  if (log_lik) {
    for (j in seq_along(distributions)) {
      if (distributions[j] == "normal") {
        sc_generated_quantities <- paste0(
          sc_generated_quantities,
          ifelse(
            repeated,
            paste0(
              "        residuals[", j, "] = y[i][", j, "] - (", lmod(j),
              " + tdrift[t,tip_id[i]][", j, "]);\n"
            ),
            paste0(
              "        tdrifts[", j, "] = y[i][", j, "] - (", lmod(j), ");\n"
            )
          )
        )
      } else {
        sc_generated_quantities <- paste0(
          sc_generated_quantities,
          ifelse(
            repeated,
            ifelse(
              "normal" %in% distributions,
              paste0("        residuals[", j, "] = residual_z[", j, ",i];\n"),
              paste0("        residuals[", j, "] = residual_v[i,", j, "];\n")
            ),
            paste0(
              "        tdrifts[", j, "] = ",
              ifelse("normal" %in% distributions, "terminal_drift", "tdrift"),
              "[t,tip_id[i]][", j, "];\n"
            )
          )
        )
      }
    }
  }
  # only calculate if there are gaussian distributions and log_lik = TRUE
  if ("normal" %in% distributions && log_lik) {
    if (repeated) {
      sc_generated_quantities <- paste0(
        sc_generated_quantities,
        ifelse(
          !is.null(measurement_error),
          paste0(
            "        matrix[J,J] residual_cov = ",
            "diag_matrix(to_vector(se[i,])) + residual_cov_base;\n"
          ),
          ""
        ),
        ifelse(
          !is.null(measurement_error),
          "        matrix[J,J] cov_inv = inverse_spd(residual_cov);\n",
          paste0(
            "        matrix[J,J] cov_inv = chol2inv(L_residual_scaled);\n"
          )
        ),
        "        mu_cond = residuals - (cov_inv * residuals) ./ ",
        "diagonal(cov_inv);\n",
        "        sigma_cond = sqrt(1 / diagonal(cov_inv));\n"
      )
    } else {
      sc_generated_quantities <- paste0(
        sc_generated_quantities,
        "        matrix[J,J] cov_inv = inverse_spd(",
        ifelse(
          !is.null(measurement_error),
          # add squared standard errors to diagonal of VCV_tips
          "add_diag(VCV_tips[t, tip_id[i]], se[i,])",
          "VCV_tips[t, tip_id[i]]"
        ),
        ");\n",
        "        mu_cond = tdrifts - (cov_inv * tdrifts) ./ ",
        "diagonal(cov_inv);\n",
        "        sigma_cond = sqrt(1 / diagonal(cov_inv));\n"
      )
    }
  }
  for (j in seq_along(distributions)) {
    if (distributions[j] == "bernoulli_logit") {
      sc_generated_quantities <- paste0(
        sc_generated_quantities,
        ifelse(
          log_lik,
          paste0(
            "        if (miss[i,", j, "] == 0) lp[t,i,", j, "] = ",
            "bernoulli_logit_lpmf(to_int(y[i,", j, "]) | ", lmod(j),
            ifelse(
              repeated,
              paste0(" + residuals[", j, "]"),
              paste0(" + tdrifts[", j, "]")
            ),
            ");\n"
          ),
          ""
        ),
        "        yrep[t,i,", j, "] = bernoulli_logit_rng(", lmod(j),
        " + terminal_drift_rep[t,tip_id[i]][", j, "]",
        ifelse(
          repeated,
          paste0(" + residuals_rep[", j, "]"),
          ""
        ),
        ");\n"
      )
    } else if (distributions[j] == "ordered_logistic") {
      sc_generated_quantities <- paste0(
        sc_generated_quantities,
        ifelse(
          log_lik,
          paste0(
            "        if (miss[i,", j, "] == 0) lp[t,i,", j, "] = ",
            "ordered_logistic_lpmf(to_int(y[i,", j, "]) | ", lmod(j),
            ifelse(
              repeated,
              paste0(" + residuals[", j, "]"),
              paste0(" + tdrifts[", j, "]")
            ),
            ", c", j, ");\n"
          ),
          ""
        ),
        "        yrep[t,i,", j, "] = ", "ordered_logistic_rng(", lmod(j),
        " + terminal_drift_rep[t,tip_id[i]][", j, "]",
        ifelse(
          repeated,
          paste0(" + residuals_rep[", j, "]"),
          ""
        ),
        ", c", j, ");\n"
      )
    } else if (distributions[j] == "poisson_softplus") {
      sc_generated_quantities <- paste0(
        sc_generated_quantities,
        ifelse(
          log_lik,
          paste0(
            "        if (miss[i,", j, "] == 0) lp[t,i,", j, "] = ",
            "poisson_lpmf(to_int(y[i,", j, "]) | mean(obs", j,
            ") * log1p_exp(", lmod(j),
            ifelse(
              repeated,
              paste0(" + residuals[", j, "]"),
              paste0(" + tdrifts[", j, "]")
            ),
            "));\n"
          ),
          ""
        ),
        "        yrep[t,i,", j, "] = poisson_rng(mean(obs", j, ") * log1p_exp(",
        lmod(j), " + terminal_drift_rep[t,tip_id[i]][", j, "]",
        ifelse(
          repeated,
          paste0(" + residuals_rep[", j, "]"),
          ""
        ),
        "));\n"
      )
    } else if (distributions[j] == "normal") {
      sc_generated_quantities <- paste0(
        sc_generated_quantities,
        ifelse(
          log_lik,
          paste0(
            "        if (miss[i,", j, "] == 0) lp[t,i,", j, "] = ",
            "normal_lpdf(",
            ifelse(
              repeated,
              paste0("residuals[", j, "]"),
              paste0("tdrifts[", j, "]")
            ),
            " | mu_cond[", j,
            "], sigma_cond[", j, "]);\n"
          ),
          ""
        ),
        "        yrep[t,i,", j, "] = ", lmod(j),
        " + terminal_drift_rep[t,tip_id[i]][", j, "]",
        ifelse(
          repeated,
          paste0(" + residuals_rep[", j, "]"),
          ""
        ),
        ";\n"
      )
    } else if (distributions[j] == "negative_binomial_softplus") {
      sc_generated_quantities <- paste0(
        sc_generated_quantities,
        ifelse(
          log_lik,
          paste0(
            "        if (miss[i,", j, "] == 0) lp[t,i,", j, "] = ",
            "neg_binomial_2_lpmf(to_int(y[i,", j, "]) | mean(obs", j,
            ") * log1p_exp(", lmod(j),
            ifelse(
              repeated,
              paste0(" + residuals[", j, "]"),
              paste0(" + tdrifts[", j, "]")
            ),
            "), phi", j, ");\n"
          ),
          ""
        ),
        "        yrep[t,i,", j, "] = neg_binomial_2_rng(mean(obs", j,
        ") * log1p_exp(", lmod(j), " + terminal_drift_rep[t,tip_id[i]][", j,
        "]",
        ifelse(
          repeated,
          paste0(" + residuals_rep[", j, "]"),
          ""
        ),
        "), phi", j, ");\n"
      )
    } else if (distributions[j] == "gamma_log") {
      sc_generated_quantities <- paste0(
        sc_generated_quantities,
        ifelse(
          log_lik,
          paste0(
            "        if (miss[i,", j, "] == 0) lp[t,i,", j, "] = ",
            "gamma_lpdf(y[i,", j, "] | shape", j, ", shape", j, " / exp(",
            lmod(j),
            ifelse(
              repeated,
              paste0(" + residuals[", j, "]"),
              paste0(" + tdrifts[", j, "]")
            ),
            "));\n"
          ),
          ""
        ),
        "        yrep[t,i,", j, "] = gamma_rng(shape", j, ", shape", j,
        " / exp(", lmod(j), " + terminal_drift_rep[t,tip_id[i]][", j, "]",
        ifelse(
          repeated,
          paste0(" + residuals_rep[", j, "]"),
          ""
        ),
        "));\n"
      )
    }
  }
  paste0(
    sc_generated_quantities,
    "      }\n",
    ifelse(
      log_lik,
      "    for (j in 1:J) log_lik_temp[i,j] += log_sum_exp(lp[,i,j]);\n",
      ""
    ),
    "    }\n",
    ifelse(
      log_lik,
      "  log_lik = to_vector(log_lik_temp);\n",
      ""
    ),
    "  }\n",
    "}"
  )
}
