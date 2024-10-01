#' Make Stan code for dynamic coevolutionary model
#'
#' Make the \pkg{Stan} code for the Bayesian dynamic coevolutionary model.
#' \pkg{Stan} code is generated, checked for syntactical errors, and then
#' returned as a character string.
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
#' @param tree A phylogenetic tree object of class \code{phylo} or
#'   \code{multiPhylo}. The tree(s) must be rooted and must include positive
#'   non-zero branch lengths.
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
#'   specified, the model uses default priors (see \code{help(coev_fit)}).
#'   Alternatively, the user can specify a named list of priors. The list must
#'   contain non-duplicated entries for any of the following parameters: the
#'   autoregressive effects (\code{A_diag}), the cross effects
#'   (\code{A_offdiag}), the Cholesky factor for the drift matrix (\code{L_R}),
#'   the drift std. dev. parameters (\code{Q_sigma}), the continuous time
#'   intercepts (\code{b}), the ancestral states for the traits
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
#' @param estimate_Q_offdiag Logical. If \code{TRUE} (default), the model
#'   estimates the off-diagonals for the \deqn{Q} drift matrix (i.e., correlated
#'   drift). If \code{FALSE}, the off-diagonals for the \deqn{Q} drift matrix
#'   are set to zero.
#' @param prior_only Logical. If \code{FALSE} (default), the model is fitted to
#'   the data and returns a posterior distribution. If \code{TRUE}, the model
#'   samples from the prior only, ignoring the likelihood.
#'
#' @return A character string containing the \pkg{Stan} code to fit the dynamic
#'   coevolutionary model.
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
#'   )
#'
#' # print Stan code
#' cat(stan_code)
#'
#' @export
coev_make_stancode <- function(data, variables, id, tree,
                               effects_mat = NULL, dist_mat = NULL,
                               prior = NULL, scale = TRUE,
                               estimate_Q_offdiag = TRUE,
                               prior_only = FALSE) {
  # check arguments
  run_checks(data, variables, id, tree, effects_mat,
             dist_mat, prior, scale, estimate_Q_offdiag, prior_only)
  # coerce data argument to data frame
  data <- as.data.frame(data)
  # extract distributions and variable names from named list
  distributions <- as.character(variables)
  variables <- names(variables)
  # get default priors
  priors <-
    list(
      b           = "std_normal()",
      eta_anc     = "std_normal()",
      A_offdiag   = "std_normal()",
      A_diag      = "std_normal()",
      L_R         = "lkj_corr_cholesky(4)",
      Q_sigma     = "std_normal()",
      c           = "normal(0, 2)",
      nu          = "gamma(2, 0.1)",
      sigma_dist  = "exponential(1)",
      rho_dist    = "exponential(1)",
      sigma_group = "exponential(1)",
      L_group     = "lkj_corr_cholesky(2)"
    )
  # note: default prior for phi (overdispersion) set within the model code
  # replace priors if user has explicitly set them
  if (!is.null(prior)) {
    for (i in names(prior)) {
      priors[[i]] <- prior[[i]]
    }
  }
  # write functions block
  sc_functions <- paste0(
    "functions {\n",
    "  // Charles Driver's optimized way of solving for the asymptotic Q matrix\n",
    "  matrix ksolve (matrix A, matrix Q) {\n",
    "    int d = rows(A);\n",
    "    int d2 = (d * d - d) %/% 2;\n",
    "    matrix [d + d2, d + d2] O;\n",
    "    vector [d + d2] triQ;\n",
    "    matrix[d,d] AQ;\n",
    "    int z = 0;         // z is row of output\n",
    "    for (j in 1:d) {   // for column reference of solution vector\n",
    "      for (i in 1:j) { // and row reference...\n",
    "        if (j >= i) {  // if i and j denote a covariance parameter (from upper tri)\n",
    "          int y = 0;   // start new output row\n",
    "          z += 1;      // shift current output row down\n",
    "          for (ci in 1:d) {   // for columns and\n",
    "            for (ri in 1:d) { // rows of solution\n",
    "              if (ci >= ri) { // when in upper tri (inc diag)\n",
    "                y += 1;       // move to next column of output\n",
    "                if (i == j) { // if output row is for a diagonal element\n",
    "                  if (ri == i) O[z, y] = 2 * A[ri, ci];\n",
    "                  if (ci == i) O[z, y] = 2 * A[ci, ri];\n",
    "                }\n",
    "                if (i != j) { // if output row is not for a diagonal element\n",
    "                  //if column of output matches row of output, sum both A diags\n",
    "                  if (y == z) O[z, y] = A[ri, ri] + A[ci, ci];\n",
    "                  if (y != z) { // otherwise...\n",
    "                    // if solution element we refer to is related to output row...\n",
    "                    if (ci == ri) { // if solution element is a variance\n",
    "                      // if variance of solution corresponds to row of our output\n",
    "                      if (ci == i) O[z, y] = A[j, ci];\n",
    "                      // if variance of solution corresponds to col of our output\n",
    "                      if (ci == j) O[z, y] = A[i, ci];\n",
    "                    }\n",
    "                    //if solution element is a related covariance\n",
    "                    if (ci != ri && (ri == i || ri == j || ci == i || ci == j )) {\n",
    "                      // for row 1,2 / 2,1 of output, if solution row ri 1 (match)\n",
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
  # write data block
  sc_data <- paste0(
    "data{\n",
    "  int<lower=1> N_tips; // number of tips\n",
    "  int<lower=1> N_tree; // number of trees\n",
    "  int<lower=1> N_obs; // number of observations\n",
    "  int<lower=2> J; // number of response traits\n",
    "  int<lower=1> N_seg; // total number of segments in the trees\n",
    "  array[N_tree, N_seg] int<lower=1> node_seq; // index of tree nodes\n",
    "  array[N_tree, N_seg] int<lower=0> parent; // index of the parent node of each descendent\n",
    "  array[N_tree, N_seg] real ts; // time since parent\n",
    "  array[N_tree, N_seg] int<lower=0,upper=1> tip; // indicator of whether a given segment ends in a tip\n",
    "  array[J,J] int<lower=0,upper=1> effects_mat; // which effects should be estimated?\n",
    "  int<lower=2> num_effects; // number of effects being estimated\n",
    "  matrix[N_obs,J] y; // observed data\n",
    "  matrix[N_obs,J] miss; // are data points missing?\n",
    "  array[N_obs] int<lower=1> tip_id; // index between 1 and N_tips that gives the group id\n"
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
  sc_data <-
    paste0(
      sc_data,
      "  int<lower=0,upper=1> prior_only; // should the likelihood be ignored?\n}"
      )
  # write transformed data block
  sc_transformed_data <- "transformed data{\n"
  for (j in 1:length(distributions)) {
    sc_transformed_data <- paste0(
      sc_transformed_data,
      "  vector[to_int(N_obs - sum(col(miss, ", j,
      ")))] obs", j, "; // observed data for variable ", j, "\n"
    )
  }
  for (j in 1:length(distributions)) {
    if (distributions[j] == "negative_binomial_softplus" & is.null(priors$phi)){
      sc_transformed_data <-
        paste0(
          sc_transformed_data,
          "  real inv_overdisp", j, "; // best guess for phi", j, "\n"
        )
    }
  }
  for (j in 1:length(distributions)) {
    sc_transformed_data <- paste0(
      sc_transformed_data,
      "  obs", j, " = col(y, ", j,
      ")[which_equal(col(miss, ", j, "), 0)];\n"
    )
  }
  for (j in 1:length(distributions)) {
    if (distributions[j] == "negative_binomial_softplus" & is.null(priors$phi)){
      sc_transformed_data <-
        paste0(
          sc_transformed_data,
          "  inv_overdisp", j, " = (mean(obs", j, ")^2) / ",
          "(sd(obs", j, ")^2 - mean(obs", j, "));\n"
        )
    }
  }
  sc_transformed_data <- paste0(sc_transformed_data, "}")
  # write parameters block
  sc_parameters <- paste0(
    "parameters{\n",
    "  vector<upper=0>[J] A_diag; // autoregressive terms of A\n",
    "  vector[num_effects - J] A_offdiag; // cross-lagged terms of A\n"
  )
  # add cholesky factor for Q matrix if estimating Q off diagonals
  if (estimate_Q_offdiag) {
    sc_parameters <- paste0(
      sc_parameters,
      "  cholesky_factor_corr[J] L_R; // lower-tri choleksy decomp corr mat, used to construct Q mat\n"
    )
  }
  sc_parameters <- paste0(
    sc_parameters,
    "  vector<lower=0>[J] Q_sigma; // std deviation parameters of the Q mat\n",
    "  vector[J] b; // SDE intercepts\n",
    "  array[N_tree] vector[J] eta_anc; // ancestral states\n",
    "  array[N_tree, N_seg - 1] vector[J] z_drift; // stochastic drift, unscaled and uncorrelated\n"
  )
  for (i in 1:length(distributions)) {
    # add cut points for ordinal_logistic distributions
    if (distributions[i] == "ordered_logistic") {
      # calculate number of cut points (number of levels - 1)
      num_cuts <- max(as.numeric(data[,variables[i]]), na.rm = TRUE) - 1
      sc_parameters <-
        paste0(
          sc_parameters,
          "  ordered[", num_cuts, "] c", i, "; ",
          "// cut points for variable ", i, "\n"
          )
    }
    # add degrees of freedom for student_t distributions
    if (distributions[i] == "student_t") {
      sc_parameters <-
        paste0(
          sc_parameters,
          "  real<lower=1> nu", i, "; ",
          "// student t degrees of freedom for variable ", i, "\n"
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
  }
  # add gaussian process parameters if distance matrix specified by user
  if (!is.null(dist_mat)) {
    sc_parameters <- paste0(
      sc_parameters,
      "  matrix[N_tips,J] dist_z; // spatial covariance random effects, unscaled and uncorrelated\n",
      "  vector<lower=0>[J] rho_dist; // how quickly does covariance decline with distance\n",
      "  vector<lower=0>[J] sigma_dist; // maximum covariance due to spatial location\n"
    )
  }
  # add multilevel parameters if there are repeated measures
  if (any(duplicated(data[,id]))) {
    sc_parameters <- paste0(
      sc_parameters,
      "  matrix[J,N_tips] group_z;\n",
      "  vector<lower=0>[J] sigma_group;\n",
      "  cholesky_factor_corr[J] L_group;\n"
    )
  }
  sc_parameters <- paste0(sc_parameters, "}")
  # write transformed parameters block
  sc_transformed_parameters <- paste0(
    "transformed parameters{\n",
    "  array[N_tree, N_seg] vector[J] eta;\n",
    "  matrix[J,J] A = diag_matrix(A_diag); // selection matrix\n",
    ifelse(
      estimate_Q_offdiag,
      "  matrix[J,J] Q = diag_matrix(Q_sigma) * (L_R * L_R') * diag_matrix(Q_sigma); // drift matrix\n",
      "  matrix[J,J] Q = diag_matrix(Q_sigma^2); // drift matrix\n"
    ),
    "  matrix[J,J] Q_inf; // asymptotic covariance matrix\n",
    "  array[N_tree, N_seg] vector[J] drift_tips; // terminal drift parameters\n",
    "  array[N_tree, N_seg] vector[J] sigma_tips; // terminal drift parameters\n"
  )
  # add distance random effects if distance matrix specified by user
  if (!is.null(dist_mat)) {
    sc_transformed_parameters <- paste0(
      sc_transformed_parameters,
      "  matrix[N_tips,J] dist_v; // distance covariance random effects\n"
    )
  }
  # add group random effects if there are repeated measures
  if (any(duplicated(data[,id]))) {
    sc_transformed_parameters <- paste0(
      sc_transformed_parameters,
      "  matrix[N_tips,J] group_v; // group random effects\n",
      "  // scale and correlate group random effects\n",
      "  group_v = (diag_pre_multiply(sigma_group, L_group) * group_z)';\n"
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
    "  // loop over phylogenetic trees\n",
    "  for (t in 1:N_tree) {\n",
    "    // setting ancestral states and placeholders\n",
    "    eta[t, node_seq[t, 1]] = eta_anc[t];\n",
    "    drift_tips[t, node_seq[t, 1]] = rep_vector(-99, J);\n",
    "    sigma_tips[t, node_seq[t, 1]] = rep_vector(-99, J);\n",
    "    for (i in 2:N_seg) {\n",
    "      matrix[J,J] A_delta; // amount of deterministic change (selection)\n",
    "      matrix[J,J] VCV; // variance-covariance matrix of stochastic change (drift)\n",
    "      vector[J] drift_seg; // accumulated drift over the segment\n",
    "      A_delta = matrix_exp(A * ts[t, i]);\n",
    "      VCV = Q_inf - quad_form_sym(Q_inf, A_delta');\n",
    "      drift_seg = cholesky_decompose(VCV) * z_drift[t, i-1];\n",
    "      // if not a tip, add the drift parameter\n",
    "      if (tip[t, i] == 0) {\n",
    "        eta[t, node_seq[t, i]] = to_vector(\n",
    "          A_delta * eta[t, parent[t, i]] + ((A \\ add_diag(A_delta, -1)) * b) + drift_seg\n",
    "        );\n",
    "        drift_tips[t, node_seq[t, i]] = rep_vector(-99, J);\n",
    "        sigma_tips[t, node_seq[t, i]] = rep_vector(-99, J);\n",
    "      }\n",
    "      // if is a tip, omit, we'll deal with it in the model block;\n",
    "      else {\n",
    "        eta[t, node_seq[t, i]] = to_vector(\n",
    "          A_delta * eta[t, parent[t, i]] + ((A \\ add_diag(A_delta, -1)) * b)\n",
    "        );\n",
    "        drift_tips[t, node_seq[t, i]] = drift_seg;\n",
    "        sigma_tips[t, node_seq[t, i]] = diagonal(Q);\n",
    "      }\n",
    "    }\n",
    "  }\n"
  )
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
      "        dist_cov[i,m] = sigma_dist[j]*exp(-(rho_dist[j]*dist_mat[i,m]));\n",
      "        dist_cov[m,i] = dist_cov[i,m];\n",
      "      }\n",
      "    for ( q in 1:N_tips )\n",
      "      dist_cov[q,q] = sigma_dist[j] + 0.01;\n",
      "    L_dist_cov = cholesky_decompose(dist_cov);\n",
      "    dist_v[,j] = L_dist_cov * dist_z[,j];\n",
      "  }\n"
    )
  }
  sc_transformed_parameters <- paste0(sc_transformed_parameters, "}")
  # write model block
  sc_model <- paste0(
    "model{\n",
    "  b ~ ", priors$b, ";\n",
    "  for (t in 1:N_tree) {\n",
    "    eta_anc[t] ~ ", priors$eta_anc, ";\n",
    "    for (i in 1:(N_seg - 1)) z_drift[t, i] ~ std_normal();\n",
    "  }\n",
    "  A_offdiag ~ ", priors$A_offdiag, ";\n",
    "  A_diag ~ ", priors$A_diag, ";\n",
    ifelse(estimate_Q_offdiag, paste0("  L_R ~ ", priors$L_R, ";\n"), ""),
    "  Q_sigma ~ ", priors$Q_sigma, ";\n"
  )
  # add priors for any cut points
  for (j in 1:length(distributions)) {
    if (distributions[j] == "ordered_logistic") {
      sc_model <- paste0(sc_model, "  c", j, " ~ ", priors$c, ";\n")
    }
  }
  # add priors for any degrees of freedom
  for (j in 1:length(distributions)) {
    if (distributions[j] == "student_t") {
      sc_model <- paste0(sc_model, "  nu", j, " ~ ", priors$nu, ";\n")
    }
  }
  # add priors for any overdispersion parameters
  for (j in 1:length(distributions)) {
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
  # add priors for any gaussian process parameters
  if (!is.null(dist_mat)) {
    sc_model <- paste0(
      sc_model,
      "  to_vector(dist_z) ~ std_normal();\n",
      "  sigma_dist ~ ", priors$sigma_dist, ";\n",
      "  rho_dist ~ ", priors$rho_dist, ";\n"
    )
  }
  # add priors for any group-level random effect parameters
  if (any(duplicated(data[,id]))) {
    sc_model <- paste0(
      sc_model,
      "  // priors for group-level random effects (non-phylogenetic)\n",
      "  to_vector(group_z) ~ std_normal();\n",
      "  sigma_group ~ ", priors$sigma_group, ";\n",
      "  L_group ~ ", priors$L_group, ";\n"
    )
  }
  # add likelihood
  sc_model <- paste0(
    sc_model,
    "  if (!prior_only) {\n",
    "    for (i in 1:N_obs) {\n",
    "      array[N_tree,J] real lp = rep_array(log(1.0 / N_tree), N_tree, J);\n",
    "      for (t in 1:N_tree) {\n"
    )
  # function to get linear model
  lmod <- function(j) {
    paste0(
      "eta[t,tip_id[i]][", j, "]",
      ifelse(
        !is.null(dist_mat),
        paste0(" + dist_v[tip_id[i],", j, "]"),
        ""
      ),
      ifelse(
        any(duplicated(data[,id])),
        paste0(" + group_v[tip_id[i],", j, "]"),
        ""
      )
    )
  }
  for (j in 1:length(distributions)) {
    if (distributions[j] == "bernoulli_logit") {
      sc_model <- paste0(
        sc_model,
        "        if (miss[i,", j, "] == 0) lp[t,", j, "] += ",
        "bernoulli_logit_lpmf(to_int(y[i,", j, "]) | ", lmod(j),
        " + drift_tips[t,tip_id[i]][", j, "]);\n"
        )
    } else if (distributions[j] == "ordered_logistic") {
      sc_model <- paste0(
        sc_model,
        "        if (miss[i,", j, "] == 0) lp[t,", j, "] += ",
        "ordered_logistic_lpmf(to_int(y[i,", j, "]) | ", lmod(j),
        " + drift_tips[t,tip_id[i]][", j, "], c", j, ");\n"
        )
    } else if (distributions[j] == "poisson_softplus") {
      sc_model <- paste0(
        sc_model,
        "        if (miss[i,", j, "] == 0) lp[t,", j, "] += ",
        "poisson_lpmf(to_int(y[i,", j, "]) | mean(obs", j, ") * log1p_exp(",
        lmod(j), " + drift_tips[t,tip_id[i]][", j, "]));\n"
        )
    } else if (distributions[j] == "normal") {
      sc_model <- paste0(
        sc_model,
        "        if (miss[i,", j, "] == 0) lp[t,", j, "] += ",
        "normal_lpdf(y[i,", j, "] | ", lmod(j),
        ", sigma_tips[t,tip_id[i]][", j, "]);\n"
        )
    } else if (distributions[j] == "student_t") {
      sc_model <- paste0(
        sc_model,
        "        if (miss[i,", j, "] == 0) lp[t,", j, "] += ",
        "student_t_lpdf(y[i,", j, "] | nu", j, ", ", lmod(j),
        ", sigma_tips[t,tip_id[i]][", j, "]);\n"
      )
    } else if (distributions[j] == "lognormal") {
      sc_model <- paste0(
        sc_model,
        "        if (miss[i,", j, "] == 0) lp[t,", j, "] += ",
        "lognormal_lpdf(y[i,", j, "] | ", lmod(j),
        ", sigma_tips[t,tip_id[i]][", j, "]);\n"
        )
    } else if (distributions[j] == "negative_binomial_softplus") {
      sc_model <- paste0(
        sc_model,
        "        if (miss[i,", j, "] == 0) lp[t,", j, "] += ",
        "neg_binomial_2_lpmf(to_int(y[i,", j,
        "]) | mean(obs", j, ") * log1p_exp(", lmod(j),
        " + drift_tips[t,tip_id[i]][", j, "]), phi", j, ");\n"
        )
    }
  }
  sc_model <- paste0(
    sc_model,
    "      }\n",
    "      for (j in 1:J) target += log_sum_exp(lp[,j]);\n",
    "    }\n",
    "  }\n",
    "}"
    )
  # generated quantities block
  sc_generated_quantities <-
    paste0(
      "generated quantities{\n",
      "  vector[N_obs*J] log_lik; // log-likelihood\n",
      "  array[N_tree,N_obs,J] real yrep; // predictive checks\n"
      )
  if (estimate_Q_offdiag) {
    sc_generated_quantities <-
      paste0(
        sc_generated_quantities,
        "  matrix[J,J] cor_R; // correlated drift\n",
        "  cor_R = multiply_lower_tri_self_transpose(L_R);\n"
      )
  }
  if (any(duplicated(data[,id]))) {
    sc_generated_quantities <-
      paste0(
        sc_generated_quantities,
        "  matrix[J,J] cor_group; // group-level correlations\n",
        "  cor_group = multiply_lower_tri_self_transpose(L_group);\n"
      )
  }
  sc_generated_quantities <-
    paste0(
      sc_generated_quantities,
      "  {\n",
      "    matrix[N_obs,J] log_lik_temp = rep_matrix(0, N_obs, J);\n",
      "    for (i in 1:N_obs) {\n",
      "      array[N_tree,N_obs,J] real lp = rep_array(log(1.0 / N_tree), N_tree, N_obs, J);\n",
      "      for (t in 1:N_tree) {\n"
      )
  for (j in 1:length(distributions)) {
    if (distributions[j] == "bernoulli_logit") {
      sc_generated_quantities <- paste0(
        sc_generated_quantities,
        "        if (miss[i,", j, "] == 0) lp[t,i,", j, "] += ",
        "bernoulli_logit_lpmf(to_int(y[i,", j, "]) | ", lmod(j),
        " + drift_tips[t,tip_id[i]][", j, "]);\n",
        "        yrep[t,i,", j, "] = ",
        "bernoulli_logit_rng(", lmod(j),
        " + drift_tips[t,tip_id[i]][", j, "]);\n"
      )
    } else if (distributions[j] == "ordered_logistic") {
      sc_generated_quantities <- paste0(
        sc_generated_quantities,
        "        if (miss[i,", j, "] == 0) lp[t,i,", j, "] += ",
        "ordered_logistic_lpmf(to_int(y[i,", j, "]) | ", lmod(j),
        " + drift_tips[t,tip_id[i]][", j, "], c", j, ");\n",
        "        yrep[t,i,", j, "] = ",
        "ordered_logistic_rng(", lmod(j),
        " + drift_tips[t,tip_id[i]][", j, "], c", j, ");\n"
      )
    } else if (distributions[j] == "poisson_softplus") {
      sc_generated_quantities <- paste0(
        sc_generated_quantities,
        "        if (miss[i,", j, "] == 0) lp[t,i,", j, "] += ",
        "poisson_lpmf(to_int(y[i,", j, "]) | mean(obs", j,
        ") * log1p_exp(", lmod(j),
        " + drift_tips[t,tip_id[i]][", j, "]));\n",
        "        yrep[t,i,", j, "] = ",
        "poisson_rng(mean(obs", j, ") * log1p_exp(", lmod(j),
        " + drift_tips[t,tip_id[i]][", j, "]));\n"
      )
    } else if (distributions[j] == "normal") {
      sc_generated_quantities <- paste0(
        sc_generated_quantities,
        "        if (miss[i,", j, "] == 0) lp[t,i,", j, "] += ",
        "normal_lpdf(y[i,", j, "] | ", lmod(j),
        ", sigma_tips[t,tip_id[i]][", j, "]);\n",
        "        yrep[t,i,", j, "] = ",
        "normal_rng(", lmod(j),
        ", sigma_tips[t,tip_id[i]][", j, "]);\n"
      )
    } else if (distributions[j] == "student_t") {
      sc_generated_quantities <- paste0(
        sc_generated_quantities,
        "        if (miss[i,", j, "] == 0) lp[t,i,", j, "] += ",
        "student_t_lpdf(y[i,", j, "] | nu", j, ", ", lmod(j),
        ", sigma_tips[t,tip_id[i]][", j, "]);\n",
        "        yrep[t,i,", j, "] = ",
        "student_t_rng(nu", j, ", ", lmod(j),
        ", sigma_tips[t,tip_id[i]][", j, "]);\n"
      )
    } else if (distributions[j] == "lognormal") {
      sc_generated_quantities <- paste0(
        sc_generated_quantities,
        "        if (miss[i,", j, "] == 0) lp[t,i,", j, "] += ",
        "lognormal_lpdf(y[i,", j, "] | ", lmod(j),
        ", sigma_tips[t,tip_id[i]][", j, "]);\n",
        "        yrep[t,i,", j, "] = ",
        "lognormal_rng(", lmod(j),
        ", sigma_tips[t,tip_id[i]][", j, "]);\n"
      )
    } else if (distributions[j] == "negative_binomial_softplus") {
      sc_generated_quantities <- paste0(
        sc_generated_quantities,
        "        if (miss[i,", j, "] == 0) lp[t,i,", j, "] += ",
        "neg_binomial_2_lpmf(to_int(y[i,", j, "]) | mean(obs", j,
        ") * log1p_exp(", lmod(j),
        " + drift_tips[t,tip_id[i]][", j, "]), phi", j, ");\n",
        "        yrep[t,i,", j, "] = ",
        "neg_binomial_2_rng(mean(obs", j, ") * log1p_exp(", lmod(j),
        " + drift_tips[t,tip_id[i]][", j, "]), phi", j, ");\n"
      )
    }
  }
  sc_generated_quantities <-
    paste0(
      sc_generated_quantities,
      "      }\n",
      "    for (j in 1:J) log_lik_temp[i,j] += log_sum_exp(lp[,i,j]);\n",
      "    }\n",
      "  log_lik = to_vector(log_lik_temp);\n",
      "  }\n",
      "}"
    )
  # put stan code together
  sc <- paste0(
    "// Generated with coevolve ", utils::packageVersion("coevolve"), "\n",
    sc_functions, "\n",
    sc_data, "\n",
    ifelse(
      is.null(sc_transformed_data),
      "",
      paste0(sc_transformed_data, "\n")
      ),
    sc_parameters, "\n",
    sc_transformed_parameters, "\n",
    sc_model, "\n",
    sc_generated_quantities
  )
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
        "distances have been included for each variable in the model."
        )
    )
  }
  if (any(duplicated(data[,id]))) {
    message(
      paste0(
        "Note: Repeated observations detected. Group-level varying effects ",
        "have been included for each variable in the model."
      )
    )
  }
  # return stan code
  return(sc)
}
