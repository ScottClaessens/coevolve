#' Make Stan code for dynamic coevolutionary model
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
#' @param prior_only Logical. If \code{FALSE} (default), the model is fitted to
#'   the data and returns a posterior distribution. If \code{TRUE}, the model
#'   samples from the prior only, ignoring the likelihood.
#'
#' @return A character string containing the \pkg{Stan} code to fit the dynamic
#'   coevolutionary model.
#' @export
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
coev_make_stancode <- function(data, variables, id, tree,
                               effects_mat = NULL, dist_mat = NULL,
                               prior = NULL, prior_only = FALSE) {
  # check arguments
  run_checks(data, variables, id, tree, effects_mat,
             dist_mat, prior, prior_only)
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
      Q_diag      = "std_normal()",
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
    "  // returns the Kronecker sum of a square matrix with itself\n",
    "  matrix kronecker_sum_self(matrix A) {\n",
    "    int n = rows(A);\n",
    "    matrix[n*n, n*n] result = rep_matrix(0, n*n, n*n);\n",
    "    // fill the diagonal blocks\n",
    "    for (i in 1:n) {\n",
    "      int start_row = (i - 1) * n + 1;\n",
    "      int end_row = i * n;\n",
    "      result[start_row:end_row, start_row:end_row] = add_diag(A, 1) * A[i, i];\n",
    "    }\n",
    "    // fill the off-diagonal blocks\n",
    "    for (i in 1:n) {\n",
    "      for (j in 1:n) {\n",
    "        if (i != j) {\n",
    "          int start_row = (i - 1) * n + 1;\n",
    "          int end_row = i * n;\n",
    "          int start_col = (j - 1) * n + 1;\n",
    "          int end_col = j * n;\n",
    "          result[start_row:end_row, start_col:end_col] = A[i, j] * identity_matrix(n);\n",
    "        }\n",
    "      }\n",
    "    }\n",
    "    return result;\n",
    "  }\n",
    "  \n",
    "  // solve SDE\n",
    "  matrix cov_drift(matrix A, matrix Q, real ts) {\n",
    "    matrix[rows(A) * rows(A), cols(A) * cols(A)] A_sharp_temp;\n",
    "    vector[rows(Q)*cols(Q)] row_Q;\n",
    "    vector[rows(A)*cols(A)] irow_vec;\n",
    "    matrix[rows(A),cols(A)] irow_mat;\n",
    "    A_sharp_temp = kronecker_sum_self(A);\n",
    "    // row operation takes elements of a matrix rowwise and puts them into a column vector\n",
    "    for (i in 1:rows(Q))\n",
    "      for (j in 1:cols(Q)) {\n",
    "        row_Q[i + (j-1)*rows(Q)] = Q[j,i];\n",
    "      }\n",
    "    irow_vec = (A_sharp_temp \\ add_diag(matrix_exp(A_sharp_temp * ts), -1)) * row_Q;\n",
    "    // irow takes elements of a column vector and puts them in a matrix rowwise\n",
    "    {\n",
    "      int row_size = rows(A);\n",
    "      int row_ticker = 1;\n",
    "      int col_ticker = 0;\n",
    "      for (i in 1:num_elements(irow_vec)) {\n",
    "        col_ticker += 1;\n",
    "        if (col_ticker > row_size) {\n",
    "          row_ticker += 1;\n",
    "          col_ticker = 1;\n",
    "        }\n",
    "        irow_mat[row_ticker,col_ticker] = irow_vec[i];\n",
    "      }\n",
    "    }\n",
    "    return irow_mat;\n",
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
    "  int N_tips; // number of tips\n",
    "  int N_obs; // number of observations\n",
    "  int J; // number of response traits\n",
    "  int N_seg; // total number of segments in the tree\n",
    "  array[N_seg] int node_seq; // index of tree nodes\n",
    "  array[N_seg] int parent; // index of the parent node of each descendent\n",
    "  array[N_seg] real ts; // time since parent\n",
    "  array[N_seg] int tip; // indicator of whether a given segment ends in a tip\n",
    "  array[J,J] int effects_mat; // which effects should be estimated?\n",
    "  int num_effects; // number of effects being estimated\n",
    "  matrix[N_obs,J] y; // observed data\n",
    "  matrix[N_obs,J] miss; // are data points missing?\n",
    "  array[N_obs] int tip_id; // index between 1 and N_tips that gives the group id\n"
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
      "  int prior_only; // should the likelihood be ignored?\n}"
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
    "  vector[num_effects - J] A_offdiag; // cross-lagged terms of A\n",
    "  vector<lower=0>[J] Q_diag; // self-drift terms\n",
    "  vector[J] b; // SDE intercepts\n",
    "  vector[J] eta_anc; // ancestral states\n",
    "  matrix[N_seg - 1,J] z_drift; // stochastic drift, unscaled and uncorrelated\n"
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
    "  matrix[N_seg,J] eta;\n",
    "  matrix[J,J] Q; // drift matrix\n",
    "  matrix[J,J] A; // selection matrix\n",
    "  vector[J*J - J] Q_offdiag = rep_vector(0.0, J*J - J);\n"
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
    "  matrix[N_seg,J] drift_tips; // terminal drift parameters\n",
    "  matrix[N_seg,J] sigma_tips; // terminal drift parameters\n",
    "  // fill A matrix //////////\n",
    "  {\n",
    "    int ticker = 1;\n",
    "    // fill upper triangle of matrix\n",
    "    for (i in 1:(J-1)) {\n",
    "      for (j in (i+1):J) {\n",
    "        if (effects_mat[i,j] == 1) {\n",
    "          A[i,j] = A_offdiag[ticker];\n",
    "          ticker += 1;\n",
    "        } else if (effects_mat[i,j] == 0) {\n",
    "          A[i,j] = 0;\n",
    "        }\n",
    "      }\n",
    "    }\n",
    "    // fill lower triangle of matrix\n",
    "    for (i in 1:(J-1)) {\n",
    "      for (j in (i+1):J) {\n",
    "        if (effects_mat[j,i] == 1) {\n",
    "          A[j,i] = A_offdiag[ticker];\n",
    "          ticker += 1;\n",
    "        } else if (effects_mat[j,i] == 0) {\n",
    "          A[j,i] = 0;\n",
    "        }\n",
    "      }\n",
    "    }\n",
    "    // fill diag of matrix\n",
    "    for (j in 1:J) A[j,j] = A_diag[j];\n",
    "  }\n",
    "  // fill Q matrix //////////\n",
    "  {\n",
    "    int ticker = 1;\n",
    "    for (i in 1:(J-1))\n",
    "      for (j in (i+1):J) {\n",
    "        Q[i,j] = Q_offdiag[ticker];\n",
    "        Q[j,i] = Q[i,j]; // symmetry of covariance\n",
    "        ticker += 1;\n",
    "      }\n",
    "    for (j in 1:J) Q[j,j] = Q_diag[j];\n",
    "  }\n",
    "  // setting ancestral states and placeholders\n",
    "  for (j in 1:J) {\n",
    "    eta[node_seq[1],j] = eta_anc[j];\n",
    "    drift_tips[node_seq[1],j] = -99;\n",
    "    sigma_tips[node_seq[1],j] = -99;\n",
    "  }\n",
    "  for (i in 2:N_seg) {\n",
    "    matrix[J,J] A_delta; // amount of deterministic change (selection)\n",
    "    matrix[J,J] VCV; // variance-covariance matrix of stochastic change (drift)\n",
    "    vector[J] drift_seg; // accumulated drift over the segment\n",
    "    A_delta = matrix_exp(A * ts[i]);\n",
    "    VCV = cov_drift(A, Q, ts[i]);\n",
    "    drift_seg = cholesky_decompose(VCV) * to_vector( z_drift[i-1,] );\n",
    "    // if not a tip, add the drift parameter\n",
    "    if (tip[i] == 0) {\n",
    "      eta[node_seq[i],] = to_row_vector(\n",
    "        A_delta * to_vector(eta[parent[i],]) + ((A \\ add_diag(A_delta, -1)) * b) + drift_seg\n",
    "      );\n",
    "      drift_tips[node_seq[i],] = to_row_vector(rep_vector(-99, J));\n",
    "      sigma_tips[node_seq[i],] = to_row_vector(rep_vector(-99, J));\n",
    "    }\n",
    "    // if is a tip, omit, we'll deal with it in the model block;\n",
    "    else {\n",
    "      eta[node_seq[i],] = to_row_vector(\n",
    "        A_delta * to_vector(eta[parent[i],]) + ((A \\ add_diag(A_delta, -1)) * b)\n",
    "      );\n",
    "      drift_tips[node_seq[i],] = to_row_vector(drift_seg);\n",
    "      sigma_tips[node_seq[i],] = to_row_vector(diagonal(Q));\n",
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
    "  eta_anc ~ ", priors$eta_anc, ";\n",
    "  to_vector(z_drift) ~ std_normal();\n",
    "  A_offdiag ~ ", priors$A_offdiag, ";\n",
    "  A_diag ~ ", priors$A_diag, ";\n",
    "  Q_diag ~ ", priors$Q_diag, ";\n"
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
    "    for (i in 1:N_obs) {\n"
    )
  for (j in 1:length(distributions)) {
    if (distributions[j] == "bernoulli_logit") {
      sc_model <- paste0(
        sc_model,
        "        if (miss[i,", j, "] == 0) to_int(y[i,", j, "]) ~ ",
        "bernoulli_logit(eta[tip_id[i],", j, "]",
        ifelse(!is.null(dist_mat), paste0(" + dist_v[tip_id[i],", j, "]"), ""),
        ifelse(any(duplicated(data[,id])), paste0(" + group_v[tip_id[i],", j, "]"), ""),
        " + drift_tips[tip_id[i],", j, "]);\n"
        )
    } else if (distributions[j] == "ordered_logistic") {
      sc_model <- paste0(
        sc_model,
        "        if (miss[i,", j, "] == 0) to_int(y[i,", j, "]) ~ ",
        "ordered_logistic(eta[tip_id[i],", j, "]",
        ifelse(!is.null(dist_mat), paste0(" + dist_v[tip_id[i],", j, "]"), ""),
        ifelse(any(duplicated(data[,id])), paste0(" + group_v[tip_id[i],", j, "]"), ""),
        " + drift_tips[tip_id[i],", j, "], c", j, ");\n"
        )
    } else if (distributions[j] == "poisson_softplus") {
      sc_model <- paste0(
        sc_model,
        "        if (miss[i,", j, "] == 0) to_int(y[i,", j, "]) ~ ",
        "poisson(mean(obs", j, ") * log1p_exp(eta[tip_id[i],", j, "]",
        ifelse(!is.null(dist_mat), paste0(" + dist_v[tip_id[i],", j, "]"), ""),
        ifelse(any(duplicated(data[,id])), paste0(" + group_v[tip_id[i],", j, "]"), ""),
        " + drift_tips[tip_id[i],", j, "]));\n"
        )
    } else if (distributions[j] == "normal") {
      sc_model <- paste0(
        sc_model,
        "        if (miss[i,", j, "] == 0) y[i,", j, "] ~ ",
        "normal(eta[tip_id[i],", j, "]",
        ifelse(!is.null(dist_mat), paste0(" + dist_v[tip_id[i],", j, "]"), ""),
        ifelse(any(duplicated(data[,id])), paste0(" + group_v[tip_id[i],", j, "]"), ""),
        ", sigma_tips[tip_id[i],", j, "]);\n"
        )
    } else if (distributions[j] == "student_t") {
      sc_model <- paste0(
        sc_model,
        "        if (miss[i,", j, "] == 0) y[i,", j, "] ~ ",
        "student_t(nu", j, ", eta[tip_id[i],", j, "]",
        ifelse(!is.null(dist_mat), paste0(" + dist_v[tip_id[i],", j, "]"), ""),
        ifelse(any(duplicated(data[,id])), paste0(" + group_v[tip_id[i],", j, "]"), ""),
        ", sigma_tips[tip_id[i],", j, "]);\n"
      )
    } else if (distributions[j] == "lognormal") {
      sc_model <- paste0(
        sc_model,
        "        if (miss[i,", j, "] == 0) y[i,", j, "] ~ ",
        "lognormal(eta[tip_id[i],", j, "]",
        ifelse(!is.null(dist_mat), paste0(" + dist_v[tip_id[i],", j, "]"), ""),
        ifelse(any(duplicated(data[,id])), paste0(" + group_v[tip_id[i],", j, "]"), ""),
        ", sigma_tips[tip_id[i],", j, "]);\n"
        )
    } else if (distributions[j] == "negative_binomial_softplus") {
      sc_model <- paste0(
        sc_model,
        "        if (miss[i,", j, "] == 0) to_int(y[i,", j, "]) ~ ",
        "neg_binomial_2(mean(obs", j, ") * log1p_exp(eta[tip_id[i],", j, "]",
        ifelse(!is.null(dist_mat), paste0(" + dist_v[tip_id[i],", j, "]"), ""),
        ifelse(any(duplicated(data[,id])), paste0(" + group_v[tip_id[i],", j, "]"), ""),
        " + drift_tips[tip_id[i],", j, "]), phi", j, ");\n"
        )
    }
  }
  sc_model <- paste0(
    sc_model,
    "    }\n",
    "  }\n",
    "}"
    )
  # generated quantities block
  sc_generated_quantities <-
    paste0(
      "generated quantities{\n",
      "  vector[N_obs*J] log_lik; // log-likelihood\n",
      "  matrix[N_obs,J] yrep; // predictive checks\n"
      )
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
      "    matrix[N_obs,J] yrep_temp;\n",
      "    for (i in 1:N_obs) {\n"
      )
  for (j in 1:length(distributions)) {
    if (distributions[j] == "bernoulli_logit") {
      sc_generated_quantities <- paste0(
        sc_generated_quantities,
        "      if (miss[i,", j, "] == 0) log_lik_temp[i,", j, "] = ",
        "bernoulli_logit_lpmf(to_int(y[i,", j, "]) | eta[tip_id[i],", j, "]",
        ifelse(!is.null(dist_mat), paste0(" + dist_v[tip_id[i],", j, "]"), ""),
        ifelse(any(duplicated(data[,id])), paste0(" + group_v[tip_id[i],", j, "]"), ""),
        " + drift_tips[tip_id[i],", j, "]);\n",
        "      yrep_temp[i,", j, "] = ",
        "bernoulli_logit_rng(eta[tip_id[i],", j, "]",
        ifelse(!is.null(dist_mat), paste0(" + dist_v[tip_id[i],", j, "]"), ""),
        ifelse(any(duplicated(data[,id])), paste0(" + group_v[tip_id[i],", j, "]"), ""),
        " + drift_tips[tip_id[i],", j, "]);\n"
      )
    } else if (distributions[j] == "ordered_logistic") {
      sc_generated_quantities <- paste0(
        sc_generated_quantities,
        "      if (miss[i,", j, "] == 0) log_lik_temp[i,", j, "] = ",
        "ordered_logistic_lpmf(to_int(y[i,", j, "]) | eta[tip_id[i],", j, "]",
        ifelse(!is.null(dist_mat), paste0(" + dist_v[tip_id[i],", j, "]"), ""),
        ifelse(any(duplicated(data[,id])), paste0(" + group_v[tip_id[i],", j, "]"), ""),
        " + drift_tips[tip_id[i],", j, "], c", j, ");\n",
        "      yrep_temp[i,", j, "] = ",
        "ordered_logistic_rng(eta[tip_id[i],", j, "]",
        ifelse(!is.null(dist_mat), paste0(" + dist_v[tip_id[i],", j, "]"), ""),
        ifelse(any(duplicated(data[,id])), paste0(" + group_v[tip_id[i],", j, "]"), ""),
        " + drift_tips[tip_id[i],", j, "], c", j, ");\n"
      )
    } else if (distributions[j] == "poisson_softplus") {
      sc_generated_quantities <- paste0(
        sc_generated_quantities,
        "      if (miss[i,", j, "] == 0) log_lik_temp[i,", j, "] = ",
        "poisson_lpmf(to_int(y[i,", j, "]) | mean(obs", j,
        ") * log1p_exp(eta[tip_id[i],", j, "]",
        ifelse(!is.null(dist_mat), paste0(" + dist_v[tip_id[i],", j, "]"), ""),
        ifelse(any(duplicated(data[,id])), paste0(" + group_v[tip_id[i],", j, "]"), ""),
        " + drift_tips[tip_id[i],", j, "]));\n",
        "      yrep_temp[i,", j, "] = ",
        "poisson_rng(mean(obs", j, ") * log1p_exp(eta[tip_id[i],", j, "]",
        ifelse(!is.null(dist_mat), paste0(" + dist_v[tip_id[i],", j, "]"), ""),
        ifelse(any(duplicated(data[,id])), paste0(" + group_v[tip_id[i],", j, "]"), ""),
        " + drift_tips[tip_id[i],", j, "]));\n"
      )
    } else if (distributions[j] == "normal") {
      sc_generated_quantities <- paste0(
        sc_generated_quantities,
        "      if (miss[i,", j, "] == 0) log_lik_temp[i,", j, "] = ",
        "normal_lpdf(y[i,", j, "] | eta[tip_id[i],", j, "]",
        ifelse(!is.null(dist_mat), paste0(" + dist_v[tip_id[i],", j, "]"), ""),
        ifelse(any(duplicated(data[,id])), paste0(" + group_v[tip_id[i],", j, "]"), ""),
        ", sigma_tips[tip_id[i],", j, "]);\n",
        "      yrep_temp[i,", j, "] = ",
        "normal_rng(eta[tip_id[i],", j, "]",
        ifelse(!is.null(dist_mat), paste0(" + dist_v[tip_id[i],", j, "]"), ""),
        ifelse(any(duplicated(data[,id])), paste0(" + group_v[tip_id[i],", j, "]"), ""),
        ", sigma_tips[tip_id[i],", j, "]);\n"
      )
    } else if (distributions[j] == "student_t") {
      sc_generated_quantities <- paste0(
        sc_generated_quantities,
        "      if (miss[i,", j, "] == 0) log_lik_temp[i,", j, "] = ",
        "student_t_lpdf(y[i,", j, "] | nu", j, ", eta[tip_id[i],", j, "]",
        ifelse(!is.null(dist_mat), paste0(" + dist_v[tip_id[i],", j, "]"), ""),
        ifelse(any(duplicated(data[,id])), paste0(" + group_v[tip_id[i],", j, "]"), ""),
        ", sigma_tips[tip_id[i],", j, "]);\n",
        "      yrep_temp[i,", j, "] = ",
        "student_t_rng(nu", j, ", eta[tip_id[i],", j, "]",
        ifelse(!is.null(dist_mat), paste0(" + dist_v[tip_id[i],", j, "]"), ""),
        ifelse(any(duplicated(data[,id])), paste0(" + group_v[tip_id[i],", j, "]"), ""),
        ", sigma_tips[tip_id[i],", j, "]);\n"
      )
    } else if (distributions[j] == "lognormal") {
      sc_generated_quantities <- paste0(
        sc_generated_quantities,
        "      if (miss[i,", j, "] == 0) log_lik_temp[i,", j, "] = ",
        "lognormal_lpdf(y[i,", j, "] | eta[tip_id[i],", j, "]",
        ifelse(!is.null(dist_mat), paste0(" + dist_v[tip_id[i],", j, "]"), ""),
        ifelse(any(duplicated(data[,id])), paste0(" + group_v[tip_id[i],", j, "]"), ""),
        ", sigma_tips[tip_id[i],", j, "]);\n",
        "      yrep_temp[i,", j, "] = ",
        "lognormal_rng(eta[tip_id[i],", j, "]",
        ifelse(!is.null(dist_mat), paste0(" + dist_v[tip_id[i],", j, "]"), ""),
        ifelse(any(duplicated(data[,id])), paste0(" + group_v[tip_id[i],", j, "]"), ""),
        ", sigma_tips[tip_id[i],", j, "]);\n"
      )
    } else if (distributions[j] == "negative_binomial_softplus") {
      sc_generated_quantities <- paste0(
        sc_generated_quantities,
        "      if (miss[i,", j, "] == 0) log_lik_temp[i,", j, "] = ",
        "neg_binomial_2_lpmf(to_int(y[i,", j, "]) | mean(obs", j,
        ") * log1p_exp(eta[tip_id[i],", j, "]",
        ifelse(!is.null(dist_mat), paste0(" + dist_v[tip_id[i],", j, "]"), ""),
        ifelse(any(duplicated(data[,id])), paste0(" + group_v[tip_id[i],", j, "]"), ""),
        " + drift_tips[tip_id[i],", j, "]), phi", j, ");\n",
        "      yrep_temp[i,", j, "] = ",
        "neg_binomial_2_rng(mean(obs", j, ") * log1p_exp(eta[tip_id[i],", j, "]",
        ifelse(!is.null(dist_mat), paste0(" + dist_v[tip_id[i],", j, "]"), ""),
        ifelse(any(duplicated(data[,id])), paste0(" + group_v[tip_id[i],", j, "]"), ""),
        " + drift_tips[tip_id[i],", j, "]), phi", j, ");\n"
      )
    }
  }
  sc_generated_quantities <-
    paste0(
      sc_generated_quantities,
      "    }\n",
      "  yrep = yrep_temp;\n",
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
  # return stan code
  return(sc)
}
