#' Make Stan code for dynamic coevolutionary model
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
#' @param effects_mat (optional) A boolean matrix with row and column names
#'   exactly matching the variables declared for the model. If not specified,
#'   all cross-lagged effects will be estimated in the model. If specified, the
#'   model will only estimate cross-lagged effects where cells in the matrix are
#'   TRUE and will ignore cross-lagged effects where cells in the matrix are FALSE.
#'   In the matrix, columns represent predictor variables and rows represent
#'   outcome variables. All autoregressive effects (e.g., X -> X) must be TRUE
#'   in the matrix.
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
#'
#' @return A character string containing the \pkg{Stan} code to fit the dynamic coevolutionary model.
#' @export
#'
#' @examples
#' # simulate data
#' n <- 20
#' tree <- ape::rcoal(n)
#' d <- data.frame(
#'    id = tree$tip.label,
#'    x = rbinom(n, size = 1, prob = 0.5),
#'    y = ordered(sample(1:4, size = n, replace = TRUE))
#' )
#' # make stan code
#' coev_make_stancode(
#'    data = d,
#'    variables = list(
#'        x = "bernoulli_logit",
#'        y = "ordered_logistic"
#'    ),
#'    id = "id",
#'    tree = tree
#' )
coev_make_stancode <- function(data, variables, id, tree,
                               effects_mat = NULL, dist_mat = NULL,
                               prior = NULL, prior_only = FALSE) {
  # check arguments
  run_checks(data, variables, id, tree, effects_mat, dist_mat, prior, prior_only)
  # extract distributions and variable names from named list
  distributions <- as.character(variables)
  variables <- names(variables)
  # get priors
  priors <-
    list(
      alpha      = "std_normal()",
      b          = "std_normal()",
      sigma      = "std_normal()",
      eta_anc    = "std_normal()",
      c          = "normal(0, 2)",
      sigma_dist = "exponential(1)",
      rho_dist   = "exponential(1)"
    )
  # replace priors if user has explicitly set them
  if (!is.null(prior)) {
    for (i in names(prior)) {
      priors[[i]] <- prior[[i]]
    }
  }
  # write functions block
  sc_functions <- paste0(
    "functions {\n",
    "  // returns the Kronecker Product\n",
    "  matrix kronecker_prod(matrix A, matrix B) {\n",
    "    matrix[rows(A) * rows(B), cols(A) * cols(B)] C;\n",
    "    int m;\n",
    "    int n;\n",
    "    int p;\n",
    "    int q;\n",
    "    m = rows(A);\n",
    "    n = cols(A);\n",
    "    p = rows(B);\n",
    "    q = cols(B);\n",
    "    for (i in 1:m)\n",
    "      for (j in 1:n)\n",
    "        for (k in 1:p)\n",
    "          for (l in 1:q)\n",
    "            C[p*(i-1)+k,q*(j-1)+l] = A[i,j]*B[k,l];\n",
    "    return C;\n",
    "  }\n",
    "  \n",
    "  // expected auto and cross effects over a discrete time t\n",
    "  matrix A_dt(matrix A, real t) {\n",
    "    return( matrix_exp(A * t) );\n",
    "  }\n",
    "  \n",
    "  // calculate A sharp matrix\n",
    "  matrix A_sharp(matrix A) {\n",
    "    matrix[rows(A) * rows(A), cols(A) * cols(A)] A_temp;\n",
    "    matrix[rows(A),cols(A)] I; // identity matrix\n",
    "    I = diag_matrix(rep_vector(1,rows(A)));\n",
    "    A_temp = kronecker_prod(A,I) + kronecker_prod(I,A);\n",
    "    return(A_temp);\n",
    "  }\n",
    "  \n",
    "  // solve SDE\n",
    "  matrix cov_drift(matrix A, matrix Q, real ts) {\n",
    "    matrix[rows(A) * rows(A), cols(A) * cols(A)] A_sharp_temp;\n",
    "    matrix[rows(A) * rows(A), cols(A) * cols(A)] I; // identity matrix\n",
    "    vector[rows(Q)*cols(Q)] row_Q;\n",
    "    vector[rows(A)*cols(A)] irow_vec;\n",
    "    matrix[rows(A),cols(A)] irow_mat;\n",
    "    I = diag_matrix(rep_vector(1,rows(A_sharp_temp)));\n",
    "    A_sharp_temp = A_sharp(A);\n",
    "    // row operation takes elements of a matrix rowwise and puts them into a column vector\n",
    "    for (i in 1:rows(Q))\n",
    "      for (j in 1:cols(Q)) {\n",
    "        row_Q[i + (j-1)*rows(Q)] = Q[j,i];\n",
    "      }\n",
    "    irow_vec = inverse(A_sharp_temp) * (matrix_exp(A_sharp_temp * ts) - I) * row_Q;\n",
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
    "    return(irow_mat);\n",
    "  }\n",
    "}"
  )
  # write data block
  sc_data <- paste0(
    "data{\n",
    "  int N; // number of taxa\n",
    "  int J; // number of variables\n",
    "  int N_seg; // number of segments in tree\n",
    "  array[N_seg] int node_seq; // sequence of nodes\n",
    "  array[N_seg] int parent; // parent node for each node\n",
    "  array[N_seg] real ts; // amount of time since parent node\n",
    "  array[N_seg] int tip; // is tip?\n",
    "  array[J,J] int effects_mat; // which effects should be estimated?\n",
    "  int num_effects; // number of effects being estimated\n"
  )
  # add observed data variables one by one depending on type
  for (i in 1:length(variables)) {
    sc_data <- paste0(
      sc_data,
      "  array[N] ",
      # continuous variables are real numbers
      # all others are integers
      ifelse(distributions[i] %in% c("normal", "lognormal"), "real", "int"),
      " y", i, "; // observed data variable ", i,
      "\n")
  }
  # add distance matrix if user has defined one
  if (!is.null(dist_mat)) {
    sc_data <- paste0(sc_data, "  array[N,N] real dist_mat; // distance matrix\n")
  }
  # add prior_only data variable
  sc_data <- paste0(sc_data, "  int prior_only; // should the likelihood be ignored?\n}")
  # write parameters block
  sc_parameters <- paste0(
    "parameters{\n",
    "  vector[num_effects] alpha_pars; // selection coefficients (actual parameters)\n",
    "  vector<lower=0>[J] sigma; // drift scale\n",
    "  vector[J] b; // SDE intercepts\n",
    "  vector[J] eta_anc; // ancestral states\n",
    "  matrix[N_seg - 1,J] z_drift; // stochastic drift, unscaled and uncorrelated\n"
  )
  for (i in 1:length(distributions)) {
    # add cut points for ordinal_logistic distributions
    if (distributions[i] == "ordered_logistic") {
      # calculate number of cut points (number of levels - 1)
      num_cuts <- max(as.numeric(data[,variables[i]])) - 1
      sc_parameters <- paste0(sc_parameters, "  ordered[", num_cuts, "] c", i, "; // cut points for variable ", i, "\n")
    }
  }
  # add gaussian process parameters if distance matrix specified by user
  if (!is.null(dist_mat)) {
    sc_parameters <- paste0(
      sc_parameters,
      "  matrix[N,J] dist_z; // spatial covariance random effects, unscaled and uncorrelated\n",
      "  vector<lower=0>[J] rho_dist; // how quickly does covariance decline with distance\n",
      "  vector<lower=0>[J] sigma_dist; // maximum covariance due to spatial location\n"
    )
  }
  sc_parameters <- paste0(sc_parameters, "}")
  # write transformed parameters block
  sc_transformed_parameters <- paste0(
    "transformed parameters{\n",
    "  matrix[N_seg,J] eta;\n",
    "  matrix[J,J] Q;\n",
    "  matrix[J,J] I;\n",
    "  matrix[J,J] A;\n",
    "  matrix[J,J] alpha;\n"
  )
  # add distance random effects if distance matrix specified by user
  if (!is.null(dist_mat)) {
    sc_transformed_parameters <- paste0(
      sc_transformed_parameters,
      "  matrix[N,J] dist_v; // distance covariance random effects\n"
    )
  }
  sc_transformed_parameters <- paste0(
    sc_transformed_parameters,
    "  matrix[N_seg,J] drift_tips; // terminal drift parameters, saved here to use in likelihood for Gaussian outcomes\n",
    "  matrix[N_seg,J] sigma_tips; // terminal drift parameters, saved here to use in likelihood for Gaussian outcomes\n",
    "  // construct alpha matrix\n",
    "  {\n",
    "    int index;\n",
    "    index = 1;\n",
    "    for (i in 1:J) {\n",
    "      for (j in 1:J) {\n",
    "        if (effects_mat[i,j] == 1) {\n",
    "          alpha[i,j] = alpha_pars[index];\n",
    "          index += 1;\n",
    "        } else if (effects_mat[i,j] == 0) {\n",
    "          alpha[i,j] = 0;\n",
    "        }\n",
    "      }\n",
    "    }\n",
    "  }\n",
    "  // calculate selection matrix\n",
    "  for (j in 1:J) {\n",
    "    for (i in 1:J) {\n",
    "      if (i == j) {\n",
    "        A[i,j] = -exp(alpha[i,j]); // autoregressive effects\n",
    "      } else {\n",
    "        A[i,j] = alpha[i,j]; // cross effects\n",
    "      }\n",
    "    }\n",
    "  }\n",
    "  // drift matrix\n",
    "  Q = diag_matrix(square(sigma));\n",
    "  // identity matrix\n",
    "  I = diag_matrix(rep_vector(1,J));\n",
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
    "    A_delta = A_dt(A, ts[i]);\n",
    "    VCV = cov_drift(A, Q, ts[i]);\n",
    "    // no drift on the interaction, bc its simply a product of vars\n",
    "    drift_seg = cholesky_decompose(VCV) * to_vector( z_drift[i-1,] );\n",
    "    // if not a tip, add the drift parameter\n",
    "    if (tip[i] == 0) {\n",
    "      eta[node_seq[i],] = to_row_vector(\n",
    "        A_delta * to_vector(eta[parent[i],]) + (inverse(A) * (A_delta - I) * b) + drift_seg\n",
    "      );\n",
    "      drift_tips[node_seq[i],] = to_row_vector(rep_vector(-99, J));\n",
    "      sigma_tips[node_seq[i],] = to_row_vector(rep_vector(-99, J));\n",
    "    }\n",
    "    // if is a tip, omit, we'll deal with it in the model block;\n",
    "    else {\n",
    "      eta[node_seq[i],] = to_row_vector(\n",
    "        A_delta * to_vector(eta[parent[i],]) + (inverse(A) * (A_delta - I) * b)\n",
    "      );\n",
    "      drift_tips[node_seq[i],] = to_row_vector(drift_seg);\n",
    "      sigma_tips[node_seq[i],] = to_row_vector(sqrt(diagonal(Q)));\n",
    "    }\n",
    "  }\n"
  )
  # add gaussian process functions if dist_mat specified by user
  if (!is.null(dist_mat)) {
    sc_transformed_parameters <- paste0(
      sc_transformed_parameters,
      "  // distance covariance functions\n",
      "  for (j in 1:J) {\n",
      "    matrix[N,N] dist_cov;\n",
      "    matrix[N,N] L_dist_cov;\n",
      "    for ( i in 1:(N-1) )\n",
      "      for ( m in (i+1):N ) {\n",
      "        dist_cov[i,m] = sigma_dist[j]*exp(-(rho_dist[j]*dist_mat[i,m]));\n",
      "        dist_cov[m,i] = dist_cov[i,m];\n",
      "      }\n",
      "    for ( q in 1:N )\n",
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
    "  to_vector(alpha_pars) ~ ", priors$alpha, ";\n",
    "  b ~ ", priors$b, ";\n",
    "  sigma ~ ", priors$sigma, ";\n",
    "  eta_anc ~ ", priors$eta_anc, ";\n",
    "  to_vector(z_drift) ~ std_normal();\n"
  )
  # add priors for any cut points
  for (j in 1:length(distributions)) {
    if (distributions[j] == "ordered_logistic") {
      sc_model <- paste0(sc_model, "  c", j, " ~ ", priors$c, ";\n")
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
  # add likelihood
  sc_model <- paste0(
    sc_model,
    "  if (!prior_only) {\n",
    "    for (i in 1:N) {\n"
    )
  for (j in 1:length(distributions)) {
    if (distributions[j] == "bernoulli_logit") {
      sc_model <- paste0(
        sc_model,
        "        y", j, "[i] ~ ",
        "bernoulli_logit(eta[i,", j, "]",
        ifelse(!is.null(dist_mat), paste0(" + dist_v[i,", j, "]"), ""),
        " + drift_tips[i,", j, "]);\n")
    } else if (distributions[j] == "ordered_logistic") {
      sc_model <- paste0(
        sc_model,
        "        y", j, "[i] ~ ",
        "ordered_logistic(eta[i,", j, "]",
        ifelse(!is.null(dist_mat), paste0(" + dist_v[i,", j, "]"), ""),
        " + drift_tips[i,", j, "], c", j, ");\n")
    } else if (distributions[j] == "poisson_softmax") {
      sc_model <- paste0(
        sc_model,
        "        y", j, "[i] ~ ",
        "poisson(log1p_exp(eta[i,", j, "]",
        ifelse(!is.null(dist_mat), paste0(" + dist_v[i,", j, "]"), ""),
        " + drift_tips[i,", j, "]));\n")
    } else if (distributions[j] == "normal") {
      sc_model <- paste0(
        sc_model,
        "        y", j, "[i] ~ ",
        "normal(eta[i,", j, "]",
        ifelse(!is.null(dist_mat), paste0(" + dist_v[i,", j, "]"), ""),
        ", sigma_tips[i,", j, "]);\n")
    } else if (distributions[j] == "lognormal") {
      sc_model <- paste0(
        sc_model,
        "        y", j, "[i] ~ ",
        "lognormal(eta[i,", j, "]",
        ifelse(!is.null(dist_mat), paste0(" + dist_v[i,", j, "]"), ""),
        ", sigma_tips[i,", j, "]);\n")
    }
  }
  sc_model <- paste0(
    sc_model,
    "    }\n",
    "  }\n",
    "}"
    )
  # put stan code together
  sc <- paste0(
    "// Generated by the coevolve package\n",
    sc_functions, "\n",
    sc_data, "\n",
    sc_parameters, "\n",
    sc_transformed_parameters, "\n",
    sc_model
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
