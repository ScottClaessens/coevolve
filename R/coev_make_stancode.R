#' Make Stan code for dynamic coevolutionary model
#'
#' @param data An object of class \code{data.frame} (or one that can be coerced
#'   to that class) containing data of all variables used in the model.
#' @param variables A named list identifying variables that should coevolve in
#'   the model and their associated response distributions as character strings (e.g.
#'   \code{list(var1 = "bernoulli_logit", var2 = "ordered_logistic")}). Must identify
#'   at least two variables. Variable names must refer to valid column names in data.
#'   Currently, the only supported response distributions are \code{bernoulli_logit}
#'   and \code{ordered_logistic}.
#' @param id A character of length one identifying the variable in the data that links rows to tips
#'   on the phylogeny. Must refer to a valid column name in the data. The id column
#'   must exactly match the tip labels in the phylogeny.
#' @param tree A phylogenetic tree object of class \code{phylo}.
#' @param prior A list of priors for the model.
#'
#' @return A character string containing the \pkg{Stan} code to fit the dynamic coevolutionary model.
#' @export
#'
#' @examples
#' # simulate data
#' set.seed(1)
#' n <- 20
#' tree <- ape::rtree(n)
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
coev_make_stancode <- function(data, variables, id, tree, prior = NULL) {
  # check arguments
  run_checks(data, variables, id, tree)
  # extract distributions and variable names from named list
  distributions <- as.character(variables)
  variables <- names(variables)
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
    "  array[N,J] int y; // observed data\n",
    "}"
  )
  # write parameters block
  sc_parameters <- paste0(
    "parameters{\n",
    "  matrix[J,J] alpha; // selection coefficients\n",
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
  sc_parameters <- paste0(sc_parameters, "}")
  # write transformed parameters block
  sc_transformed_parameters <- paste0(
    "transformed parameters{\n",
    "  matrix[N_seg,J] eta;\n",
    "  matrix[J,J] Q;\n",
    "  matrix[J,J] I;\n",
    "  matrix[J,J] A;\n",
    "  // selection matrix\n",
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
    "  // setting ancestral states\n",
    "  for (j in 1:J) {\n",
    "    eta[node_seq[1],j] = eta_anc[j];\n",
    "  }\n",
    "  for (i in 2:N_seg) {\n",
    "    matrix[J,J] A_delta; // amount of deterministic change (selection)\n",
    "    matrix[J,J] VCV; // variance-covariance matrix of stochastic change (drift)\n",
    "    vector[J] drift_seg; // accumulated drift over the segment\n",
    "    A_delta = A_dt(A, ts[i]);\n",
    "    VCV = cov_drift(A, Q, ts[i]);\n",
    "    // No drift on the interaction, bc its simply a product of vars\n",
    "    drift_seg = cholesky_decompose(VCV) * to_vector( z_drift[i-1,] );\n",
    "    eta[node_seq[i],] = to_row_vector(\n",
    "      A_delta * to_vector(eta[parent[i],]) + (inverse(A) * (A_delta - I) * b) + drift_seg\n",
    "    );\n",
    "  }\n",
    "}"
  )
  # write model block
  sc_model <- paste0(
    "model{\n",
    "  to_vector(alpha) ~ std_normal();\n",
    "  b ~ std_normal();\n",
    "  sigma ~ std_normal();\n",
    "  eta_anc ~ std_normal();\n",
    "  to_vector(z_drift) ~ std_normal();\n"
  )
  # add priors for any cut points
  for (j in 1:length(distributions)) {
    if (distributions[j] == "ordered_logistic") {
      sc_model <- paste0(sc_model, "  c", j, " ~ normal(0, 2);\n")
    }
  }
  # add response distributions
  sc_model <- paste0(sc_model, "  for (i in 1:N) {\n")
  for (j in 1:length(distributions)) {
    if (distributions[j] == "bernoulli_logit") {
      sc_model <- paste0(sc_model, "      y[i,", j, "] ~ bernoulli_logit(eta[i,", j, "]);\n")
    } else if (distributions[j] == "ordered_logistic") {
      sc_model <- paste0(sc_model, "      y[i,", j, "] ~ ordered_logistic(eta[i,", j, "], c", j, ");\n")
    }
  }
  sc_model <- paste0(sc_model, "  }\n}")
  # put stan code together
  sc <- paste0(
    "// Generated by the coevolve package\n",
    sc_functions, "\n",
    sc_data, "\n",
    sc_parameters, "\n",
    sc_transformed_parameters, "\n",
    sc_model
  )
  return(sc)
}
