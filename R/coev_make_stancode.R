coev_make_stancode <- function(data, variables, distributions, tree, prior = NULL) {
  ## 0. Check arguments

  # stop if variables are not valid column names in data
  if (!all(variables %in% colnames(data))) {
    stop("Some variable names are not valid column names in data.")
  }
  # stop if response distributions are not valid
  if (!all(distributions %in% c("bernoulli_logit", "ordered_logistic"))) {
    stop("Response distributions other than 'bernoulli_logit' and 'ordered_logistic' are not yet supported.")
  }
  # stop if not at least two variables
  if (length(variables) < 2) {
    stop("Must be at least two variables.")
  }
  # stop if number of variables and response distributions do not match
  if (length(variables) != length(distributions)) {
    stop("Must be the same number of variables and response distributions.")
  }
  # stop if any binary variables are not 0/1 integers
  for (i in 1:length(distributions)) {
    if (distributions[i] == "bernoulli_logit" & (!is.integer(data[,variables[i]]) | !all(data[,variables[i]] %in% 0:1))) {
      stop("Variables following the 'bernoulli_logit' response distribution must be integers with values of 0/1.")
    }
  }
  # stop if any ordinal variables are not ordered factors in data
  for (i in 1:length(distributions)) {
    if (distributions[i] == "ordered_logistic" & !is.ordered(data[,variables[i]])) {
      stop("Variables following the 'ordered_logistic' response distribution must be ordered factors.")
    }
  }

  ## 1. Write Stan code for model

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
    "  int node_seq[N_seg]; // sequence of nodes\n",
    "  int parent[N_seg]; // parent node for each node\n",
    "  vector[N_seg] ts; // amount of time since parent node\n",
    "  int y[N,J]; // observed data\n",
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
