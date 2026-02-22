transformed parameters{

  array[N_tree, N_seg] vector[J] eta;
  matrix[J,J] A = diag_matrix(A_diag); // selection matrix
  {{#estimate_correlated_drift}}
  matrix[J,J] Q = diag_matrix(Q_sigma) * (L_R * L_R') * diag_matrix(Q_sigma); // drift matrix
  {{/estimate_correlated_drift}}
  {{#no_correlated_drift}}
  matrix[J,J] Q = diag_matrix(Q_sigma^2); // drift matrix
  {{/no_correlated_drift}}
  matrix[J,J] Q_inf; // asymptotic covariance matrix
  array[N_tree, N_seg] matrix[J,J] VCV_tips; // vcov matrix for drift
  array[N_tree, N_seg] matrix[J,J] L_VCV_tips; // Cholesky factor of VCV_tips
  {{#dist_mat}}
  matrix[N_tips,J] dist_v; // distance covariance random effects
  {{/dist_mat}}
  {{#tdrift}}
  array[N_tree,N_tips] vector[J] tdrift; // terminal drift
  {{/tdrift}}
  {{#residual}}
  matrix[N_obs,J] residual_v; // residual pars
  // scale and correlate residual pars
  residual_v = (diag_pre_multiply(sigma_residual, L_residual) * residual_z)';
  {{/residual}}
  // fill off diagonal of A matrix
  {
    int ticker = 1;
    for (i in 1:J) {
      for (j in 1:J) {
        if (i != j) {
          if (effects_mat[i,j] == 1) {
            A[i,j] = A_offdiag[ticker];
            ticker += 1;
          } else if (effects_mat[i,j] == 0) {
            A[i,j] = 0;
          }
        }
      }
    }
  }
  // calculate asymptotic covariance
  Q_inf = ksolve(A, Q);
  {
    array[N_unique_lengths] matrix[J,J] A_delta_cache;
    array[N_unique_lengths] matrix[J,J] VCV_cache;
    array[N_unique_lengths] matrix[J,J] L_VCV_cache;
    array[N_unique_lengths] matrix[J,J] A_solve_cache;
    for (u in 1:N_unique_lengths) {
      A_delta_cache[u] = matrix_exp(A * unique_lengths[u]);
      VCV_cache[u] = Q_inf - quad_form_sym(Q_inf, A_delta_cache[u]');
      L_VCV_cache[u] = cholesky_decompose(VCV_cache[u]);
      A_solve_cache[u] = A \ add_diag(A_delta_cache[u], -1);
      for (i in 1:J) {
        for (j in 1:i) {
          real val = 0.5 * (A_solve_cache[u][i, j] + A_solve_cache[u][j, i]);
          A_solve_cache[u][i, j] = val;
          A_solve_cache[u][j, i] = val;
        }
      }
    }
    for (t in 1:N_tree) {
      // setting ancestral states and placeholders
      eta[t, node_seq[t, 1]] = eta_anc[t];
      VCV_tips[t, node_seq[t, 1]] = diag_matrix(rep_vector(-99, J));
      L_VCV_tips[t, node_seq[t, 1]] = diag_matrix(rep_vector(1.0, J));
      for (i in 2:N_seg) {
        matrix[J,J] A_delta;
        matrix[J,J] VCV;
        vector[J] drift_seg;
        matrix[J,J] L_VCV;
        matrix[J,J] A_solve;
        if (length_index[t, i] > 0) {
          A_delta = A_delta_cache[length_index[t, i]];
          VCV = VCV_cache[length_index[t, i]];
          L_VCV = L_VCV_cache[length_index[t, i]];
          A_solve = A_solve_cache[length_index[t, i]];
        } else {
          A_delta = matrix_exp(A * ts[t, i]);
          VCV = Q_inf - quad_form_sym(Q_inf, A_delta');
          L_VCV = cholesky_decompose(VCV);
          A_solve = A \ add_diag(A_delta, -1);
        }
        drift_seg = L_VCV * z_drift[t, i-1];
        // if not a tip, add the drift parameter
        if (tip[t, i] == 0) {
          eta[t, node_seq[t, i]] = to_vector(
            A_delta * eta[t, parent[t, i]] + (A_solve * b) + drift_seg
          );
          VCV_tips[t, node_seq[t, i]] = diag_matrix(rep_vector(-99, J));
          L_VCV_tips[t, node_seq[t, i]] = diag_matrix(rep_vector(1.0, J));
        }
        // if is a tip, omit, we'll deal with it in the model block;
        else {
          eta[t, node_seq[t, i]] = to_vector(
            A_delta * eta[t, parent[t, i]] + (A_solve * b)
          );
          VCV_tips[t, node_seq[t, i]] = VCV;
          L_VCV_tips[t, node_seq[t, i]] = L_VCV;
        }
      }
    }
  }
  {{#tdrift}}
  for (t in 1:N_tree) {
    for (i in 1:N_tips) {
      tdrift[t,i] = L_VCV_tips[t, i] * to_vector(terminal_drift[t][i,]);
    }
  }
  {{/tdrift}}
  {{#dist_mat}}
  // distance covariance functions
  for (j in 1:J) {
    matrix[N_tips,N_tips] dist_cov;
    matrix[N_tips,N_tips] L_dist_cov;
    for ( i in 1:(N_tips-1) )
      for ( m in (i+1):N_tips ) {
        dist_cov[i,m] = {{dist_cov_code}};
        dist_cov[m,i] = dist_cov[i,m];
      }
    for ( q in 1:N_tips )
      dist_cov[q,q] = sigma_dist[j] + 0.01;
    L_dist_cov = cholesky_decompose(dist_cov);
    dist_v[,j] = L_dist_cov * dist_z[,j];
  }
  {{/dist_mat}}

}
