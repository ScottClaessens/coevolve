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
  {{#gps}}
  matrix[N_tips,J] dist_v; // distance covariance random effects
  {{/gps}}
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
    {{#j_equals_2}}
    // Fast cache (J=2): Cayley-Hamilton exp(At) = c0*I + c1t*A
    {
      real tr_A = A[1,1] + A[2,2];
      real det_A = A[1,1]*A[2,2] - A[1,2]*A[2,1];
      matrix[J,J] A_inv;
      A_inv[1,1] = A[2,2] / det_A;
      A_inv[1,2] = -A[1,2] / det_A;
      A_inv[2,1] = -A[2,1] / det_A;
      A_inv[2,2] = A[1,1] / det_A;
      // precomputed products for VCV = Q_inf - A_delta * Q_inf * A_delta'
      matrix[J,J] AQ_pre = A * Q_inf;
      matrix[J,J] S_pre = AQ_pre + Q_inf * A';
      matrix[J,J] T_pre = AQ_pre * A';
      for (u in 1:N_unique_lengths) {
        real t_u = unique_lengths[u];
        vector[2] ch = cayley_hamilton_coeffs(tr_A * t_u, det_A * t_u * t_u);
        real c0 = ch[1];
        real c1t = ch[2] * t_u;
        A_delta_cache[u] = add_diag(c1t * A, c0);
        real c0sq = c0 * c0;
        VCV_cache[u] = (1 - c0sq) * Q_inf - c0 * c1t * S_pre
                        - c1t * c1t * T_pre;
        // analytical Cholesky for 2x2 SPD
        {
          real v11 = VCV_cache[u][1,1];
          real v12 = VCV_cache[u][1,2];
          real l11 = sqrt(v11);
          real l21 = v12 / l11;
          L_VCV_cache[u][1,1] = l11;
          L_VCV_cache[u][1,2] = 0;
          L_VCV_cache[u][2,1] = l21;
          L_VCV_cache[u][2,2] = sqrt(VCV_cache[u][2,2] - l21 * l21);
        }
        // A_solve = A^{-1}*(A_delta - I) = (c0-1)*A^{-1} + c1t*I
        A_solve_cache[u] = add_diag((c0 - 1) * A_inv, c1t);
        for (i in 1:J) {
          for (j in 1:i) {
            real val = 0.5 * (A_solve_cache[u][i, j]
                              + A_solve_cache[u][j, i]);
            A_solve_cache[u][i, j] = val;
            A_solve_cache[u][j, i] = val;
          }
        }
      }
    }
    {{/j_equals_2}}
    {{#j_equals_3}}
    // Fast cache (J=3): Cayley-Hamilton exp(At) = c0*I + c1*At + c2*(At)^2
    {
      matrix[J,J] A_inv = inverse(A);
      // precomputed products for VCV = Q_inf - Ad * Q_inf * Ad'
      // Ad = c0*I + c1t*A + c2t2*A2, with A2 = A*A
      matrix[J,J] A2 = A * A;
      matrix[J,J] AQ = A * Q_inf;
      matrix[J,J] A2Q = A2 * Q_inf;
      // 6 precomputed terms for the quadratic form expansion
      matrix[J,J] P_00 = Q_inf;              // I*Q*I
      matrix[J,J] P_01 = AQ + Q_inf * A';    // A*Q*I + I*Q*A'
      matrix[J,J] P_02 = A2Q + Q_inf * A2';  // A2*Q*I + I*Q*A2'
      matrix[J,J] P_11 = AQ * A';            // A*Q*A'
      matrix[J,J] P_12 = A2Q * A' + AQ * A2';// A2*Q*A' + A*Q*A2'
      matrix[J,J] P_22 = A2Q * A2';          // A2*Q*A2'
      for (u in 1:N_unique_lengths) {
        real t_u = unique_lengths[u];
        matrix[J,J] M_u = A * t_u;
        vector[3] ch = cayley_hamilton_coeffs_3(M_u);
        real c0 = ch[1];
        real c1t = ch[2] * t_u;
        real c2t2 = ch[3] * t_u * t_u;
        // A_delta = c0*I + c1t*A + c2t2*A2
        A_delta_cache[u] = add_diag(c1t * A + c2t2 * A2, c0);
        // VCV via precomputed products (expand quadratic form)
        VCV_cache[u] = (1 - c0*c0) * P_00
                        - c0 * c1t * P_01
                        - c0 * c2t2 * P_02
                        - c1t * c1t * P_11
                        - c1t * c2t2 * P_12
                        - c2t2 * c2t2 * P_22;
        L_VCV_cache[u] = cholesky_decompose(VCV_cache[u]);
        // A_solve = A^{-1}*(A_delta - I) = (c0-1)*A^{-1} + c1t*I + c2t2*A
        A_solve_cache[u] = (c0 - 1) * A_inv + c2t2 * A;
        A_solve_cache[u][1,1] += c1t;
        A_solve_cache[u][2,2] += c1t;
        A_solve_cache[u][3,3] += c1t;
        for (i in 1:J) {
          for (j in 1:i) {
            real val = 0.5 * (A_solve_cache[u][i, j]
                              + A_solve_cache[u][j, i]);
            A_solve_cache[u][i, j] = val;
            A_solve_cache[u][j, i] = val;
          }
        }
      }
    }
    {{/j_equals_3}}
    {{#j_general}}
    // General cache computation (J >= 4)
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
    {{/j_general}}
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
  {{#gps}}
  // distance covariance functions
  {{#exact_gps}}
  for (j in 1:J) {
    matrix[N_tips, N_tips] dist_cov;
    dist_cov = {{dist_cov_function}}(coords, sigma_dist[j], rho_dist[j]);
    for (n in 1:N_tips) {
      dist_cov[n, n] += 1e-12;
    }
    dist_v[, j] = cholesky_decompose(dist_cov) * dist_z[, j];
  }
  {{/exact_gps}}
  {{#approximate_gps}}
  for (j in 1:J) {
    vector[NBgp] rgp = sqrt({{dist_cov_function}}(slambda, sigma_dist[j], rho_dist[j])) .* dist_z[, j];
    dist_v[, j] = Xgp * rgp;
  }
  {{/approximate_gps}}
  {{/gps}}

}
