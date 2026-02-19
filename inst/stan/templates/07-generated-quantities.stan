generated quantities{

  {{#log_lik}}
  vector[N_obs*J] log_lik; // log-likelihood
  {{/log_lik}}
  array[N_tree,N_obs,J] real yrep; // predictive checks
  {{#estimate_correlated_drift}}
  matrix[J,J] cor_R; // correlated drift
  cor_R = multiply_lower_tri_self_transpose(L_R);
  {{/estimate_correlated_drift}}
  {{#repeated}}
  matrix[J,J] cor_residual; // residual correlations
  cor_residual = multiply_lower_tri_self_transpose(L_residual);
  {{/repeated}}
  {
    {{#log_lik}}
    matrix[N_obs,J] log_lik_temp = rep_matrix(0.0, N_obs, J);
    {{/log_lik}}
    array[N_tree,N_tips] vector[J] terminal_drift_rep;
    for (i in 1:N_tips) {
      for (t in 1:N_tree) {
        for (j in 1:J) terminal_drift_rep[t,i][j] = normal_rng(0, 1);
        terminal_drift_rep[t,i] = cholesky_decompose({{terminal_drift_cov_matrix}}) * terminal_drift_rep[t,i];
      }
    }
    for (i in 1:N_obs) {
      {{#log_lik}}
      array[N_tree,N_obs,J] real lp = rep_array(0.0, N_tree, N_obs, J);
      {{/log_lik}}
      for (t in 1:N_tree) {
        {{#normal_present_and_log_lik}}
        vector[J] mu_cond;
        vector[J] sigma_cond;
        {{/normal_present_and_log_lik}}
        {{#init_residuals}}
        vector[J] residuals;
        {{/init_residuals}}
        {{#init_tdrifts}}
        vector[J] tdrifts;
        {{/init_tdrifts}}
        {{#repeated}}
        vector[J] residuals_rep;
        for (j in 1:J) residuals_rep[j] = normal_rng(0, 1);
        residuals_rep = diag_pre_multiply(sigma_residual, L_residual) * residuals_rep;
        {{/repeated}}
        {{#set_residuals}}
        {{#is_normal}}
        residuals[{{j}}] = y[i][{{j}}] - ({{lmod}} + tdrift[t,tip_id[i]][{{j}}]);
        {{/is_normal}}
        {{#is_not_normal_and_normal_present}}
        residuals[{{j}}] = residual_z[{{j}},i];
        {{/is_not_normal_and_normal_present}}
        {{#is_not_normal_and_normal_absent}}
        residuals[{{j}}] = residual_v[i,{{j}}];
        {{/is_not_normal_and_normal_absent}}
        {{/set_residuals}}
        {{#set_tdrifts}}
        {{#is_normal}}
        tdrifts[{{j}}] = y[i][{{j}}] - ({{lmod}});
        {{/is_normal}}
        {{#is_not_normal_and_normal_present}}
        tdrifts[{{j}}] = terminal_drift[t,tip_id[i]][{{j}}];
        {{/is_not_normal_and_normal_present}}
        {{#is_not_normal_and_normal_absent}}
        tdrifts[{{j}}] = tdrift[t,tip_id[i]][{{j}}];
        {{/is_not_normal_and_normal_absent}}
        {{/set_tdrifts}}
        {{#set_mu_cond_and_sigma_cond}}
        {{#repeated}}
        {{#measurement_error}}
        matrix[J,J] residual_cov = diag_matrix(to_vector(se[i,])) + quad_form_diag(L_residual * L_residual', sigma_residual);
        matrix[J,J] cov_inv = inverse_spd(residual_cov);
        {{/measurement_error}}
        {{#no_measurement_error}}
        matrix[J,J] cov_inv = chol2inv(diag_pre_multiply(sigma_residual, L_residual));
        {{/no_measurement_error}}
        mu_cond = residuals - (cov_inv * residuals) ./ diagonal(cov_inv);
        sigma_cond = sqrt(1 / diagonal(cov_inv));
        {{/repeated}}
        {{#no_repeated}}
        {{#measurement_error}}
        matrix[J,J] cov_inv = inverse_spd(add_diag(VCV_tips[t, tip_id[i]], se[i,]));
        {{/measurement_error}}
        {{#no_measurement_error}}
        matrix[J,J] cov_inv = inverse_spd(VCV_tips[t, tip_id[i]]);
        {{/no_measurement_error}}
        mu_cond = tdrifts - (cov_inv * tdrifts) ./ diagonal(cov_inv);
        sigma_cond = sqrt(1 / diagonal(cov_inv));
        {{/no_repeated}}
        {{/set_mu_cond_and_sigma_cond}}
        {{#posterior_predictions_and_log_likelihoods}}
        {{#log_lik}}
        if (miss[i,{{j}}] == 0) lp[t,i,{{j}}] = {{log_lik_statement}};
        {{/log_lik}}
        {{#yrep}}
        yrep[t,i,{{j}}] = {{yrep_statement}};
        {{/yrep}}
        {{/posterior_predictions_and_log_likelihoods}}
      }
      {{#log_lik}}
      for (j in 1:J) log_lik_temp[i,j] += log_sum_exp(lp[,i,j]);
      {{/log_lik}}
    }
    {{#log_lik}}
    log_lik = to_vector(log_lik_temp);
    {{/log_lik}}
  }

}
