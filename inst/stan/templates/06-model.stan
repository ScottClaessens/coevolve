model{

  // priors
  b ~ {{prior_b}};
  for (t in 1:N_tree) {
    eta_anc[t] ~ {{prior_eta_anc}};
    for (i in 1:(N_seg - 1)) z_drift[t, i] ~ std_normal();
    {{#add_terminal_drift_prior}}
    to_vector(terminal_drift[t]) ~ std_normal();
    {{/add_terminal_drift_prior}}
  }
  A_offdiag ~ {{prior_A_offdiag}};
  A_diag ~ {{prior_A_diag}};
  Q_sigma ~ {{prior_Q_sigma}};
  {{#estimate_correlated_drift}}
  L_R ~ {{prior_L_R}};
  {{/estimate_correlated_drift}}
  {{#ordered_seq}}
  c{{j}} ~ {{prior_c}};
  {{/ordered_seq}}
  {{#prior_phi_default}}
  phi{{j}} ~ normal(inv_overdisp{{j}}, inv_overdisp{{j}});
  {{/prior_phi_default}}
  {{#prior_phi_manual}}
  phi{{j}} ~ {{prior_phi}};
  {{/prior_phi_manual}}
  {{#gamma_seq}}
  shape{{j}} ~ {{prior_shape}};
  {{/gamma_seq}}
  {{#dist_mat}}
  to_vector(dist_z) ~ std_normal();
  sigma_dist ~ {{prior_sigma_dist}};
  rho_dist ~ {{prior_rho_dist}};
  {{/dist_mat}}
  {{#add_priors_residual_sds_cors}}
  {{#normal_present}}
  for (i in 1:N_obs) {
    {{#is_normal}}
    if (miss[i,{{j}}] == 0) residual_z[{{j}},i] ~ std_normal();
    {{/is_normal}}
    {{#is_not_normal}}
    residual_z[{{j}},i] ~ std_normal();
    {{/is_not_normal}}
  }
  {{/normal_present}}
  {{#normal_absent}}
  to_vector(residual_z) ~ std_normal();
  {{/normal_absent}}
  sigma_residual ~ {{prior_sigma_residual}};
  L_residual ~ {{prior_L_residual}};
  {{/add_priors_residual_sds_cors}}

  // likelihood
  if (!prior_only) {
    for (i in 1:N_obs) {
      vector[N_tree] lp = rep_vector(0.0, N_tree);
      for (t in 1:N_tree) {
        {{#init_residuals}}
        // initialise residuals
        vector[J] residuals;
        {{/init_residuals}}
        {{#init_tdrift}}
        // initialise tdrift
        vector[J] tdrift;
        {{/init_tdrift}}
        {{#set_residuals}}
        // set residuals
        {{#is_normal}}
        if (miss[i,{{j}}] == 0) {
          residuals[{{j}}] = y[i,{{j}}] - ({{lmod}} + tdrift[t,tip_id[i]][{{j}}]);
        } else {
          residuals[{{j}}] = residual_z[{{j}},i];
        }
        {{/is_normal}}
        {{#is_not_normal}}
        residuals[{{j}}] = residual_z[{{j}},i];
        {{/is_not_normal}}
        {{#measurement_error}}
        matrix[J,J] residual_cov = diag_matrix(to_vector(se[i,])) + quad_form_diag(L_residual * L_residual', sigma_residual);
        {{/measurement_error}}
        lp[t] = multi_normal_cholesky_lpdf(residuals | rep_vector(0.0, J), {{cov_matrix}});
        {{/set_residuals}}
        {{#set_tdrift}}
        // set tdrift
        {{#is_normal}}
        if (miss[i,{{j}}] == 0) {
          tdrift[{{j}}] = y[i,{{j}}] - ({{lmod}});
          terminal_drift[t][tip_id[i],{{j}}] ~ std_normal();
        } else {
          tdrift[{{j}}] = terminal_drift[t][tip_id[i],{{j}}];
        }
        {{/is_normal}}
        {{#is_not_normal}}
        tdrift[{{j}}] = terminal_drift[t][tip_id[i],{{j}}];
        {{/is_not_normal}}
        lp[t] = multi_normal_cholesky_lpdf(tdrift | rep_vector(0.0, J), {{cov_matrix}});
        {{/set_tdrift}}
        {{#likelihoods}}
        if (miss[i,{{j}}] == 0) lp[t] += {{likelihood}};
        {{/likelihoods}}
      }
      target += log_sum_exp(lp);
    }
  }

}
