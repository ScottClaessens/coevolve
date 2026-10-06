transformed data {

  {{#variable_seq}}
  vector[to_int(N_obs - sum(col(miss, {{j}})))] obs{{j}}; // observed data for variable {{j}}
  {{/variable_seq}}
  {{#neg_binomial_seq}}
  real inv_overdisp{{j}}; // best guess for phi{{j}}
  {{/neg_binomial_seq}}
  {{#conditional_ncp}}
  array[N_obs, J] int tdrift_perm; // observed gaussian variables first
  array[N_obs] int n_tdrift_obs; // number of observed gaussian variables
  array[N_obs] int tdrift_perm_identity; // is tdrift_perm the identity?
  {{/conditional_ncp}}
  {{#variable_seq}}
  obs{{j}} = col(y, {{j}})[which_equal(col(miss, {{j}}), 0)];
  {{/variable_seq}}
  {{#neg_binomial_seq}}
  inv_overdisp{{j}} = (mean(obs{{j}})^2) / (sd(obs{{j}})^2 - mean(obs{{j}}));
  {{/neg_binomial_seq}}
  {{#conditional_ncp}}
  {
    array[J] int is_normal = { {{is_normal_flags}} };
    for (i in 1:N_obs) {
      int k = 0;
      for (j in 1:J) {
        if (is_normal[j] == 1 && miss[i,j] == 0) {
          k += 1;
          tdrift_perm[i,k] = j;
        }
      }
      n_tdrift_obs[i] = k;
      for (j in 1:J) {
        if (is_normal[j] == 0 || miss[i,j] == 1) {
          k += 1;
          tdrift_perm[i,k] = j;
        }
      }
      tdrift_perm_identity[i] = 1;
      for (j in 1:J) {
        if (tdrift_perm[i,j] != j) tdrift_perm_identity[i] = 0;
      }
    }
  }
  {{/conditional_ncp}}

}
