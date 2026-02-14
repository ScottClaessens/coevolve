parameters{

  vector<upper=0>[J] A_diag; // autoregressive terms of A
  vector[num_effects - J] A_offdiag; // cross-lagged terms of A
{{#estimate_correlated_drift}}
  cholesky_factor_corr[J] L_R; // lower-tri choleksy decomp corr mat
{{/estimate_correlated_drift}}
  vector<lower=0>[J] Q_sigma; // std deviation parameters of the Q mat
  vector[J] b; // SDE intercepts
  array[N_tree] vector[J] eta_anc; // ancestral states
  array[N_tree, N_seg - 1] vector[J] z_drift; // stochastic drift
  array[N_tree] matrix[N_tips, J] terminal_drift; // drift for the tips
{{#ordered_seq}}
  ordered[{{num_cuts}}] c{{j}}; // cut points for variable {{j}}
{{/ordered_seq}}
{{#neg_binomial_seq}}
  real<lower=0> phi{{j}}; // neg binom inv overdispersion par for variable {{j}}
{{/neg_binomial_seq}}
{{#gamma_seq}}
  real<lower=0> shape{{j}}; // gamma shape par for variable {{j}}
{{/gamma_seq}}
{{#dist_mat}}
  matrix[N_tips,J] dist_z; // spatial covariance random effects
  vector<lower=0>[J] rho_dist; // covariance declining with distance
  vector<lower=0>[J] sigma_dist; // maximum covariance
{{/dist_mat}}
{{#repeated_measures}}
  matrix[J,N_obs] residual_z;
  vector<lower=0>[J] sigma_residual;
  cholesky_factor_corr[J] L_residual;
{{/repeated_measures}}

}
