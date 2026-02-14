transformed data {
{{#variable_seq}}
  vector[to_int(N_obs - sum(col(miss, {{j}})))] obs{{j}}; // observed data {{j}}
{{/variable_seq}}
{{#neg_binomial_seq}}
  real inv_overdisp{{j}}; // best guess for phi{{j}}
{{/neg_binomial_seq}}
{{#variable_seq}}
  obs{{j}} = col(y, {{j}})[which_equal(col(miss, {{j}}), 0)];
{{/variable_seq}}
{{#neg_binomial_seq}}
  inv_overdisp{{j}} = (mean(obs{{j}})^2) / (sd(obs{{j}})^2 - mean(obs{{j}}));
{{/neg_binomial_seq}}

}
