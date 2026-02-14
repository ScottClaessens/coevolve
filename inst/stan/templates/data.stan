data{
  int<lower=1> N_tips; // number of tips
  int<lower=1> N_tree; // number of trees
  int<lower=1> N_obs; // number of observations
  int<lower=2> J; // number of response traits
  int<lower=1> N_seg; // total number of segments in the trees
  array[N_tree, N_seg] int<lower=1> node_seq; // index of tree nodes
  array[N_tree, N_seg] int<lower=0> parent; // index of parent nodes
  array[N_tree, N_seg] real ts; // time since parent
  array[N_tree, N_seg] int<lower=0,upper=1> tip; // segment ends in tip
  array[J,J] int<lower=0,upper=1> effects_mat; // effects matrix
  int<lower=2> num_effects; // number of effects being estimated
  matrix[N_obs,J] y; // observed data
  matrix[N_obs,J] miss; // are data points missing?
{{#measurement_error}}
  matrix[N_obs,J] se; // squared standard errors
{{/measurement_error}}
  array[N_obs] int<lower=1> tip_id; // group index between 1 and N_tips
  int<lower=1> N_unique_lengths; // number of unique branch lengths
  array[N_unique_lengths] real unique_lengths; // unique branch lengths
  array[N_tree, N_seg] int<lower=0> length_index; // map segments to lengths
  array[N_tree, N_tips] int<lower=0> tip_to_seg; // map tips to segments
{{#dist_mat}}
  matrix[N_tips,N_tips] dist_mat; // distance matrix
{{/dist_mat}}
  int<lower=0,upper=1> prior_only; // should likelihood be ignored?
}
