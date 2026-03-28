#' Make model configuration for dynamic coevolutionary model
#'
#' Computes model structure flags and prior specifications consumed by
#' the JAX/NumPyro model builder. The returned list is embedded into the
#' data dict (via \code{embed_model_config}) so the model builder can
#' construct the model without any R-level string generation.
#'
#' @inheritParams coev_make_stancode
#'
#' @returns A named list of model configuration flags and prior
#'   specifications. Merged into the data list via
#'   \code{embed_model_config}.
#'
#' @author Erik Ringen \email{erikjacob.ringen@@uzh.ch}
#'
#' @seealso \code{\link{coev_make_stancode}}, \code{\link{coev_make_standata}},
#'   \code{\link{coev_fit}}
#'
#' @examples
#' cfg <- coev_make_model_config(
#'   data = authority$data,
#'   variables = list(
#'     political_authority = "ordered_logistic",
#'     religious_authority = "ordered_logistic"
#'   ),
#'   id = "language",
#'   tree = authority$phylogeny
#' )
#'
#' @export
coev_make_model_config <- function(data, variables, id, tree,
                                   effects_mat = NULL, complete_cases = FALSE,
                                   dist_mat = NULL, dist_cov = "exp_quad",
                                   measurement_error = NULL,
                                   prior = NULL, scale = TRUE,
                                   estimate_correlated_drift = TRUE,
                                   estimate_residual = TRUE,
                                   prior_only = FALSE) {

  run_checks(data, variables, id, tree, effects_mat, complete_cases, dist_mat,
             dist_cov, measurement_error, prior, scale,
             estimate_correlated_drift, estimate_residual,
             log_lik = FALSE, prior_only = prior_only)

  data <- as.data.frame(data)
  distributions <- as.character(variables)
  variables <- names(variables)
  J <- length(variables)

  priors <- list(
    b              = "std_normal()",
    eta_anc        = "std_normal()",
    A_offdiag      = "std_normal()",
    A_diag         = "std_normal()",
    L_R            = "lkj_corr_cholesky(4)",
    Q_sigma        = "std_normal()",
    c              = "normal(0, 2)",
    shape          = "gamma(0.01, 0.01)",
    sigma_dist     = "exponential(1)",
    rho_dist       = "exponential(5)",
    sigma_residual = "exponential(1)",
    L_residual     = "lkj_corr_cholesky(4)"
  )
  if (!is.null(prior)) {
    for (nm in names(prior)) priors[[nm]] <- prior[[nm]]
  }

  if (is.null(effects_mat)) {
    effects_mat <- matrix(TRUE, nrow = J, ncol = J,
                          dimnames = list(variables, variables))
  }
  effects_mat <- effects_mat[variables, variables]
  effects_mat_int <- +effects_mat
  n_offdiag <- sum(effects_mat_int) - J

  repeated <- any(duplicated(data[, id])) && estimate_residual
  tdrift <- repeated || !("normal" %in% distributions)
  residual_v <- repeated && !("normal" %in% distributions)
  normal_vars <- variables[distributions == "normal"]
  has_non_normal <- any(distributions != "normal")
  needs_terminal_drift <-
    tdrift ||
    (length(normal_vars) > 0 &&
     (has_non_normal || any(is.na(data[, normal_vars, drop = FALSE]))))

  ordered_j    <- integer(0)
  ordered_ncuts <- integer(0)
  if ("ordered_logistic" %in% distributions) {
    for (j in which(distributions == "ordered_logistic")) {
      num_cuts <- max(as.numeric(data[, variables[j]]), na.rm = TRUE) - 1L
      ordered_j     <- c(ordered_j, j)
      ordered_ncuts <- c(ordered_ncuts, num_cuts)
    }
  }

  lkj_eta_drift    <- as.numeric(parse_stan_prior(priors$L_R)$args[1])
  lkj_eta_residual <- as.numeric(parse_stan_prior(priors$L_residual)$args[1])

  prior_specs <- list(
    A_diag         = prior_spec_from_stan(priors$A_diag,         "upper_zero"),
    A_offdiag      = prior_spec_from_stan(priors$A_offdiag,      "none"),
    Q_sigma        = prior_spec_from_stan(priors$Q_sigma,         "lower_zero"),
    b              = prior_spec_from_stan(priors$b,               "none"),
    eta_anc        = prior_spec_from_stan(priors$eta_anc,         "none"),
    c              = prior_spec_from_stan(priors$c,               "none"),
    shape          = prior_spec_from_stan(priors$shape,           "lower_zero"),
    sigma_dist     = prior_spec_from_stan(priors$sigma_dist,      "lower_zero"),
    rho_dist       = prior_spec_from_stan(priors$rho_dist,        "lower_zero"),
    sigma_residual = prior_spec_from_stan(priors$sigma_residual,  "lower_zero")
  )
  if (!is.null(priors$phi)) {
    prior_specs$phi <- prior_spec_from_stan(priors$phi, "lower_zero")
  }

  list(
    distributions           = distributions,
    tdrift                  = as.integer(tdrift),
    repeated                = as.integer(repeated),
    needs_terminal_drift    = as.integer(needs_terminal_drift),
    residual_v              = as.integer(residual_v),
    estimate_correlated_drift = as.integer(estimate_correlated_drift),
    n_offdiag               = as.integer(n_offdiag),
    has_dist_mat            = as.integer(!is.null(dist_mat)),
    dist_cov_type           = if (!is.null(dist_mat)) dist_cov else "",
    has_measurement_error   = as.integer(!is.null(measurement_error)),
    ordered_j               = ordered_j,
    ordered_ncuts           = ordered_ncuts,
    nb_j0                   = which(
      distributions == "negative_binomial_softplus"
    ) - 1L,
    gamma_j0                = which(distributions == "gamma_log") - 1L,
    nonnormal_j0            = which(distributions != "normal") - 1L,
    normal_j0               = which(distributions == "normal") - 1L,
    lkj_eta_drift           = lkj_eta_drift,
    lkj_eta_residual        = lkj_eta_residual,
    prior_specs             = prior_specs
  )
}
