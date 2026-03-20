#' Make PyMC/PyTensor code for dynamic coevolutionary model
#'
#' Generate a Python string containing a PyMC model function equivalent to
#' the Stan model produced by \code{\link{coev_make_stancode}}. The generated
#' code can be executed via \pkg{reticulate} and sampled with PyMC's NUTS.
#' When sampled with nutpie's JAX backend, the graph compiles to XLA for
#' efficient CPU execution.
#'
#' @inheritParams coev_make_stancode
#'
#' @returns A character string containing Python code that defines a
#'   \code{build_model(data)} function returning a PyMC model.
#'
#' @author Erik Ringen \email{erikjacob.ringen@@uzh.ch}
#'
#' @seealso \code{\link{coev_make_stancode}}, \code{\link{coev_make_standata}},
#'   \code{\link{coev_fit}}
#'
#' @examples
#' pymc_code <- coev_make_pymc(
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
coev_make_pymc <- function(data, variables, id, tree,
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
  num_effects <- sum(effects_mat_int)
  n_offdiag <- num_effects - J

  repeated <- any(duplicated(data[, id])) && estimate_residual
  tdrift <- repeated || !("normal" %in% distributions)
  residual_v <- repeated && !("normal" %in% distributions)
  normal_vars <- variables[distributions == "normal"]
  needs_terminal_drift <-
    tdrift || (length(normal_vars) > 0 &&
      any(is.na(data[, normal_vars, drop = FALSE])))

  ordered_info <- list()
  if ("ordered_logistic" %in% distributions) {
    for (j in which(distributions == "ordered_logistic")) {
      num_cuts <- max(as.numeric(data[, variables[j]]), na.rm = TRUE) - 1
      ordered_info[[as.character(j)]] <- list(j = j, num_cuts = num_cuts)
    }
  }

  code <- c(
    pymc_header(J),
    "",
    pymc_model_fn(
      J = J,
      distributions = distributions,
      variables = variables,
      data = data,
      id = id,
      effects_mat_int = effects_mat_int,
      n_offdiag = n_offdiag,
      priors = priors,
      estimate_correlated_drift = estimate_correlated_drift,
      tdrift = tdrift,
      residual_v = residual_v,
      needs_terminal_drift = needs_terminal_drift,
      repeated = repeated,
      dist_mat = dist_mat,
      dist_cov = dist_cov,
      measurement_error = measurement_error,
      ordered_info = ordered_info
    )
  )

  paste(code, collapse = "\n")
}

# ---------------------------------------------------------------------------
# Internal helpers for Python code generation (PyMC)
# ---------------------------------------------------------------------------

#' Append a keyword argument before the closing `)` of a pm.Dist(...) call
#'
#' PyMC's default `initval` for vectors of Normals is often all zeros. For
#' ordered cutpoints we apply `pt.sort(c_raw)`; identical raw values yield tied
#' cutpoints and break `OrderedLogistic` (-inf logp). Constrained OU parameters
#' should not start exactly on boundaries (e.g. A_diag = 0, Q_sigma = 0).
#'
#' @noRd
append_pymc_kwarg <- function(call_str, kwarg) {
  sub("\\)$", paste0(", ", kwarg, ")"), call_str)
}

#' @noRd
pymc_header <- function(J) {
  c(
    sprintf("# Generated with coevolve %s (PyMC backend)",
            utils::packageVersion("coevolve")),
    "import pymc as pm",
    "import pytensor.tensor as pt",
    "import numpy as np",
    "",
    "",
    "def diag_mat(v):",
    "    return v[:, None] * pt.eye(v.shape[0])",
    "",
    "",
    "def mvn_chol_logp(value, chol):",
    "    z = pt.linalg.solve_triangular(chol, value[..., None], lower=True)[..., 0]",
    "    logdet = pt.sum(pt.log(pt.diagonal(chol, axis1=-2, axis2=-1)), axis=-1)",
    "    quad = pt.sum(z ** 2, axis=-1)",
    "    jdim = pt.cast(value.shape[-1], value.dtype)",
    "    return -0.5 * (jdim * np.log(2.0 * np.pi) + quad) - logdet",
    "",
    "",
    "def ksolve(A, Q):",
    "    \"\"\"Solve AX + XA^T = -Q for X (continuous Lyapunov).\"\"\"",
    sprintf("    _I = pt.eye(%d)", J),
    "    M = (A[:, :, None, None] * _I[None, None, :, :] +",
    "         _I[:, :, None, None] * A[None, None, :, :])",
    sprintf("    M = M.transpose(0, 2, 1, 3).reshape((%d, %d))", J * J, J * J),
    sprintf("    X = pt.linalg.solve(M, -Q.reshape((-1,))).reshape((%d, %d))", J, J),
    "    return 0.5 * (X + X.T)",
    "",
    "",
    pymc_lkj_corr_chol_fn(J)
  )
}


#' Generate LKJ correlation Cholesky via stick-breaking (JAX-compatible)
#'
#' Uses \code{pt.stacklists} for matrix assembly instead of
#' \code{pt.set_subtensor} loops, producing a much smaller graph.
#' @noRd
pymc_lkj_corr_chol_fn <- function(J) {
  if (J < 2) return(character(0))

  n_corr <- J * (J - 1L) / 2L
  entry <- function(i, j) sprintf("_e_%d_%d", i, j)

  lines <- c(
    "def sample_lkj_cholesky(name, n, eta):",
    "    \"\"\"Sample Cholesky of correlation matrix with LKJ(eta) prior.\"\"\"",
    "    if n == 1:",
    "        return pt.eye(1)",
    sprintf("    raw = pm.Flat(f'_{name}_raw', shape=%d)", n_corr),
    "    z = pt.tanh(raw)"
  )

  # Emit scalar entries for L (stick-breaking parameterization)
  idx <- 0L
  for (i in seq_len(J) - 1L) {
    for (j in 0L:i) {
      if (i == 0L && j == 0L) next
      if (i == j) {
        if (i == 1L) {
          sq_terms <- sprintf("%s**2", entry(i, 0L))
        } else {
          sq_terms <- paste(
            vapply(0L:(j - 1L), function(k) sprintf("%s**2", entry(i, k)),
                   character(1)),
            collapse = " + ")
        }
        lines <- c(lines, sprintf(
          "    %s = pt.sqrt(pt.maximum(1.0 - (%s), 1e-10))", entry(i, j), sq_terms))
      } else if (j == 0L) {
        lines <- c(lines, sprintf("    %s = z[%d]", entry(i, j), idx))
        idx <- idx + 1L
      } else {
        sq_terms <- paste(
          vapply(0L:(j - 1L), function(k) sprintf("%s**2", entry(i, k)),
                 character(1)),
          collapse = " + ")
        lines <- c(lines, sprintf(
          "    %s = z[%d] * pt.sqrt(pt.maximum(1.0 - (%s), 1e-10))",
          entry(i, j), idx, sq_terms))
        idx <- idx + 1L
      }
    }
  }

  # Build matrix via stacklists (single op, no set_subtensor)
  row_strs <- character(J)
  for (i in seq_len(J) - 1L) {
    elems <- vapply(seq_len(J) - 1L, function(j) {
      if (j > i) "pt.constant(0.0)"
      else if (i == 0L && j == 0L) "pt.constant(1.0)"
      else entry(i, j)
    }, character(1))
    row_strs[i + 1L] <- sprintf("[%s]", paste(elems, collapse = ", "))
  }
  lines <- c(lines, sprintf("    L = pt.stacklists([%s])",
                             paste(row_strs, collapse = ", ")))

  # LKJ log prior + Jacobian
  lkj_terms <- character(0)
  for (i in 1L:(J - 1L)) {
    coeff <- J - i - 1L + 2L  # (n - i - 1 + 2*(eta - 1)) simplified at eta=eta
    lkj_terms <- c(lkj_terms, sprintf(
      "(%d + 2*(eta - 1)) * pt.log(%s)", J - i - 1L, entry(i, i)))
  }
  lines <- c(lines, sprintf("    lkj_logp = %s",
                             paste(lkj_terms, collapse = " + ")))
  lines <- c(lines,
    "    jac = pt.sum(pt.log(1.0 - z**2 + 1e-10))",
    "    pm.Potential(f'_{name}_lkj_prior', lkj_logp + jac)",
    "    return L",
    ""
  )
  lines
}

#' Build the complete build_model() function body
#' @noRd
pymc_model_fn <- function(J, distributions, variables, data, id,
                          effects_mat_int, n_offdiag, priors,
                          estimate_correlated_drift, tdrift, residual_v,
                          needs_terminal_drift, repeated, dist_mat, dist_cov,
                          measurement_error, ordered_info) {

  I1 <- "    "
  I2 <- "        "
  I3 <- "            "
  I4 <- "                "
  I5 <- "                    "

  lines <- c()
  a <- function(...) lines <<- c(lines, ...)

  a("def build_model(data, compile_mode='cpu'):")
  a(paste0(I1, "with pm.Model() as model:"))

  # ---- Unpack data ----
  a(paste0(I2, "J = data['J']"))
  a(paste0(I2, "N_tips = data['N_tips']"))
  a(paste0(I2, "N_tree = data['N_tree']"))
  a(paste0(I2, "N_obs = data['N_obs']"))
  a(paste0(I2, "N_seg = data['N_seg']"))
  a(paste0(I2, "node_seq = data['node_seq']"))
  a(paste0(I2, "parent = data['parent']"))
  a(paste0(I2, "ts = data['ts']"))
  a(paste0(I2, "tip = data['tip']"))
  a(paste0(I2, "effects_mat = data['effects_mat']"))
  a(paste0(I2, "y = data['y']"))
  a(paste0(I2, "miss = data['miss']"))
  a(paste0(I2, "tip_id = data['tip_id']"))
  a(paste0(I2, "N_unique_lengths = data['N_unique_lengths']"))
  a(paste0(I2, "unique_lengths = data['unique_lengths']"))
  a(paste0(I2, "length_index = data['length_index']"))
  a(paste0(I2, "prior_only = data['prior_only']"))
  a(paste0(I2, "tip_to_seg = data['tip_to_seg']"))
  a(paste0(I2, "internal_mask = (np.asarray(tip[:, 1:], dtype=np.int32) == 0).astype(np.int32)"))
  a(paste0(I2, "internal_slot = np.cumsum(internal_mask, axis=1) - 1"))
  a(paste0(I2, "internal_slot[internal_mask == 0] = -1"))
  a(paste0(I2, "N_internal = int(np.max(np.sum(internal_mask, axis=1))) if N_seg > 1 else 0"))
  if (!is.null(dist_mat)) {
    a(paste0(I2, "dist_mat_data = data['dist_mat']"))
  }
  if (!is.null(measurement_error)) {
    a(paste0(I2, "se = data['se']"))
  }
  has_count <- any(distributions %in%
                     c("poisson_softplus", "negative_binomial_softplus"))
  if (has_count) {
    a(paste0(I2, "obs_means = data['obs_means']"))
    a(paste0(I2, "inv_overdisp = data['inv_overdisp']"))
  }
  a("")

  # ---- Parameters ----
  a(paste0(I2, "# === Parameters ==="))

  # A_diag (upper=0): must start strictly negative for stable OU / Lyapunov
  a_diag_tmpl <- stan_prior_to_pymc(priors$A_diag, "upper_zero", "J")
  a_diag_call <- append_pymc_kwarg(
    sprintf(a_diag_tmpl, "'A_diag'"),
    "initval=-0.5 * np.ones(J, dtype=np.float64)"
  )
  a(paste0(I2, sprintf("A_diag = %s", a_diag_call)))

  # A_offdiag
  if (n_offdiag > 0) {
    a_offdiag_tmpl <- stan_prior_to_pymc(priors$A_offdiag, "none",
                                          as.character(n_offdiag))
    a(paste0(I2, sprintf("A_offdiag = %s",
                          sprintf(a_offdiag_tmpl, "'A_offdiag'"))))
  }

  # L_R (Cholesky correlation)
  if (estimate_correlated_drift) {
    lkj_eta <- parse_stan_prior(priors$L_R)$args[1]
    a(paste0(I2, sprintf(
      "L_R = sample_lkj_cholesky('L_R', J, %s)", lkj_eta)))
  }

  # Q_sigma (lower=0): keep away from zero so Q is positive definite
  q_sigma_tmpl <- stan_prior_to_pymc(priors$Q_sigma, "lower_zero", "J")
  q_sigma_call <- append_pymc_kwarg(
    sprintf(q_sigma_tmpl, "'Q_sigma'"),
    "initval=0.5 * np.ones(J, dtype=np.float64)"
  )
  a(paste0(I2, sprintf("Q_sigma = %s", q_sigma_call)))

  # b
  b_tmpl <- stan_prior_to_pymc(priors$b, "none", "J")
  a(paste0(I2, sprintf("b = %s", sprintf(b_tmpl, "'b'"))))

  # eta_anc
  eta_anc_tmpl <- stan_prior_to_pymc(priors$eta_anc, "none", "(N_tree, J)")
  a(paste0(I2, sprintf("eta_anc = %s",
                        sprintf(eta_anc_tmpl, "'eta_anc'"))))

  # z_drift
  a(paste0(I2,
    "z_drift = pm.Normal('z_drift', mu=0.0, sigma=1.0, shape=(N_tree, N_internal, J))"))

  # terminal_drift
  if (needs_terminal_drift) {
    a(paste0(I2,
      "terminal_drift = pm.Normal('terminal_drift', mu=0.0, sigma=1.0, shape=(N_tree, N_tips, J))"))
  }

  # Ordered cutpoints: raw must not initialize to identical values (sort → ties
  # breaks OrderedLogistic); use strictly increasing initval before sort
  for (info in ordered_info) {
    j <- info$j
    nc <- info$num_cuts
    c_tmpl <- stan_prior_to_pymc(priors$c, "none", as.character(nc))
    c_code <- sprintf(c_tmpl, sprintf("'c%d_raw'", j))
    c_code <- append_pymc_kwarg(
      c_code,
      sprintf("initval=np.linspace(-1.0, 1.0, %d, dtype=np.float64)", nc)
    )
    a(paste0(I2, sprintf(
      "c%d_raw = %s", j, c_code)))
    a(paste0(I2, sprintf(
      "c%d = pm.Deterministic('c%d', pt.sort(c%d_raw))", j, j, j)))
  }

  # Negative binomial phi
  nb_idx <- which(distributions == "negative_binomial_softplus")
  for (j in nb_idx) {
    if (is.null(priors$phi)) {
      a(paste0(I2, sprintf(
        "phi%d = pm.TruncatedNormal('phi%d', mu=inv_overdisp[%d], sigma=inv_overdisp[%d], lower=0.0)",
        j, j, j - 1, j - 1)))
    } else {
      phi_tmpl <- stan_prior_to_pymc(priors$phi, "lower_zero")
      a(paste0(I2, sprintf("phi%d = %s", j,
                            sprintf(phi_tmpl, sprintf("'phi%d'", j)))))
    }
  }

  # Gamma shape
  gamma_idx <- which(distributions == "gamma_log")
  for (j in gamma_idx) {
    shape_tmpl <- stan_prior_to_pymc(priors$shape, "lower_zero")
    a(paste0(I2, sprintf("shape%d = %s", j,
                          sprintf(shape_tmpl, sprintf("'shape%d'", j)))))
  }

  # GP parameters
  if (!is.null(dist_mat)) {
    a(paste0(I2,
      "dist_z = pm.Normal('dist_z', mu=0.0, sigma=1.0, shape=(N_tips, J))"))
    sigma_dist_tmpl <- stan_prior_to_pymc(priors$sigma_dist, "lower_zero", "J")
    a(paste0(I2, sprintf("sigma_dist = %s",
                          sprintf(sigma_dist_tmpl, "'sigma_dist'"))))
    rho_dist_tmpl <- stan_prior_to_pymc(priors$rho_dist, "lower_zero", "J")
    a(paste0(I2, sprintf("rho_dist = %s",
                          sprintf(rho_dist_tmpl, "'rho_dist'"))))
  }

  # Residual parameters
  if (repeated) {
    a(paste0(I2,
      "residual_z = pm.Normal('residual_z', mu=0.0, sigma=1.0, shape=(J, N_obs))"))
    sigma_res_tmpl <- stan_prior_to_pymc(
      priors$sigma_residual, "lower_zero", "J")
    a(paste0(I2, sprintf("sigma_residual = %s",
                          sprintf(sigma_res_tmpl, "'sigma_residual'"))))
    lkj_eta_res <- parse_stan_prior(priors$L_residual)$args[1]
    a(paste0(I2, sprintf(
      "L_residual = sample_lkj_cholesky('L_residual', J, %s)",
      lkj_eta_res)))
  }
  a("")

  # ---- Transformed Parameters ----
  a(paste0(I2, "# === Transformed Parameters ==="))

  # Build A matrix
  a(paste0(I2, "A_mat = diag_mat(A_diag)"))
  if (n_offdiag > 0) {
    a(paste0(I2, "ticker = 0"))
    a(paste0(I2, "for i in range(J):"))
    a(paste0(I3, "for j in range(J):"))
    a(paste0(I4, "if i != j:"))
    a(paste0(I5, "if effects_mat[i, j] == 1:"))
    a(paste0(I5, "    A_mat = pt.set_subtensor(A_mat[i, j], A_offdiag[ticker])"))
    a(paste0(I5, "    ticker += 1"))
  }
  a(paste0(I2, "pm.Deterministic('A', A_mat)"))
  a("")

  # Build Q matrix
  if (estimate_correlated_drift) {
    a(paste0(I2,
      "Q = pt.dot(diag_mat(Q_sigma), pt.dot(pt.dot(L_R, L_R.T), diag_mat(Q_sigma)))"))
    a(paste0(I2, "pm.Deterministic('cor_R', pt.dot(L_R, L_R.T))"))
  } else {
    a(paste0(I2, "Q = diag_mat(Q_sigma ** 2)"))
  }
  a(paste0(I2, "pm.Deterministic('Q', Q)"))
  a("")

  # Q_inf
  a(paste0(I2, "Q_inf = ksolve(A_mat, Q)"))
  a("")

  # Cache for unique branch lengths (batched — single vectorized op)
  a(paste0(I2, "eye_J = pt.eye(J)"))
  a(paste0(I2, "_A_dt = A_mat[None, :, :] * (unique_lengths[:, None, None] / 256.0)"))
  a(paste0(I2, "_eye_b = pt.broadcast_to(eye_J[None, :, :], _A_dt.shape)"))
  a(paste0(I2, "_S = _eye_b"))
  a(paste0(I2, "_T = _eye_b"))
  a(paste0(I2, "for _k in range(1, 13):"))
  a(paste0(I3, "_T = pt.matmul(_T, _A_dt) * (1.0 / _k)"))
  a(paste0(I3, "_S = _S + _T"))
  a(paste0(I2, "for _ in range(8):"))
  a(paste0(I3, "_S = pt.matmul(_S, _S)"))
  a(paste0(I2, "A_delta_cache = _S"))
  a(paste0(I2, "_Qi = Q_inf[None, :, :]"))
  a(paste0(I2, "_V = _Qi - pt.matmul(A_delta_cache, pt.matmul(_Qi, A_delta_cache.transpose(0, 2, 1)))"))
  a(paste0(I2, "VCV_cache = 0.5 * (_V + _V.transpose(0, 2, 1))"))
  a(paste0(I2, "L_VCV_cache = pt.linalg.cholesky(VCV_cache)"))
  a(paste0(I2, "_A_inv = pt.linalg.solve(A_mat, eye_J)"))
  a(paste0(I2, "_As = pt.matmul(_A_inv[None, :, :], A_delta_cache - _eye_b)"))
  a(paste0(I2, "A_solve_cache = 0.5 * (_As + _As.transpose(0, 2, 1))"))
  a(paste0(I2, "b_delta_cache = pt.matmul(A_solve_cache, b[None, :, None])[:, :, 0]"))
  a("")

  # Tree traversal — level-batched (O(depth) graph ops instead of O(N_seg))
  a(paste0(I2, "# === Tree traversal (level-batched) ==="))
  a(paste0(I2, "n_levels = data['n_levels']"))
  a(paste0(I2, "max_level_size = data['max_level_size']"))
  a(paste0(I2, "root_ids = np.atleast_1d(np.asarray(data['root_ids'], dtype=np.int32))"))
  a(paste0(I2, "level_node_ids = np.asarray(data['level_node_ids'], dtype=np.int32).reshape((N_tree, n_levels, max_level_size))"))
  a(paste0(I2, "level_parent_ids = np.asarray(data['level_parent_ids'], dtype=np.int32).reshape((N_tree, n_levels, max_level_size))"))
  a(paste0(I2, "level_length_idx_data = np.asarray(data['level_length_idx'], dtype=np.int32).reshape((N_tree, n_levels, max_level_size))"))
  a(paste0(I2, "level_is_internal = np.asarray(data['level_is_internal'], dtype=np.int32).reshape((N_tree, n_levels, max_level_size))"))
  a(paste0(I2, "level_drift_idx = np.asarray(data['level_drift_idx'], dtype=np.int32).reshape((N_tree, n_levels, max_level_size))"))
  a(paste0(I2, "eta_trees = []"))
  a(paste0(I2, "tip_L_VCV_trees = []"))
  a(paste0(I2, "for t in range(int(N_tree)):"))
  a(paste0(I3, "eta = pt.zeros((N_seg, J))"))
  a(paste0(I3, "eta = pt.set_subtensor(eta[int(root_ids[t])], eta_anc[t])"))
  a(paste0(I3, "for _l in range(int(n_levels)):"))
  a(paste0(I4, "_nids = np.asarray(level_node_ids[t, _l], dtype=np.int32)"))
  a(paste0(I4, "_pids = np.asarray(level_parent_ids[t, _l], dtype=np.int32)"))
  a(paste0(I4, "_li = np.asarray(level_length_idx_data[t, _l], dtype=np.int32)"))
  a(paste0(I4, "_is_int = pt.as_tensor_variable(np.asarray(level_is_internal[t, _l], dtype=np.float64))"))
  a(paste0(I4, "_didx = np.asarray(level_drift_idx[t, _l], dtype=np.int32)"))
  a(paste0(I4, "_pv = eta[_pids]"))
  a(paste0(I4, "_Ad = A_delta_cache[_li]"))
  a(paste0(I4, "_bd = b_delta_cache[_li]"))
  a(paste0(I4, "_base = pt.matmul(_Ad, _pv[:, :, None])[:, :, 0] + _bd"))
  a(paste0(I4, "_Lv = L_VCV_cache[_li]"))
  a(paste0(I4, "_dv = z_drift[t, _didx]"))
  a(paste0(I4, "_noise = pt.matmul(_Lv, _dv[:, :, None])[:, :, 0]"))
  a(paste0(I4, "_vals = _base + _is_int[:, None] * _noise"))
  a(paste0(I4, "eta = pt.set_subtensor(eta[_nids], _vals)"))
  a(paste0(I3, "eta = pt.set_subtensor(eta[int(root_ids[t])], eta_anc[t])"))
  a(paste0(I3, "eta_trees.append(eta)"))
  a(paste0(I3, "# Vectorized tip L_VCV"))
  a(paste0(I3, "_tip_li = np.asarray(length_index[t][tip_to_seg[t]], dtype=np.int32)"))
  a(paste0(I3, "tip_L_VCV_trees.append(L_VCV_cache[_tip_li])"))
  a("")

  # tdrift (vectorised batched mat-vec)
  if (tdrift) {
    a(paste0(I2, "tdrift_trees = []"))
    a(paste0(I2, "for t in range(N_tree):"))
    a(paste0(I3, "L_tips = tip_L_VCV_trees[t]"))
    a(paste0(I3, "td = terminal_drift[t]"))
    a(paste0(I3,
      "tdrift_trees.append(pt.matmul(L_tips, td[:, :, None])[:, :, 0])"))
    a("")
  }

  # residual_v
  if (residual_v) {
    a(paste0(I2,
      "residual_v = pt.dot(pt.dot(diag_mat(sigma_residual), L_residual), residual_z).T"))
    a("")
  }

  # GP distance covariance (vectorized kernel matrix)
  if (!is.null(dist_mat)) {
    a(paste0(I2, "dist_v_list = []"))
    a(paste0(I2, "_dm = pt.as_tensor_variable(np.array(dist_mat_data, dtype=np.float64))"))
    a(paste0(I2, "for j_gp in range(J):"))
    gp_kernel_full <- switch(dist_cov,
      "exp_quad" = "sigma_dist[j_gp] * pt.exp(-(_dm ** 2) / rho_dist[j_gp])",
      "exponential" = "sigma_dist[j_gp] * pt.exp(-_dm / rho_dist[j_gp])",
      "matern32" = paste0(
        "sigma_dist[j_gp] * (1.0 + (pt.sqrt(3.0) * _dm) / rho_dist[j_gp]) * ",
        "pt.exp(-(pt.sqrt(3.0) * _dm) / rho_dist[j_gp])")
    )
    a(paste0(I3, sprintf("dc = %s", gp_kernel_full)))
    a(paste0(I3, "dc = pt.fill_diagonal(dc, sigma_dist[j_gp] + 0.01)"))
    a(paste0(I3, "L_dc = pt.linalg.cholesky(dc)"))
    a(paste0(I3, "dist_v_list.append(pt.dot(L_dc, dist_z[:, j_gp]))"))
    a(paste0(I2, "dist_v = pt.stack(dist_v_list, axis=1)"))
    a("")
  }

  # Residual correlation deterministic
  if (repeated) {
    a(paste0(I2,
      "pm.Deterministic('cor_residual', pt.dot(L_residual, L_residual.T))"))
    a("")
  }

  # ---- Vectorised Likelihood ----
  a(paste0(I2, "# === Likelihood (vectorised) ==="))
  a(paste0(I2, "if not prior_only:"))
  a(paste0(I3, "y_pt = pt.as_tensor_variable(np.array(y, dtype=np.float64))"))
  a(paste0(I3, "miss_pt = pt.as_tensor_variable(np.array(miss, dtype=np.int32))"))
  a(paste0(I3, "tid = np.array(tip_id, dtype=np.int32)"))
  a("")

  has_normal <- "normal" %in% distributions

  a(paste0(I3, "tree_lps = []"))
  a(paste0(I3, "for t in range(N_tree):"))
  a(paste0(I4, "eta_t = eta_trees[t]"))
  a(paste0(I4, "eta_obs = eta_t[tid]"))

  # Build vectorised linear model (base) per variable
  lmod_base_fn <- function(j0) {
    base <- sprintf("eta_obs[:, %d]", j0)
    if (!is.null(dist_mat)) {
      base <- sprintf("(%s + dist_v[tid, %d])", base, j0)
    }
    base
  }

  # --- Normal MvNormal block (vectorised) ---
  if (has_normal && !repeated) {
    normal_idx <- which(distributions == "normal")
    not_normal_idx <- which(distributions != "normal")

    a(paste0(I4, "L_cov_raw = tip_L_VCV_trees[t][tid]"))
    if (!is.null(measurement_error)) {
      a(paste0(I4,
        "VCV_obs = pt.matmul(L_cov_raw, L_cov_raw.transpose(0, 2, 1))"))
      a(paste0(I4,
        "se_diag = pt.as_tensor_variable(np.array(se, dtype=np.float64))[:, :, None] * pt.eye(J)[None, :, :]"))
      a(paste0(I4, "L_cov_obs = pt.linalg.cholesky(VCV_obs + se_diag)"))
    } else {
      a(paste0(I4, "L_cov_obs = L_cov_raw"))
    }

    if (needs_terminal_drift) {
      a(paste0(I4, "tdrift_vec = terminal_drift[t][tid]"))
    } else {
      a(paste0(I4, "tdrift_vec = pt.zeros((N_obs, J))"))
    }
    for (j in normal_idx) {
      j0 <- j - 1
      missing_expr <- if (needs_terminal_drift) {
        sprintf("terminal_drift[t][tid][:, %d]", j0)
      } else {
        "0.0"
      }
      a(paste0(I4, sprintf(
        "tdrift_vec = pt.set_subtensor(tdrift_vec[:, %d], pt.switch(pt.eq(miss_pt[:, %d], 0), y_pt[:, %d] - %s, %s))",
        j0, j0, j0, lmod_base_fn(j0), missing_expr)))
    }
    a(paste0(I4, "lp_mv = mvn_chol_logp(tdrift_vec, L_cov_obs)"))
    a(paste0(I4, "obs_lp = lp_mv"))
  } else if (has_normal && repeated) {
    normal_idx <- which(distributions == "normal")
    not_normal_idx <- which(distributions != "normal")

    a(paste0(I4, "residuals = residual_z.T"))
    for (j in normal_idx) {
      j0 <- j - 1
      tdrift_expr <- sprintf("tdrift_trees[t][tid][:, %d]", j0)
      a(paste0(I4, sprintf(
        "residuals = pt.set_subtensor(residuals[:, %d], pt.switch(pt.eq(miss_pt[:, %d], 0), y_pt[:, %d] - (%s + %s), residual_z[%d]))",
        j0, j0, j0, lmod_base_fn(j0), tdrift_expr, j0)))
    }
    if (!is.null(measurement_error)) {
      a(paste0(I4,
        "res_cov = diag_mat(se[0]) + pt.dot(diag_mat(sigma_residual), pt.dot(pt.dot(L_residual, L_residual.T), diag_mat(sigma_residual)))"))
      a(paste0(I4, "L_cov_res = pt.linalg.cholesky(res_cov)"))
    } else {
      a(paste0(I4,
        "L_cov_res = pt.dot(diag_mat(sigma_residual), L_residual)"))
    }
    a(paste0(I4, "lp_mv = mvn_chol_logp(residuals, L_cov_res)"))
    a(paste0(I4, "obs_lp = lp_mv"))
  } else {
    a(paste0(I4, "obs_lp = pt.zeros(N_obs)"))
  }

  # --- Non-normal variable likelihoods (vectorised) ---
  for (j in seq_along(distributions)) {
    j0 <- j - 1
    d <- distributions[j]
    if (d == "normal") next

    lmod <- lmod_base_fn(j0)
    if (!repeated && has_normal) {
      lmod <- sprintf("(%s + tdrift_vec[:, %d])", lmod, j0)
    } else if (tdrift && !repeated) {
      lmod <- sprintf("(%s + tdrift_trees[t][tid][:, %d])", lmod, j0)
    } else if (repeated && !has_normal) {
      lmod <- sprintf("(%s + residual_v[:, %d])", lmod, j0)
    } else if (repeated && has_normal) {
      lmod <- sprintf("(%s + residuals[:, %d])", lmod, j0)
    } else {
      lmod <- sprintf("(%s + tdrift_trees[t][tid][:, %d])", lmod, j0)
    }

    ll_expr <- switch(d,
      "bernoulli_logit" =
        sprintf("pm.logp(pm.Bernoulli.dist(logit_p=%s), y_pt[:, %d])", lmod, j0),
      "ordered_logistic" =
        sprintf("pm.logp(pm.OrderedLogistic.dist(eta=%s, cutpoints=c%d), pt.cast(y_pt[:, %d], 'int32') - 1)",
                lmod, j, j0),
      "poisson_softplus" =
        sprintf("pm.logp(pm.Poisson.dist(mu=obs_means[%d] * pt.softplus(%s)), pt.cast(y_pt[:, %d], 'int32'))",
                j0, lmod, j0),
      "negative_binomial_softplus" =
        sprintf("pm.logp(pm.NegativeBinomial.dist(mu=obs_means[%d] * pt.softplus(%s), alpha=phi%d), pt.cast(y_pt[:, %d], 'int32'))",
                j0, lmod, j, j0),
      "gamma_log" =
        sprintf("pm.logp(pm.Gamma.dist(alpha=shape%d, beta=shape%d / pt.exp(%s)), y_pt[:, %d])",
                j, j, lmod, j0),
      stop2("Unsupported distribution for PyMC: ", d)
    )

    a(paste0(I4, sprintf(
      "obs_lp = obs_lp + pt.switch(pt.eq(miss_pt[:, %d], 0), %s, 0.0)",
      j0, ll_expr)))
  }

  a(paste0(I4, "tree_lps.append(obs_lp)"))
  a(paste0(I3,
    "all_tree_lps = pt.stack(tree_lps, axis=0)"))
  a(paste0(I3,
    "pm.Potential('log_lik', pt.sum(pt.logsumexp(all_tree_lps, axis=0)))"))
  a("")
  a(paste0(I2, "return model"))

  lines
}
