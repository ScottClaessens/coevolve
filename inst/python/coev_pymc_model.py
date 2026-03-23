"""PyMC model builder for the coevolve dynamic coevolutionary model.

Constructs the PyMC model from a data_dict containing numeric arrays (from
standata_to_pymc) and model config flags (from coev_make_pymc).
"""

import numpy as np
import pymc as pm
import pytensor.tensor as pt


# ---------------------------------------------------------------------------
# Standalone helpers
# ---------------------------------------------------------------------------


def diag_mat(v):
    return v[:, None] * pt.eye(v.shape[0])


def mvn_chol_logp(value, chol):
    z = pt.linalg.solve_triangular(chol, value[..., None], lower=True)[..., 0]
    logdet = pt.sum(pt.log(pt.diagonal(chol, axis1=-2, axis2=-1)), axis=-1)
    quad = pt.sum(z**2, axis=-1)
    jdim = pt.cast(value.shape[-1], value.dtype)
    return -0.5 * (jdim * np.log(2.0 * np.pi) + quad) - logdet


def ksolve(A, Q, J):
    """Solve AX + XA^T = -Q for X (continuous Lyapunov equation)."""
    _I = pt.eye(J)
    M = (
        A[:, :, None, None] * _I[None, None, :, :]
        + _I[:, :, None, None] * A[None, None, :, :]
    )
    M = M.transpose(0, 2, 1, 3).reshape((J * J, J * J))
    X = pt.linalg.solve(M, -Q.reshape((-1,))).reshape((J, J))
    return 0.5 * (X + X.T)


def sample_lkj_cholesky(name, n, eta):
    """LKJ Cholesky via stick-breaking parameterization."""
    if n == 1:
        return pt.eye(1)

    n_corr = n * (n - 1) // 2
    raw = pm.Flat(f"_{name}_raw", shape=n_corr)
    z = pt.tanh(raw)

    entries = {(0, 0): pt.constant(1.0)}
    idx = 0
    for i in range(1, n):
        for j in range(i + 1):
            if j == 0:
                entries[(i, j)] = z[idx]
                idx += 1
            elif j == i:
                sq_sum = sum(entries[(i, k)] ** 2 for k in range(j))
                entries[(i, j)] = pt.sqrt(pt.maximum(1.0 - sq_sum, 1e-10))
            else:
                sq_sum = sum(entries[(i, k)] ** 2 for k in range(j))
                entries[(i, j)] = z[idx] * pt.sqrt(pt.maximum(1.0 - sq_sum, 1e-10))
                idx += 1

    rows = [[entries.get((i, j), pt.constant(0.0)) for j in range(n)] for i in range(n)]
    L = pt.stacklists(rows)

    lkj_logp = sum(
        (n - i - 1 + 2 * (eta - 1)) * pt.log(entries[(i, i)]) for i in range(1, n)
    )
    jac = pt.sum(pt.log(1.0 - z**2 + 1e-10))
    pm.Potential(f"_{name}_lkj_prior", lkj_logp + jac)
    return L


# ---------------------------------------------------------------------------
# Data helpers
# ---------------------------------------------------------------------------


def _as_int_list(val):
    """Convert a scalar or array-like to a Python list of ints (empty-safe)."""
    arr = np.atleast_1d(np.asarray(val, dtype=np.int32)).ravel()
    if arr.size == 0:
        return []
    return arr.tolist()


def _flag(data, key):
    """Read an int flag from data_dict as a Python bool."""
    return bool(int(data[key]))


# ---------------------------------------------------------------------------
# CoevPymcModel
# ---------------------------------------------------------------------------


class CoevPymcModel:
    """Build a PyMC coevolve model from a data_dict."""

    def build(self, data, compile_mode="cpu"):
        with pm.Model() as model:
            J = int(data["J"])
            params = self._add_priors(data, J)
            trans = self._build_tp(data, J, params)
            self._add_likelihood(data, J, params, trans)
        return model

    def _make_rv(self, name, spec, shape=None, initval=None):
        """Create a PyMC RV from a prior spec dict with keys dist, args, constraint."""
        dist = str(spec["dist"])
        args = [float(a) for a in (spec.get("args") or [])]
        constraint = str(spec.get("constraint") or "none")
        kw = {}
        if shape is not None:
            kw["shape"] = shape
        if initval is not None:
            kw["initval"] = initval

        if dist == "std_normal":
            if constraint == "upper_zero":
                return pm.TruncatedNormal(name, mu=0.0, sigma=1.0, upper=0.0, **kw)
            elif constraint == "lower_zero":
                return pm.HalfNormal(name, sigma=1.0, **kw)
            return pm.Normal(name, mu=0.0, sigma=1.0, **kw)

        if dist == "normal":
            mu, sigma = args[0], args[1]
            if constraint == "upper_zero":
                return pm.TruncatedNormal(name, mu=mu, sigma=sigma, upper=0.0, **kw)
            elif constraint == "lower_zero":
                return pm.TruncatedNormal(name, mu=mu, sigma=sigma, lower=0.0, **kw)
            return pm.Normal(name, mu=mu, sigma=sigma, **kw)

        if dist == "exponential":
            return pm.Exponential(name, lam=args[0], **kw)

        if dist == "gamma":
            return pm.Gamma(name, alpha=args[0], beta=args[1], **kw)

        raise ValueError(f"Unknown prior distribution: {dist!r}")

    def _add_priors(self, data, J):
        specs = data["prior_specs"]
        N_tree = int(data["N_tree"])
        N_tips = int(data["N_tips"])
        N_seg = int(data["N_seg"])

        tdrift = _flag(data, "tdrift")
        repeated = _flag(data, "repeated")
        needs_terminal_drift = _flag(data, "needs_terminal_drift")
        estimate_correlated_drift = _flag(data, "estimate_correlated_drift")
        has_dist_mat = _flag(data, "has_dist_mat")
        n_offdiag = int(data["n_offdiag"])
        lkj_eta_drift = float(data["lkj_eta_drift"])
        lkj_eta_residual = float(data["lkj_eta_residual"])
        nb_j0 = [int(x) for x in np.atleast_1d(data["nb_j0"])]
        gamma_j0 = [int(x) for x in np.atleast_1d(data["gamma_j0"])]
        ordered_j = [int(x) for x in np.atleast_1d(data["ordered_j"])]
        ordered_ncuts = [int(x) for x in np.atleast_1d(data["ordered_ncuts"])]
        nonnormal_j0 = [int(x) for x in np.atleast_1d(data["nonnormal_j0"])]
        normal_j0 = [int(x) for x in np.atleast_1d(data["normal_j0"])]
        prior_only = _flag(data, "prior_only")
        inv_overdisp = np.asarray(data["inv_overdisp"], dtype=np.float64)

        tip_np = np.asarray(data["tip"], dtype=np.int32)
        internal_mask_np = (tip_np[:, 1:] == 0).astype(np.int32)
        N_internal = int(np.max(np.sum(internal_mask_np, axis=1))) if N_seg > 1 else 0

        p = {}

        p["A_diag"] = self._make_rv(
            "A_diag", specs["A_diag"], shape=J,
            initval=-0.5 * np.ones(J, dtype=np.float64),
        )

        if n_offdiag > 0:
            p["A_offdiag"] = self._make_rv("A_offdiag", specs["A_offdiag"], shape=n_offdiag)

        if estimate_correlated_drift:
            p["L_R"] = sample_lkj_cholesky("L_R", J, lkj_eta_drift)

        p["Q_sigma"] = self._make_rv(
            "Q_sigma", specs["Q_sigma"], shape=J,
            initval=0.5 * np.ones(J, dtype=np.float64),
        )
        p["b"] = self._make_rv("b", specs["b"], shape=J)
        p["eta_anc"] = self._make_rv("eta_anc", specs["eta_anc"], shape=(N_tree, J))
        p["z_drift"] = pm.Normal("z_drift", mu=0.0, sigma=1.0, shape=(N_tree, N_internal, J))

        if needs_terminal_drift:
            p["terminal_drift"] = pm.Normal(
                "terminal_drift", mu=0.0, sigma=1.0, shape=(N_tree, N_tips, J)
            )
            # Cancel the implicit N(0,1) prior for slots that Stan omits.
            has_normal = len(normal_j0) > 0
            if has_normal and not repeated:
                for j0 in nonnormal_j0:
                    pm.Potential(
                        f"_td_nonnorm_cancel_{j0}",
                        pt.sum(0.5 * p["terminal_drift"][:, :, j0] ** 2),
                    )
                for j0 in normal_j0:
                    if prior_only:
                        pm.Potential(
                            f"_td_norm_cancel_{j0}",
                            pt.sum(0.5 * p["terminal_drift"][:, :, j0] ** 2),
                        )

        for j1, nc in zip(ordered_j, ordered_ncuts):
            c_raw = self._make_rv(
                f"c{j1}_raw", specs["c"], shape=nc,
                initval=np.linspace(-1.0, 1.0, nc, dtype=np.float64),
            )
            p[f"c{j1}_raw"] = c_raw
            p[f"c{j1}"] = pm.Deterministic(f"c{j1}", pt.sort(c_raw))

        for j0 in nb_j0:
            j1 = j0 + 1
            phi_spec = specs.get("phi")
            if phi_spec is None:
                p[f"phi{j1}"] = pm.TruncatedNormal(
                    f"phi{j1}", mu=inv_overdisp[j0], sigma=inv_overdisp[j0], lower=0.0,
                )
            else:
                p[f"phi{j1}"] = self._make_rv(f"phi{j1}", phi_spec)

        for j0 in gamma_j0:
            p[f"shape{j0 + 1}"] = self._make_rv(f"shape{j0 + 1}", specs["shape"])

        if has_dist_mat:
            p["dist_z"] = pm.Normal("dist_z", mu=0.0, sigma=1.0, shape=(N_tips, J))
            p["sigma_dist"] = self._make_rv("sigma_dist", specs["sigma_dist"], shape=J)
            p["rho_dist"] = self._make_rv("rho_dist", specs["rho_dist"], shape=J)

        if repeated:
            N_obs = int(data["N_obs"])
            p["residual_z"] = pm.Normal("residual_z", mu=0.0, sigma=1.0, shape=(J, N_obs))
            p["sigma_residual"] = self._make_rv("sigma_residual", specs["sigma_residual"], shape=J)
            p["L_residual"] = sample_lkj_cholesky("L_residual", J, lkj_eta_residual)

        return p

    def _build_tp(self, data, J, p):
        """Build A matrix, drift covariance, branch caches, and tree traversal."""
        effects_mat = data["effects_mat"]
        n_offdiag = int(data["n_offdiag"])
        tdrift = _flag(data, "tdrift")
        repeated = _flag(data, "repeated")
        residual_v_flag = _flag(data, "residual_v")
        estimate_correlated_drift = _flag(data, "estimate_correlated_drift")
        has_dist_mat = _flag(data, "has_dist_mat")
        dist_cov_type = str(data.get("dist_cov_type") or "")
        has_measurement_error = _flag(data, "has_measurement_error")

        N_tree = int(data["N_tree"])
        N_seg = int(data["N_seg"])
        unique_lengths = np.asarray(data["unique_lengths"], dtype=np.float64)
        length_index = np.asarray(data["length_index"], dtype=np.int32)
        tip_to_seg = np.asarray(data["tip_to_seg"], dtype=np.int32)

        A_mat = diag_mat(p["A_diag"])
        if n_offdiag > 0:
            ticker = 0
            for i in range(J):
                for j in range(J):
                    if i != j and effects_mat[i, j] == 1:
                        A_mat = pt.set_subtensor(A_mat[i, j], p["A_offdiag"][ticker])
                        ticker += 1
        pm.Deterministic("A", A_mat)

        if estimate_correlated_drift:
            Q = pt.dot(
                diag_mat(p["Q_sigma"]),
                pt.dot(pt.dot(p["L_R"], p["L_R"].T), diag_mat(p["Q_sigma"])),
            )
            pm.Deterministic("cor_R", pt.dot(p["L_R"], p["L_R"].T))
        else:
            Q = diag_mat(p["Q_sigma"] ** 2)
        pm.Deterministic("Q", Q)

        Q_inf = ksolve(A_mat, Q, J)

        # Matrix exponential via scaling-and-squaring: 12-term Taylor of
        # exp(A*dt/2^8), then 8 repeated squarings to recover exp(A*dt).
        eye_J = pt.eye(J)
        _A_dt = A_mat[None, :, :] * (unique_lengths[:, None, None] / 256.0)
        _eye_b = pt.broadcast_to(eye_J[None, :, :], _A_dt.shape)
        _S = _eye_b
        _T = _eye_b
        for _k in range(1, 13):
            _T = pt.matmul(_T, _A_dt) * (1.0 / _k)
            _S = _S + _T
        for _ in range(8):
            _S = pt.matmul(_S, _S)
        A_delta_cache = _S

        _Qi = Q_inf[None, :, :]
        _V = _Qi - pt.matmul(A_delta_cache, pt.matmul(_Qi, A_delta_cache.transpose(0, 2, 1)))
        VCV_cache = 0.5 * (_V + _V.transpose(0, 2, 1))
        L_VCV_cache = pt.linalg.cholesky(VCV_cache)

        _A_inv = pt.linalg.solve(A_mat, eye_J)
        _As = pt.matmul(_A_inv[None, :, :], A_delta_cache - _eye_b)
        A_solve_cache = 0.5 * (_As + _As.transpose(0, 2, 1))
        b_delta_cache = pt.matmul(A_solve_cache, p["b"][None, :, None])[:, :, 0]

        n_levels = int(data["n_levels"])
        max_level_size = int(data["max_level_size"])
        root_ids = np.atleast_1d(np.asarray(data["root_ids"], dtype=np.int32))
        _lshape = (N_tree, n_levels, max_level_size)
        level_node_ids = np.asarray(data["level_node_ids"], dtype=np.int32).reshape(_lshape)
        level_parent_ids = np.asarray(data["level_parent_ids"], dtype=np.int32).reshape(_lshape)
        level_length_idx = np.asarray(data["level_length_idx"], dtype=np.int32).reshape(_lshape)
        level_is_internal = np.asarray(data["level_is_internal"], dtype=np.int32).reshape(_lshape)
        level_drift_idx = np.asarray(data["level_drift_idx"], dtype=np.int32).reshape(_lshape)

        eta_trees = []
        tip_L_VCV_trees = []
        for t in range(N_tree):
            eta = pt.zeros((N_seg, J))
            eta = pt.set_subtensor(eta[int(root_ids[t])], p["eta_anc"][t])
            for _l in range(n_levels):
                _nids = level_node_ids[t, _l]
                _pids = level_parent_ids[t, _l]
                _li = level_length_idx[t, _l]
                _is_int = pt.as_tensor_variable(level_is_internal[t, _l].astype(np.float64))
                _didx = level_drift_idx[t, _l]
                _base = pt.matmul(A_delta_cache[_li], eta[_pids][:, :, None])[:, :, 0] + b_delta_cache[_li]
                _noise = pt.matmul(L_VCV_cache[_li], p["z_drift"][t, _didx][:, :, None])[:, :, 0]
                eta = pt.set_subtensor(eta[_nids], _base + _is_int[:, None] * _noise)
            eta = pt.set_subtensor(eta[int(root_ids[t])], p["eta_anc"][t])
            eta_trees.append(eta)
            _tip_li = length_index[t][tip_to_seg[t]]
            tip_L_VCV_trees.append(L_VCV_cache[_tip_li])

        trans = {
            "eta_trees": eta_trees,
            "tip_L_VCV_trees": tip_L_VCV_trees,
            "A_delta_cache": A_delta_cache,
            "L_VCV_cache": L_VCV_cache,
        }

        if tdrift:
            trans["tdrift_trees"] = [
                pt.matmul(tip_L_VCV_trees[t], p["terminal_drift"][t][:, :, None])[:, :, 0]
                for t in range(N_tree)
            ]

        if residual_v_flag:
            trans["residual_v"] = pt.dot(
                pt.dot(diag_mat(p["sigma_residual"]), p["L_residual"]), p["residual_z"]
            ).T

        if repeated:
            pm.Deterministic("cor_residual", pt.dot(p["L_residual"], p["L_residual"].T))

        if has_dist_mat:
            dist_mat_arr = np.array(data["dist_mat"], dtype=np.float64)
            _dm = pt.as_tensor_variable(dist_mat_arr)
            dist_v_list = []
            for j_gp in range(J):
                if dist_cov_type == "exp_quad":
                    dc = p["sigma_dist"][j_gp] * pt.exp(-(_dm**2) / p["rho_dist"][j_gp])
                elif dist_cov_type == "exponential":
                    dc = p["sigma_dist"][j_gp] * pt.exp(-_dm / p["rho_dist"][j_gp])
                elif dist_cov_type == "matern32":
                    dc = (
                        p["sigma_dist"][j_gp]
                        * (1.0 + pt.sqrt(3.0) * _dm / p["rho_dist"][j_gp])
                        * pt.exp(-pt.sqrt(3.0) * _dm / p["rho_dist"][j_gp])
                    )
                else:
                    raise ValueError(f"Unknown dist_cov type: {dist_cov_type!r}")
                dc = pt.fill_diagonal(dc, p["sigma_dist"][j_gp] + 0.01)
                dist_v_list.append(pt.dot(pt.linalg.cholesky(dc), p["dist_z"][:, j_gp]))
            trans["dist_v"] = pt.stack(dist_v_list, axis=1)

        return trans

    def _add_likelihood(self, data, J, p, trans):
        if _flag(data, "prior_only"):
            return

        distributions = list(data["distributions"])
        repeated = _flag(data, "repeated")
        tdrift = _flag(data, "tdrift")
        residual_v_flag = _flag(data, "residual_v")
        has_dist_mat = _flag(data, "has_dist_mat")
        has_measurement_error = _flag(data, "has_measurement_error")
        needs_terminal_drift = _flag(data, "needs_terminal_drift")

        N_obs = int(data["N_obs"])
        N_tree = int(data["N_tree"])
        y_pt = pt.as_tensor_variable(np.array(data["y"], dtype=np.float64))
        miss_pt = pt.as_tensor_variable(np.array(data["miss"], dtype=np.int32))
        tid = np.array(data["tip_id"], dtype=np.int32)
        obs_means = np.asarray(data["obs_means"], dtype=np.float64)

        has_normal = "normal" in distributions
        normal_idx = [j for j, d in enumerate(distributions) if d == "normal"]
        nonnormal_j0 = [int(x) for x in np.atleast_1d(data["nonnormal_j0"])]

        def base_lmod(j0, eta_obs):
            lm = eta_obs[:, j0]
            if has_dist_mat:
                lm = lm + trans["dist_v"][tid, j0]
            return lm

        se = np.array(data["se"], dtype=np.float64) if has_measurement_error else None

        tree_lps = []
        for t in range(N_tree):
            eta_obs = trans["eta_trees"][t][tid]
            tdrift_vec = None
            residuals = None

            if has_normal and not repeated:
                L_cov_raw = trans["tip_L_VCV_trees"][t][tid]
                if has_measurement_error:
                    VCV_obs = pt.matmul(L_cov_raw, L_cov_raw.transpose(0, 2, 1))
                    se_diag = pt.as_tensor_variable(se)[:, :, None] * pt.eye(J)[None, :, :]
                    L_cov_obs = pt.linalg.cholesky(VCV_obs + se_diag)
                else:
                    L_cov_obs = L_cov_raw

                if needs_terminal_drift:
                    tdrift_vec = p["terminal_drift"][t][tid]  # shape (N_obs, J)
                else:
                    tdrift_vec = pt.zeros((N_obs, J))

                for j in normal_idx:
                    j0 = j
                    missing_fill = p["terminal_drift"][t][tid][:, j0] if needs_terminal_drift else pt.constant(np.float64(0.0))
                    tdrift_vec = pt.set_subtensor(
                        tdrift_vec[:, j0],
                        pt.switch(
                            pt.eq(miss_pt[:, j0], 0),
                            y_pt[:, j0] - base_lmod(j0, eta_obs),
                            missing_fill,
                        ),
                    )
                obs_lp = mvn_chol_logp(tdrift_vec, L_cov_obs)

            elif has_normal and repeated:
                residuals = p["residual_z"].T  # (N_obs, J)
                for j in normal_idx:
                    j0 = j
                    tdrift_expr = trans["tdrift_trees"][t][tid][:, j0]
                    residuals = pt.set_subtensor(
                        residuals[:, j0],
                        pt.switch(
                            pt.eq(miss_pt[:, j0], 0),
                            y_pt[:, j0] - (base_lmod(j0, eta_obs) + tdrift_expr),
                            p["residual_z"][j0],
                        ),
                    )
                if has_measurement_error:
                    # Per-observation SE: se is (N_obs, J), build (N_obs, J, J)
                    # diag matrices, matching Stan's diag_matrix(se[i,]) per obs.
                    se_diag = pt.as_tensor_variable(se)[:, :, None] * pt.eye(J)[None, :, :]
                    base_cov = pt.dot(
                        diag_mat(p["sigma_residual"]),
                        pt.dot(pt.dot(p["L_residual"], p["L_residual"].T), diag_mat(p["sigma_residual"])),
                    )
                    L_cov_res = pt.linalg.cholesky(se_diag + base_cov[None, :, :])
                else:
                    L_cov_res = pt.dot(diag_mat(p["sigma_residual"]), p["L_residual"])
                obs_lp = mvn_chol_logp(residuals, L_cov_res)

            else:
                obs_lp = pt.zeros(N_obs)

            for j0, d in enumerate(distributions):
                if d == "normal":
                    continue
                j1 = j0 + 1

                base = base_lmod(j0, eta_obs)
                if not repeated and has_normal:
                    lmod = base + tdrift_vec[:, j0]
                elif tdrift and not repeated:
                    lmod = base + trans["tdrift_trees"][t][tid][:, j0]
                elif repeated and not has_normal:
                    lmod = base + trans["residual_v"][:, j0]
                elif repeated and has_normal:
                    lmod = base + residuals[:, j0]
                else:
                    raise RuntimeError(
                        f"Unreachable lmod branch: repeated={repeated}, "
                        f"has_normal={has_normal}, tdrift={tdrift}"
                    )

                # JAX traces both branches of pt.switch, so missing entries
                # (filled with -9999) must be replaced with valid values before
                # calling pm.logp to avoid NaN gradients.
                miss_j = pt.eq(miss_pt[:, j0], 1)
                y_obs = y_pt[:, j0]

                if d == "bernoulli_logit":
                    y_safe = pt.where(miss_j, pt.zeros_like(y_obs), y_obs)
                    ll = pm.logp(pm.Bernoulli.dist(logit_p=lmod), y_safe)
                elif d == "ordered_logistic":
                    y_safe = pt.where(miss_j, pt.ones_like(y_obs), y_obs)
                    ll = pm.logp(
                        pm.OrderedLogistic.dist(eta=lmod, cutpoints=p[f"c{j1}"]),
                        pt.cast(y_safe, "int32") - 1,
                    )
                elif d == "poisson_softplus":
                    y_safe = pt.where(miss_j, pt.zeros_like(y_obs), y_obs)
                    ll = pm.logp(
                        pm.Poisson.dist(mu=obs_means[j0] * pt.softplus(lmod)),
                        pt.cast(y_safe, "int32"),
                    )
                elif d == "negative_binomial_softplus":
                    y_safe = pt.where(miss_j, pt.zeros_like(y_obs), y_obs)
                    ll = pm.logp(
                        pm.NegativeBinomial.dist(
                            mu=obs_means[j0] * pt.softplus(lmod),
                            alpha=p[f"phi{j1}"],
                        ),
                        pt.cast(y_safe, "int32"),
                    )
                elif d == "gamma_log":
                    alpha = p[f"shape{j1}"]
                    y_safe = pt.where(miss_j, pt.ones_like(y_obs), y_obs)
                    ll = pm.logp(
                        pm.Gamma.dist(alpha=alpha, beta=alpha / pt.exp(lmod)),
                        y_safe,
                    )
                else:
                    raise ValueError(f"Unsupported distribution for PyMC: {d!r}")

                obs_lp = obs_lp + pt.switch(pt.eq(miss_pt[:, j0], 0), ll, 0.0)

            tree_lps.append(obs_lp)

        all_tree_lps = pt.stack(tree_lps, axis=0)
        pm.Potential("log_lik", pt.sum(pt.logsumexp(all_tree_lps, axis=0)))


def build_model(data, compile_mode="cpu"):
    """Build and return a PyMC coevolve model."""
    return CoevPymcModel().build(data, compile_mode)
