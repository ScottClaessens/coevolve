"""Pure JAX model for the coevolve dynamic coevolutionary model.

Implements the log-density as a pure JAX function with parameter transforms,
using numpyro.distributions for distribution log-probs. No PyMC or PyTensor.

Designed to be called via nutpie.compiled_pyfunc.from_pyfunc().
"""

import jax
import jax.numpy as jnp
import jax.scipy.linalg
import numpy as np
import numpyro.distributions as dist


# --------------------------------------------------------------------------
# Math helpers
# --------------------------------------------------------------------------


def ksolve(A, Q, J):
    """Solve AX + XA^T = -Q for X (continuous Lyapunov equation)."""
    I = jnp.eye(J)
    M = (
        A[:, :, None, None] * I[None, None, :, :]
        + I[:, :, None, None] * A[None, None, :, :]
    )
    M = M.transpose(0, 2, 1, 3).reshape((J * J, J * J))
    X = jnp.linalg.solve(M, -Q.reshape((-1,))).reshape((J, J))
    return 0.5 * (X + X.T)


def matrix_exp_batch(A_scaled, n_terms=12, n_squarings=8):
    """Matrix exponential via scaling-and-squaring (batched).

    Computes exp(A_scaled * 2^n_squarings) by Taylor-expanding exp(A_scaled)
    then squaring n_squarings times.
    """
    shape = A_scaled.shape
    eye = jnp.broadcast_to(jnp.eye(shape[-1])[None, :, :], shape)
    S = eye
    T = eye
    for k in range(1, n_terms + 1):
        T = jnp.matmul(T, A_scaled) * (1.0 / k)
        S = S + T
    for _ in range(n_squarings):
        S = jnp.matmul(S, S)
    return S


def mvn_chol_logp(value, chol):
    """Multivariate normal log-density given Cholesky factor (batched)."""
    z = jax.scipy.linalg.solve_triangular(
        chol, value[..., None], lower=True
    )[..., 0]
    logdet = jnp.sum(jnp.log(jnp.diagonal(chol, axis1=-2, axis2=-1)), axis=-1)
    quad = jnp.sum(z**2, axis=-1)
    J = value.shape[-1]
    return -0.5 * (J * jnp.log(2.0 * jnp.pi) + quad) - logdet


# --------------------------------------------------------------------------
# Parameter transforms (unconstrained -> constrained)
# --------------------------------------------------------------------------


def transform_upper_zero(raw):
    """Unconstrained -> upper=0: use -softplus."""
    return -jax.nn.softplus(raw)


def transform_upper_zero_logdet(raw):
    """Log abs det jacobian for upper=0 transform."""
    return jnp.sum(jax.nn.log_sigmoid(raw))


def transform_lower_zero(raw):
    """Unconstrained -> lower=0: use softplus."""
    return jax.nn.softplus(raw)


def transform_lower_zero_logdet(raw):
    """Log abs det jacobian for lower=0 transform."""
    return jnp.sum(jax.nn.log_sigmoid(raw))


def transform_ordered(raw):
    """Unconstrained -> ordered vector via cumulative softplus."""
    # first element is unconstrained, rest are positive increments
    increments = jax.nn.softplus(raw[1:])
    return jnp.concatenate([raw[:1], raw[:1] + jnp.cumsum(increments)])


def transform_ordered_logdet(raw):
    """Log abs det jacobian for ordered transform."""
    return jnp.sum(jax.nn.log_sigmoid(raw[1:]))


def raw_to_cholesky(raw_vec, n):
    """Convert unconstrained vector to Cholesky factor via stick-breaking.

    Returns (L, log_det_jacobian) where L is lower-triangular with positive
    diagonal, and log_det_jacobian accounts for the tanh transform.
    """
    if n == 1:
        return jnp.eye(1), 0.0

    z = jnp.tanh(raw_vec)
    jac = jnp.sum(jnp.log(1.0 - z**2 + 1e-10))

    L = jnp.zeros((n, n))
    idx = 0
    for i in range(n):
        if i == 0:
            L = L.at[0, 0].set(1.0)
        else:
            cum_sq = 0.0
            for j in range(i):
                remainder = jnp.sqrt(jnp.maximum(1.0 - cum_sq, 1e-10))
                L = L.at[i, j].set(z[idx] * remainder)
                cum_sq = cum_sq + L[i, j] ** 2
                idx += 1
            L = L.at[i, i].set(jnp.sqrt(jnp.maximum(1.0 - cum_sq, 1e-10)))

    return L, jac


def lkj_cholesky_logp(L, eta, n):
    """LKJ correlation Cholesky log-density (unnormalized)."""
    if n == 1:
        return 0.0
    logp = 0.0
    for i in range(1, n):
        logp = logp + (n - i - 1 + 2 * (eta - 1)) * jnp.log(L[i, i])
    return logp


# --------------------------------------------------------------------------
# Prior log-densities
# --------------------------------------------------------------------------


def prior_logp(value, spec):
    """Evaluate prior log-density from a spec dict.

    spec has keys: dist, args, constraint.
    This evaluates the density of the CONSTRAINED value under the named prior.
    """
    d = spec["dist"]
    args = [float(a) for a in (spec.get("args") or [])]

    if d == "std_normal":
        return dist.Normal(0.0, 1.0).log_prob(value).sum()
    elif d == "normal":
        return dist.Normal(args[0], args[1]).log_prob(value).sum()
    elif d == "exponential":
        return dist.Exponential(args[0]).log_prob(value).sum()
    elif d == "gamma":
        return dist.Gamma(args[0], args[1]).log_prob(value).sum()
    else:
        raise ValueError(f"Unknown prior distribution: {d!r}")


# --------------------------------------------------------------------------
# Model builder
# --------------------------------------------------------------------------


class CoevJaxModel:
    """Build a pure JAX log-density function from a data_dict."""

    def build(self, data):
        """Parse data dict and return (log_density_fn, param_info).

        param_info is a list of (name, shape, transform_type) tuples describing
        the unconstrained parameter layout.
        """
        self.data = data
        self.J = int(data["J"])
        self.N_tree = int(data["N_tree"])
        self.N_tips = int(data["N_tips"])
        self.N_seg = int(data["N_seg"])
        self.N_obs = int(data["N_obs"])

        self.distributions = list(data["distributions"])
        self.prior_specs = data["prior_specs"]
        self.effects_mat = np.asarray(data["effects_mat"], dtype=np.int32)
        self.n_offdiag = int(data["n_offdiag"])
        self.prior_only = bool(int(data["prior_only"]))

        self.tdrift = bool(int(data["tdrift"]))
        self.repeated = bool(int(data["repeated"]))
        self.needs_terminal_drift = bool(int(data["needs_terminal_drift"]))
        self.residual_v_flag = bool(int(data["residual_v"]))
        self.estimate_correlated_drift = bool(
            int(data["estimate_correlated_drift"])
        )
        self.has_dist_mat = bool(int(data["has_dist_mat"]))
        self.has_measurement_error = bool(int(data["has_measurement_error"]))

        self.ordered_j = [int(x) for x in np.atleast_1d(data["ordered_j"])]
        self.ordered_ncuts = [
            int(x) for x in np.atleast_1d(data["ordered_ncuts"])
        ]
        self.nb_j0 = [int(x) for x in np.atleast_1d(data["nb_j0"])]
        self.gamma_j0 = [int(x) for x in np.atleast_1d(data["gamma_j0"])]
        self.nonnormal_j0 = [
            int(x) for x in np.atleast_1d(data["nonnormal_j0"])
        ]
        self.normal_j0 = [int(x) for x in np.atleast_1d(data["normal_j0"])]

        self.lkj_eta_drift = float(data["lkj_eta_drift"])
        self.lkj_eta_residual = float(data["lkj_eta_residual"])
        self.inv_overdisp = np.asarray(data["inv_overdisp"], dtype=np.float64)
        self.obs_means = np.asarray(data["obs_means"], dtype=np.float64)

        # Precompute data arrays as JAX arrays
        self.y = jnp.array(data["y"], dtype=jnp.float64)
        self.miss = jnp.array(data["miss"], dtype=jnp.int32)
        self.tip_id = np.asarray(data["tip_id"], dtype=np.int32)
        self.unique_lengths = jnp.array(
            data["unique_lengths"], dtype=jnp.float64
        )
        self.length_index = np.asarray(data["length_index"], dtype=np.int32)
        self.tip_to_seg = np.asarray(data["tip_to_seg"], dtype=np.int32)

        if self.has_measurement_error:
            self.se = jnp.array(data["se"], dtype=jnp.float64)

        # Tree traversal level data
        self.n_levels = int(data["n_levels"])
        self.max_level_size = int(data["max_level_size"])
        self.root_ids = np.atleast_1d(
            np.asarray(data["root_ids"], dtype=np.int32)
        )
        _lshape = (self.N_tree, self.n_levels, self.max_level_size)
        self.level_node_ids = np.asarray(
            data["level_node_ids"], dtype=np.int32
        ).reshape(_lshape)
        self.level_parent_ids = np.asarray(
            data["level_parent_ids"], dtype=np.int32
        ).reshape(_lshape)
        self.level_length_idx = np.asarray(
            data["level_length_idx"], dtype=np.int32
        ).reshape(_lshape)
        self.level_is_internal = np.asarray(
            data["level_is_internal"], dtype=np.int32
        ).reshape(_lshape)
        self.level_drift_idx = np.asarray(
            data["level_drift_idx"], dtype=np.int32
        ).reshape(_lshape)

        # Compute N_internal for z_drift shape
        tip_np = np.asarray(data["tip"], dtype=np.int32)
        internal_mask_np = (tip_np[:, 1:] == 0).astype(np.int32)
        self.N_internal = (
            int(np.max(np.sum(internal_mask_np, axis=1)))
            if self.N_seg > 1
            else 0
        )

        # Build parameter layout
        self.param_info = self._build_param_layout()
        self.ndim = sum(int(np.prod(shape)) for _, shape, _ in self.param_info)

        # Build expand info for nutpie
        self.expand_info = self._build_expand_info()

        return self

    def _build_param_layout(self):
        """Define the unconstrained parameter vector layout.

        Returns list of (name, shape, transform_type) where transform_type
        is one of: 'none', 'upper_zero', 'lower_zero', 'ordered', 'cholesky'.
        """
        J = self.J
        info = []

        # A_diag: J params, upper=0
        info.append(("A_diag", (J,), "upper_zero"))

        # A_offdiag: n_offdiag params, unconstrained
        if self.n_offdiag > 0:
            info.append(("A_offdiag", (self.n_offdiag,), "none"))

        # L_R: cholesky factor for drift correlation
        if self.estimate_correlated_drift:
            n_corr = J * (J - 1) // 2
            if n_corr > 0:
                info.append(("_L_R_raw", (n_corr,), "none"))

        # Q_sigma: J params, lower=0
        info.append(("Q_sigma", (J,), "lower_zero"))

        # b: J params, unconstrained
        info.append(("b", (J,), "none"))

        # eta_anc: (N_tree, J), unconstrained
        info.append(("eta_anc", (self.N_tree, J), "none"))

        # z_drift: (N_tree, N_internal, J), unconstrained
        info.append(("z_drift", (self.N_tree, self.N_internal, J), "none"))

        # terminal_drift: (N_tree, N_tips, J), unconstrained
        if self.needs_terminal_drift:
            info.append(
                ("terminal_drift", (self.N_tree, self.N_tips, J), "none")
            )

        # ordered cutpoints
        for j1, nc in zip(self.ordered_j, self.ordered_ncuts):
            info.append((f"c{j1}", (nc,), "ordered"))

        # neg binom overdispersion
        for j0 in self.nb_j0:
            info.append((f"phi{j0 + 1}", (1,), "lower_zero"))

        # gamma shape
        for j0 in self.gamma_j0:
            info.append((f"shape{j0 + 1}", (1,), "lower_zero"))

        # GP parameters
        if self.has_dist_mat:
            info.append(("dist_z", (self.N_tips, J), "none"))
            info.append(("sigma_dist", (J,), "lower_zero"))
            info.append(("rho_dist", (J,), "lower_zero"))

        # Residual parameters
        if self.repeated:
            info.append(("residual_z", (J, self.N_obs), "none"))
            info.append(("sigma_residual", (J,), "lower_zero"))
            n_corr_res = J * (J - 1) // 2
            if n_corr_res > 0:
                info.append(("_L_residual_raw", (n_corr_res,), "none"))

        return info

    def _expand_var_list(self):
        """Ordered list of (name, shape) for expanded output variables."""
        J = self.J
        vs = []
        vs.append(("A", (J, J)))
        vs.append(("Q", (J, J)))
        vs.append(("A_diag", (J,)))
        if self.n_offdiag > 0:
            vs.append(("A_offdiag", (self.n_offdiag,)))
        vs.append(("Q_sigma", (J,)))
        vs.append(("b", (J,)))
        vs.append(("eta_anc", (self.N_tree, J)))
        if self.estimate_correlated_drift:
            vs.append(("cor_R", (J, J)))
        for j1, nc in zip(self.ordered_j, self.ordered_ncuts):
            vs.append((f"c{j1}", (nc,)))
        for j0 in self.nb_j0:
            vs.append((f"phi{j0 + 1}", ()))
        for j0 in self.gamma_j0:
            vs.append((f"shape{j0 + 1}", ()))
        if self.repeated:
            vs.append(("sigma_residual", (J,)))
            vs.append(("cor_residual", (J, J)))
        return vs

    def _build_expand_info(self):
        """Build the expanded variable info for nutpie output."""
        vs = self._expand_var_list()
        names = [n for n, _ in vs]
        shapes = [s for _, s in vs]
        dtypes = [np.dtype("float64")] * len(vs)
        return names, shapes, dtypes

    def unpack_params(self, x):
        """Unpack flat unconstrained vector into named constrained params.

        Returns (params_dict, log_det_jacobian).
        """
        params = {}
        logdet = 0.0
        offset = 0

        for name, shape, transform in self.param_info:
            size = int(np.prod(shape))
            raw = x[offset : offset + size].reshape(shape)
            offset += size

            if transform == "none":
                params[name] = raw
            elif transform == "upper_zero":
                params[name] = transform_upper_zero(raw)
                logdet = logdet + transform_upper_zero_logdet(raw)
            elif transform == "lower_zero":
                params[name] = transform_lower_zero(raw)
                logdet = logdet + transform_lower_zero_logdet(raw)
            elif transform == "ordered":
                params[name] = transform_ordered(raw)
                logdet = logdet + transform_ordered_logdet(raw)
            elif transform == "cholesky":
                # Should not appear directly; cholesky is handled via _raw
                params[name] = raw

        return params, logdet

    def log_density(self, x):
        """Compute log-density at unconstrained parameter vector x.

        This is the function that gets JIT-compiled and differentiated by JAX.
        """
        params, logdet_jac = self.unpack_params(x)
        J = self.J
        lp = logdet_jac

        # --- Build A matrix ---
        A_diag = params["A_diag"]
        A_mat = jnp.diag(A_diag)
        if self.n_offdiag > 0:
            ticker = 0
            for i in range(J):
                for j in range(J):
                    if i != j and self.effects_mat[i, j] == 1:
                        A_mat = A_mat.at[i, j].set(
                            params["A_offdiag"][ticker]
                        )
                        ticker += 1

        # --- Build Q matrix ---
        if self.estimate_correlated_drift:
            n_corr = J * (J - 1) // 2
            if n_corr > 0:
                L_R, lkj_jac = raw_to_cholesky(params["_L_R_raw"], J)
                lp = lp + lkj_jac  # tanh jacobian
                lp = lp + lkj_cholesky_logp(L_R, self.lkj_eta_drift, J)
            else:
                L_R = jnp.eye(J)
            D = jnp.diag(params["Q_sigma"])
            Q = D @ (L_R @ L_R.T) @ D
        else:
            L_R = None
            Q = jnp.diag(params["Q_sigma"] ** 2)

        # --- Priors ---
        lp = lp + prior_logp(params["A_diag"], self.prior_specs["A_diag"])
        if self.n_offdiag > 0:
            lp = lp + prior_logp(
                params["A_offdiag"], self.prior_specs["A_offdiag"]
            )
        lp = lp + prior_logp(params["Q_sigma"], self.prior_specs["Q_sigma"])
        lp = lp + prior_logp(params["b"], self.prior_specs["b"])
        lp = lp + prior_logp(params["eta_anc"], self.prior_specs["eta_anc"])
        lp = lp + dist.Normal(0.0, 1.0).log_prob(params["z_drift"]).sum()

        if self.needs_terminal_drift:
            has_normal = len(self.normal_j0) > 0
            if has_normal and not self.repeated:
                # Only add std_normal prior for non-normal variables
                # (normal variables get their prior from the likelihood)
                for j0 in self.nonnormal_j0:
                    lp = lp + (
                        dist.Normal(0.0, 1.0)
                        .log_prob(params["terminal_drift"][:, :, j0])
                        .sum()
                    )
                if self.prior_only:
                    for j0 in self.normal_j0:
                        lp = lp + (
                            dist.Normal(0.0, 1.0)
                            .log_prob(params["terminal_drift"][:, :, j0])
                            .sum()
                        )
            else:
                lp = lp + (
                    dist.Normal(0.0, 1.0)
                    .log_prob(params["terminal_drift"])
                    .sum()
                )

        for j1, nc in zip(self.ordered_j, self.ordered_ncuts):
            lp = lp + prior_logp(params[f"c{j1}"], self.prior_specs["c"])

        for j0 in self.nb_j0:
            j1 = j0 + 1
            phi_spec = self.prior_specs.get("phi")
            if phi_spec is None:
                lp = lp + (
                    dist.TruncatedNormal(
                        self.inv_overdisp[j0], self.inv_overdisp[j0], low=0.0
                    )
                    .log_prob(params[f"phi{j1}"])
                    .sum()
                )
            else:
                lp = lp + prior_logp(params[f"phi{j1}"], phi_spec)

        for j0 in self.gamma_j0:
            lp = lp + prior_logp(
                params[f"shape{j0 + 1}"], self.prior_specs["shape"]
            )

        if self.has_dist_mat:
            lp = lp + dist.Normal(0.0, 1.0).log_prob(params["dist_z"]).sum()
            lp = lp + prior_logp(
                params["sigma_dist"], self.prior_specs["sigma_dist"]
            )
            lp = lp + prior_logp(
                params["rho_dist"], self.prior_specs["rho_dist"]
            )

        if self.repeated:
            lp = lp + (
                dist.Normal(0.0, 1.0).log_prob(params["residual_z"]).sum()
            )
            lp = lp + prior_logp(
                params["sigma_residual"], self.prior_specs["sigma_residual"]
            )
            n_corr_res = J * (J - 1) // 2
            if n_corr_res > 0:
                L_residual, lkj_jac_res = raw_to_cholesky(
                    params["_L_residual_raw"], J
                )
                lp = lp + lkj_jac_res
                lp = lp + lkj_cholesky_logp(
                    L_residual, self.lkj_eta_residual, J
                )
            else:
                L_residual = jnp.eye(J)

        # --- Transformed parameters ---
        Q_inf = ksolve(A_mat, Q, J)

        # Matrix exponential cache (batched over unique branch lengths)
        eye_J = jnp.eye(J)
        A_dt = A_mat[None, :, :] * (
            self.unique_lengths[:, None, None] / 256.0
        )
        A_delta_cache = matrix_exp_batch(A_dt)

        Qi = Q_inf[None, :, :]
        V = Qi - jnp.matmul(
            A_delta_cache, jnp.matmul(Qi, A_delta_cache.transpose(0, 2, 1))
        )
        VCV_cache = 0.5 * (V + V.transpose(0, 2, 1))
        L_VCV_cache = jnp.linalg.cholesky(VCV_cache)

        A_inv = jnp.linalg.solve(A_mat, eye_J)
        As = jnp.matmul(
            A_inv[None, :, :],
            A_delta_cache - jnp.broadcast_to(eye_J[None, :, :], A_delta_cache.shape),
        )
        A_solve_cache = 0.5 * (As + As.transpose(0, 2, 1))
        b_delta_cache = jnp.matmul(
            A_solve_cache, params["b"][None, :, None]
        )[:, :, 0]

        # --- Tree traversal ---
        eta_trees = []
        tip_L_VCV_trees = []
        for t in range(self.N_tree):
            eta = jnp.zeros((self.N_seg, J))
            eta = eta.at[int(self.root_ids[t])].set(params["eta_anc"][t])
            for _l in range(self.n_levels):
                _nids = self.level_node_ids[t, _l]
                _pids = self.level_parent_ids[t, _l]
                _li = self.level_length_idx[t, _l]
                _is_int = jnp.array(
                    self.level_is_internal[t, _l], dtype=jnp.float64
                )
                _didx = self.level_drift_idx[t, _l]
                _base = (
                    jnp.matmul(A_delta_cache[_li], eta[_pids][:, :, None])[
                        :, :, 0
                    ]
                    + b_delta_cache[_li]
                )
                _noise = jnp.matmul(
                    L_VCV_cache[_li],
                    params["z_drift"][t, _didx][:, :, None],
                )[:, :, 0]
                eta = eta.at[_nids].set(_base + _is_int[:, None] * _noise)
            eta = eta.at[int(self.root_ids[t])].set(params["eta_anc"][t])
            eta_trees.append(eta)
            _tip_li = self.length_index[t][self.tip_to_seg[t]]
            tip_L_VCV_trees.append(L_VCV_cache[_tip_li])

        # Compute tdrift if needed
        tdrift_trees = None
        if self.tdrift and self.needs_terminal_drift:
            tdrift_trees = []
            for t in range(self.N_tree):
                td = jnp.matmul(
                    tip_L_VCV_trees[t],
                    params["terminal_drift"][t][:, :, None],
                )[:, :, 0]
                tdrift_trees.append(td)

        # Compute residual_v if needed
        residual_v = None
        if self.residual_v_flag and self.repeated:
            D_res = jnp.diag(params["sigma_residual"])
            residual_v = (D_res @ L_residual @ params["residual_z"]).T

        # --- Likelihood ---
        if not self.prior_only:
            lp = lp + self._likelihood(
                params,
                eta_trees,
                tip_L_VCV_trees,
                tdrift_trees,
                residual_v,
                L_residual if self.repeated else None,
            )

        return lp

    def _likelihood(
        self,
        params,
        eta_trees,
        tip_L_VCV_trees,
        tdrift_trees,
        residual_v,
        L_residual,
    ):
        """Compute the log-likelihood contribution."""
        J = self.J
        tid = self.tip_id
        has_normal = "normal" in self.distributions
        normal_idx = [
            j for j, d in enumerate(self.distributions) if d == "normal"
        ]

        tree_lps = []
        for t in range(self.N_tree):
            eta_obs = eta_trees[t][tid]
            obs_lp = jnp.zeros(self.N_obs)

            # Build base linear model accessor
            def base_lmod(j0):
                lm = eta_obs[:, j0]
                # GP spatial effects would be added here
                return lm

            if has_normal and not self.repeated:
                L_cov_obs = tip_L_VCV_trees[t][tid]
                if self.has_measurement_error:
                    VCV_obs = jnp.matmul(
                        L_cov_obs, L_cov_obs.transpose(0, 2, 1)
                    )
                    se_diag = (
                        self.se[:, :, None] * jnp.eye(J)[None, :, :]
                    )
                    L_cov_obs = jnp.linalg.cholesky(VCV_obs + se_diag)

                if self.needs_terminal_drift:
                    tdrift_vec = params["terminal_drift"][t][tid]
                else:
                    tdrift_vec = jnp.zeros((self.N_obs, J))

                for j in normal_idx:
                    missing_fill = (
                        params["terminal_drift"][t][tid][:, j]
                        if self.needs_terminal_drift
                        else 0.0
                    )
                    tdrift_vec = tdrift_vec.at[:, j].set(
                        jnp.where(
                            self.miss[:, j] == 0,
                            self.y[:, j] - base_lmod(j),
                            missing_fill,
                        )
                    )
                obs_lp = obs_lp + mvn_chol_logp(tdrift_vec, L_cov_obs)

            elif has_normal and self.repeated:
                residuals = params["residual_z"].T  # (N_obs, J)
                for j in normal_idx:
                    tdrift_expr = tdrift_trees[t][tid][:, j]
                    residuals = residuals.at[:, j].set(
                        jnp.where(
                            self.miss[:, j] == 0,
                            self.y[:, j]
                            - (base_lmod(j) + tdrift_expr),
                            params["residual_z"][j],
                        )
                    )
                if self.has_measurement_error:
                    se_diag = (
                        self.se[:, :, None] * jnp.eye(J)[None, :, :]
                    )
                    D_res = jnp.diag(params["sigma_residual"])
                    base_cov = D_res @ (L_residual @ L_residual.T) @ D_res
                    L_cov_res = jnp.linalg.cholesky(
                        se_diag + base_cov[None, :, :]
                    )
                else:
                    L_cov_res = jnp.diag(params["sigma_residual"]) @ L_residual
                obs_lp = obs_lp + mvn_chol_logp(residuals, L_cov_res)

            # Non-normal likelihoods
            for j0, d in enumerate(self.distributions):
                if d == "normal":
                    continue
                j1 = j0 + 1

                base = base_lmod(j0)
                if not self.repeated and has_normal:
                    lmod = base + tdrift_vec[:, j0]
                elif self.tdrift and not self.repeated:
                    lmod = base + tdrift_trees[t][tid][:, j0]
                elif self.repeated and not has_normal:
                    lmod = base + residual_v[:, j0]
                elif self.repeated and has_normal:
                    lmod = base + residuals[:, j0]
                else:
                    lmod = base

                miss_j = self.miss[:, j0] == 1
                y_obs = self.y[:, j0]

                if d == "bernoulli_logit":
                    y_safe = jnp.where(miss_j, 0.0, y_obs)
                    ll = dist.Bernoulli(logits=lmod).log_prob(y_safe)
                elif d == "ordered_logistic":
                    y_safe = jnp.where(miss_j, 1.0, y_obs)
                    cutpoints = params[f"c{j1}"]
                    # OrderedLogistic: P(Y=k) via cumulative logistic
                    y_int = (y_safe - 1).astype(jnp.int32)
                    ll = dist.OrderedLogistic(
                        predictor=lmod, cutpoints=cutpoints
                    ).log_prob(y_int)
                elif d == "poisson_softplus":
                    y_safe = jnp.where(miss_j, 0.0, y_obs)
                    mu = self.obs_means[j0] * jax.nn.softplus(lmod)
                    ll = dist.Poisson(rate=mu).log_prob(y_safe.astype(jnp.int32))
                elif d == "negative_binomial_softplus":
                    y_safe = jnp.where(miss_j, 0.0, y_obs)
                    mu = self.obs_means[j0] * jax.nn.softplus(lmod)
                    phi = params[f"phi{j1}"].squeeze()
                    # NegBinomial2: mu parameterization
                    ll = dist.NegativeBinomial2(
                        mean=mu, concentration=phi
                    ).log_prob(y_safe.astype(jnp.int32))
                elif d == "gamma_log":
                    y_safe = jnp.where(miss_j, 1.0, y_obs)
                    alpha = params[f"shape{j1}"].squeeze()
                    beta = alpha / jnp.exp(lmod)
                    ll = dist.Gamma(concentration=alpha, rate=beta).log_prob(
                        y_safe
                    )
                else:
                    raise ValueError(f"Unsupported distribution: {d!r}")

                obs_lp = obs_lp + jnp.where(self.miss[:, j0] == 0, ll, 0.0)

            tree_lps.append(obs_lp)

        all_tree_lps = jnp.stack(tree_lps, axis=0)
        return jnp.sum(jax.scipy.special.logsumexp(all_tree_lps, axis=0))

    def make_expand_fn(self, seed1=0, seed2=0, chain=0):
        """Return a function that maps flat unconstrained -> named params.

        Keys are returned in the exact order of _expand_var_list().
        The JAX computation is JIT-compiled; only the final
        jax->numpy conversion runs in Python.
        """
        J = self.J
        var_order = [n for n, _ in self._expand_var_list()]
        expand_shapes = {n: s for n, s in self._expand_var_list()}

        # JIT-compile the pure JAX part
        @jax.jit
        def _expand_jax(x):
            params, _ = self.unpack_params(x)
            vals = []

            # Build A matrix
            A_mat = jnp.diag(params["A_diag"])
            if self.n_offdiag > 0:
                ticker = 0
                for i in range(J):
                    for j in range(J):
                        if i != j and self.effects_mat[i, j] == 1:
                            A_mat = A_mat.at[i, j].set(
                                params["A_offdiag"][ticker]
                            )
                            ticker += 1
            vals.append(A_mat)  # A

            # Build Q matrix and cor_R
            if self.estimate_correlated_drift:
                n_corr = J * (J - 1) // 2
                if n_corr > 0:
                    L_R, _ = raw_to_cholesky(
                        params["_L_R_raw"], J
                    )
                else:
                    L_R = jnp.eye(J)
                D = jnp.diag(params["Q_sigma"])
                vals.append(D @ (L_R @ L_R.T) @ D)  # Q
            else:
                vals.append(jnp.diag(params["Q_sigma"] ** 2))  # Q

            vals.append(params["A_diag"])
            if self.n_offdiag > 0:
                vals.append(params["A_offdiag"])
            vals.append(params["Q_sigma"])
            vals.append(params["b"])
            vals.append(params["eta_anc"])

            if self.estimate_correlated_drift:
                vals.append(L_R @ L_R.T)  # cor_R

            for j1, nc in zip(self.ordered_j, self.ordered_ncuts):
                vals.append(params[f"c{j1}"])

            for j0 in self.nb_j0:
                vals.append(params[f"phi{j0 + 1}"].squeeze())

            for j0 in self.gamma_j0:
                vals.append(params[f"shape{j0 + 1}"].squeeze())

            if self.repeated:
                vals.append(params["sigma_residual"])
                n_corr_res = J * (J - 1) // 2
                if n_corr_res > 0:
                    L_res, _ = raw_to_cholesky(
                        params["_L_residual_raw"], J
                    )
                else:
                    L_res = jnp.eye(J)
                vals.append(L_res @ L_res.T)  # cor_residual

            return vals

        # Warm up JIT
        _x0 = jnp.zeros(self.ndim)
        _expand_jax(_x0)

        def expand(x):
            jax_vals = _expand_jax(jnp.asarray(x))
            result = {}
            for i, name in enumerate(var_order):
                result[name] = np.asarray(
                    jax_vals[i], dtype=np.float64
                )
            return result

        return expand


def build_nutpie_model(data):
    """Build a compiled nutpie model from data dict.

    Returns a nutpie CompiledModel ready for nutpie.sample().
    """
    import nutpie.compiled_pyfunc

    # Force JAX to CPU
    jax.config.update("jax_platforms", "cpu")
    jax.config.update("jax_enable_x64", True)

    model = CoevJaxModel()
    model.build(data)

    def make_logp_fn():
        logp_and_grad = jax.jit(jax.value_and_grad(model.log_density))
        # Warm up JIT (includes tracing + XLA compilation)
        _x0 = jnp.zeros(model.ndim)
        logp_and_grad(_x0)

        def logp(x):
            val, grad = logp_and_grad(x)
            # np.asarray on a JAX DeviceArray is a zero-copy
            # view when both are on CPU
            return val.item(), np.asarray(grad)

        return logp

    def make_expand_fn(seed1=0, seed2=0, chain=0):
        return model.make_expand_fn(seed1, seed2, chain)

    expand_names, expand_shapes, expand_dtypes = model.expand_info

    compiled = nutpie.compiled_pyfunc.from_pyfunc(
        ndim=model.ndim,
        make_logp_fn=make_logp_fn,
        make_expand_fn=make_expand_fn,
        expanded_names=expand_names,
        expanded_shapes=expand_shapes,
        expanded_dtypes=expand_dtypes,
    )

    return compiled, model
