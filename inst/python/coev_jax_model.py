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

LOG_2PI = jnp.log(2.0 * jnp.pi)


def _ordered_logistic_logp(predictor, cutpoints, y):
    """Ordered logistic log-PMF. y is 0-indexed category.

    Uses log-space cumulative differences for numerical stability,
    avoiding concat/clip.
    """
    # log-cumulative: log P(Y <= k) = log_sigmoid(c_k - eta)
    # With boundaries: log P(Y <= -1) = -inf, log P(Y <= K-1) = 0
    logcdf = jax.nn.log_sigmoid(cutpoints[None, :] - predictor[:, None])

    K = cutpoints.shape[0] + 1
    n = predictor.shape[0]
    idx = jnp.arange(n)

    # log P(Y = k) = log(P(Y<=k) - P(Y<=k-1))
    # Upper bound: P(Y <= y) — for y == K-1 (last cat), this is 1 → log=0
    log_upper = jnp.where(
        y == K - 1, 0.0, logcdf[idx, y]
    )
    # Lower bound: P(Y <= y-1) — for y == 0 (first cat), this is 0 → log=-inf
    log_lower = jnp.where(
        y == 0, -jnp.inf, logcdf[idx, y - 1]
    )
    # log(exp(a) - exp(b)) = a + log1p(-exp(b-a))  for a > b
    return log_upper + jnp.log1p(-jnp.exp(log_lower - log_upper))


# --------------------------------------------------------------------------
# Math helpers  (cf. Stan 01-functions.stan)
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


def matrix_exp_batch(A_scaled, n_terms=3, n_squarings=4):
    """Matrix exponential via scaling-and-squaring (batched).

    Computes exp(A_scaled * 2^n_squarings) using Horner's method
    for the Taylor polynomial, then squaring n_squarings times.
    """
    shape = A_scaled.shape
    eye = jnp.broadcast_to(jnp.eye(shape[-1])[None, :, :], shape)
    # Horner: I + A*(1/1! + A*(1/2! + A*(1/3!)))
    S = eye * (1.0 / 6.0)
    S = jnp.matmul(A_scaled, S) + eye * 0.5
    S = jnp.matmul(A_scaled, S) + eye
    S = jnp.matmul(A_scaled, S) + eye
    for _ in range(n_squarings):
        S = jnp.matmul(S, S)
    return S


def mvn_chol_logp(value, chol):
    """Multivariate normal log-density given Cholesky factor (batched).

    Computes -0.5*(J*log(2pi) + v'(LL')^{-1}v) - log|L|
    using cho_solve for numerical stability.
    """
    z = jax.lax.linalg.triangular_solve(
        chol, value[..., None], left_side=True, lower=True
    )[..., 0]
    logdet = jnp.sum(jnp.log(jnp.diagonal(chol, axis1=-2, axis2=-1)), axis=-1)
    quad = jnp.sum(z**2, axis=-1)
    J = value.shape[-1]
    return -0.5 * (J * LOG_2PI + quad) - logdet


# --------------------------------------------------------------------------
# GP kernel and spectral density functions
# --------------------------------------------------------------------------


def gp_cov_exp_quad(coords, sigma, rho):
    """Exponentiated-quadratic (RBF) covariance matrix."""
    sq_dist = jnp.sum((coords[:, None, :] - coords[None, :, :]) ** 2, axis=-1)
    return sigma**2 * jnp.exp(-0.5 * sq_dist / rho**2)


def gp_cov_exponential(coords, sigma, rho):
    """Exponential covariance matrix."""
    d = jnp.sqrt(
        jnp.sum((coords[:, None, :] - coords[None, :, :]) ** 2, axis=-1)
    )
    return sigma**2 * jnp.exp(-d / rho)


def gp_cov_matern32(coords, sigma, rho):
    """Matern 3/2 covariance matrix."""
    d = jnp.sqrt(
        jnp.sum((coords[:, None, :] - coords[None, :, :]) ** 2, axis=-1)
    )
    s3 = jnp.sqrt(3.0) * d / rho
    return sigma**2 * (1.0 + s3) * jnp.exp(-s3)


def spd_exp_quad(slambda, sigma, rho):
    """Spectral density for exp_quad kernel (HSGP)."""
    D = slambda.shape[1]
    constant = sigma**2 * jnp.sqrt(2 * jnp.pi) ** D * rho**D
    return constant * jnp.exp(-0.5 * rho**2 * jnp.sum(slambda**2, axis=1))


def spd_exponential(slambda, sigma, rho):
    """Spectral density for exponential kernel (HSGP)."""
    D = slambda.shape[1]
    constant = (
        sigma**2
        * 2**D
        * jnp.pi ** (D / 2.0)
        * jnp.exp(jax.lax.lgamma((D + 1.0) / 2))
        / jnp.sqrt(jnp.pi)
        * rho**D
    )
    expo = -(D + 1.0) / 2
    return constant * (1 + rho**2 * jnp.sum(slambda**2, axis=1)) ** expo


def spd_matern32(slambda, sigma, rho):
    """Spectral density for Matern 3/2 kernel (HSGP)."""
    D = slambda.shape[1]
    constant = (
        sigma**2
        * 2**D
        * jnp.pi ** (D / 2.0)
        * jnp.exp(jax.lax.lgamma((D + 3.0) / 2))
        * 3 ** (3.0 / 2)
        / (0.5 * jnp.sqrt(jnp.pi))
        * rho**D
    )
    expo = -(D + 3.0) / 2
    return constant * (3 + rho**2 * jnp.sum(slambda**2, axis=1)) ** expo


GP_COV_FNS = {
    "exp_quad": gp_cov_exp_quad,
    "exponential": gp_cov_exponential,
    "matern32": gp_cov_matern32,
}

SPD_FNS = {
    "exp_quad": spd_exp_quad,
    "exponential": spd_exponential,
    "matern32": spd_matern32,
}


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
    """Build a pure JAX log-density function from a data_dict.

    Method structure mirrors the Stan template blocks:
      _build_A_matrix    -> 05-transformed-parameters: fill off-diagonal of A
      _build_Q_matrix    -> 05-transformed-parameters: Q = D * (L_R L_R') * D
      _compute_priors    -> 06-model: priors
      _compute_caches    -> 05-transformed-parameters: branch-length caches
      _tree_traversal    -> 05-transformed-parameters: tree loop over segments
      _transformed_params -> 05-transformed-parameters: tdrift, residual_v, dist_v
      _likelihood        -> 06-model: likelihood
    """

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
        self.has_lon_lat = bool(int(data["has_lon_lat"]))
        self.use_hsgp = bool(int(data["use_hsgp"]))
        self.dist_cov_type = str(data["dist_cov_type"])
        self.has_measurement_error = bool(int(data["has_measurement_error"]))

        if self.has_lon_lat:
            if self.use_hsgp:
                self.NBgp = int(data["NBgp"])
                self.Xgp = jnp.array(data["Xgp"], dtype=jnp.float64)
                self.slambda = jnp.array(
                    data["slambda"], dtype=jnp.float64
                )
            else:
                self.coords = jnp.array(
                    data["coords"], dtype=jnp.float64
                )

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

        # Pre-convert level indexing arrays to JAX for JIT closure capture
        self.jax_level_node_ids = jnp.array(self.level_node_ids)
        self.jax_level_parent_ids = jnp.array(self.level_parent_ids)
        self.jax_level_length_idx = jnp.array(self.level_length_idx)
        self.jax_level_drift_idx = jnp.array(self.level_drift_idx)
        self.jax_level_is_internal = jnp.array(
            self.level_is_internal, dtype=jnp.float64
        )

        # Build parameter layout
        self.param_info = self._build_param_layout()
        self.ndim = sum(int(np.prod(shape)) for _, shape, _ in self.param_info)

        # Build expand info for nutpie
        self.expand_info = self._build_expand_info()

        return self

    # ------------------------------------------------------------------
    # Parameter layout  (cf. Stan 04-parameters.stan)
    # ------------------------------------------------------------------

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
        if self.has_lon_lat:
            dz_dim = self.NBgp if self.use_hsgp else self.N_tips
            info.append(("dist_z", (dz_dim, J), "none"))
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
        if self.has_lon_lat:
            vs.append(("sigma_dist", (J,)))
            vs.append(("rho_dist", (J,)))
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

    # ------------------------------------------------------------------
    # Unpack flat vector -> named constrained params
    # ------------------------------------------------------------------

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

    # ------------------------------------------------------------------
    # Transformed parameters  (cf. Stan 05-transformed-parameters.stan)
    # ------------------------------------------------------------------

    def _build_A_matrix(self, params):
        """Build selection matrix A from diagonal and off-diagonal params.

        Mirrors Stan 05-transformed-parameters.stan lines 26-40:
          A = diag_matrix(A_diag);
          { int ticker = 1;
            for (i in 1:J) for (j in 1:J) if (i != j)
              if (effects_mat[i,j] == 1) { A[i,j] = A_offdiag[ticker]; ... }
          }
        """
        J = self.J
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
        return A_mat

    def _build_Q_matrix(self, params):
        """Build drift matrix Q and return (Q, L_R).

        Mirrors Stan 05-transformed-parameters.stan lines 5-10:
          Q = diag_matrix(Q_sigma) * (L_R * L_R') * diag_matrix(Q_sigma);
        or:
          Q = diag_matrix(Q_sigma^2);
        """
        J = self.J
        lp = 0.0
        if self.estimate_correlated_drift:
            n_corr = J * (J - 1) // 2
            if n_corr > 0:
                L_R, lkj_jac = raw_to_cholesky(params["_L_R_raw"], J)
                lp = lp + lkj_jac
                lp = lp + lkj_cholesky_logp(
                    L_R, self.lkj_eta_drift, J
                )
            else:
                L_R = jnp.eye(J)
            D = jnp.diag(params["Q_sigma"])
            Q = D @ (L_R @ L_R.T) @ D
        else:
            L_R = None
            Q = jnp.diag(params["Q_sigma"] ** 2)
        return Q, L_R, lp

    def _compute_caches(self, A_mat, Q_inf, params):
        """Compute branch-length caches for matrix exponentials & drift VCV.

        Mirrors Stan 05-transformed-parameters.stan lines 43-60:
          for (u in 1:N_unique_lengths) {
            A_delta_cache[u] = matrix_exp(A * unique_lengths[u]);
            VCV_cache[u] = Q_inf - quad_form_sym(Q_inf, A_delta_cache[u]');
            L_VCV_cache[u] = cholesky_decompose(VCV_cache[u]);
            A_solve_cache[u] = A \\ add_diag(A_delta_cache[u], -1);
            // symmetrize A_solve_cache[u]
          }

        Returns (A_delta_cache, L_VCV_cache, A_solve_cache, b_delta_cache).
        """
        J = self.J
        eye_J = jnp.eye(J)

        # Batched matrix exponential over unique branch lengths
        A_dt = A_mat[None, :, :] * (
            self.unique_lengths[:, None, None] / 16.0
        )
        A_delta_cache = matrix_exp_batch(A_dt)

        # VCV = Q_inf - A_delta @ Q_inf @ A_delta'  (symmetrized)
        Qi = Q_inf[None, :, :]
        V = Qi - jnp.matmul(
            A_delta_cache, jnp.matmul(Qi, A_delta_cache.transpose(0, 2, 1))
        )
        VCV_cache = 0.5 * (V + V.transpose(0, 2, 1))
        L_VCV_cache = jnp.linalg.cholesky(VCV_cache)

        # A_solve = A^{-1} (A_delta - I)  (symmetrized)
        A_inv = jnp.linalg.solve(A_mat, eye_J)
        As = jnp.matmul(
            A_inv[None, :, :],
            A_delta_cache - jnp.broadcast_to(
                eye_J[None, :, :], A_delta_cache.shape
            ),
        )
        A_solve_cache = 0.5 * (As + As.transpose(0, 2, 1))

        # b_delta = A_solve @ b  (precomputed for tree traversal)
        b_delta_cache = jnp.matmul(
            A_solve_cache, params["b"][None, :, None]
        )[:, :, 0]

        return A_delta_cache, L_VCV_cache, A_solve_cache, b_delta_cache

    def _tree_traversal(self, params, A_delta_cache, L_VCV_cache,
                        b_delta_cache):
        """Propagate ancestral states down each tree.

        Mirrors Stan 05-transformed-parameters.stan lines 61-101:
          for (t in 1:N_tree) {
            eta[t, node_seq[t,1]] = eta_anc[t];
            for (i in 2:N_seg) {
              // eta[t, node_seq[t,i]] = A_delta * eta[parent] + A_solve*b
              //   + L_VCV * z_drift  (internal only)
            }
          }

        Returns (eta_trees, tip_L_VCV_trees) — lists of length N_tree.
        """
        J = self.J
        eta_trees = []
        tip_L_VCV_trees = []
        for t in range(self.N_tree):
            eta = jnp.zeros((self.N_seg, J))
            eta = eta.at[int(self.root_ids[t])].set(params["eta_anc"][t])
            for _l in range(self.n_levels):
                _nids = self.jax_level_node_ids[t, _l]
                _pids = self.jax_level_parent_ids[t, _l]
                _li = self.jax_level_length_idx[t, _l]
                _is_int = self.jax_level_is_internal[t, _l]
                _didx = self.jax_level_drift_idx[t, _l]
                _base = (
                    jnp.einsum('bij,bj->bi', A_delta_cache[_li], eta[_pids])
                    + b_delta_cache[_li]
                )
                _noise = jnp.einsum(
                    'bij,bj->bi', L_VCV_cache[_li],
                    params["z_drift"][t, _didx]
                )
                eta = eta.at[_nids].set(_base + _is_int[:, None] * _noise)
            eta = eta.at[int(self.root_ids[t])].set(params["eta_anc"][t])
            eta_trees.append(eta)
            _tip_li = self.length_index[t][self.tip_to_seg[t]]
            tip_L_VCV_trees.append(L_VCV_cache[_tip_li])

        return eta_trees, tip_L_VCV_trees

    def _transformed_params(self, params, tip_L_VCV_trees, L_residual):
        """Compute tdrift, residual_v, and dist_v.

        Mirrors Stan 05-transformed-parameters.stan lines 103-128:
          tdrift[t,i] = L_VCV_tips[t,i] * to_vector(terminal_drift[t][i,]);
          residual_v = (diag_pre_multiply(sigma_residual, L_residual)
                        * residual_z)';
          dist_v[,j] = cholesky_decompose(dist_cov) * dist_z[,j];
        """
        J = self.J

        # Terminal drift: tdrift[t,i] = L_VCV_tips[t,i] * terminal_drift[t,i]
        tdrift_trees = None
        if self.tdrift and self.needs_terminal_drift:
            tdrift_trees = []
            for t in range(self.N_tree):
                td = jnp.matmul(
                    tip_L_VCV_trees[t],
                    params["terminal_drift"][t][:, :, None],
                )[:, :, 0]
                tdrift_trees.append(td)

        # Residual effects: residual_v = (D * L_residual * residual_z)'
        residual_v = None
        if self.residual_v_flag and self.repeated:
            D_res = jnp.diag(params["sigma_residual"])
            residual_v = (D_res @ L_residual @ params["residual_z"]).T

        # Spatial GP effects (N_tips x J)
        dist_v = None
        if self.has_lon_lat:
            dist_v = jnp.zeros((self.N_tips, J))
            if self.use_hsgp:
                spd_fn = SPD_FNS[self.dist_cov_type]
                for j in range(J):
                    spd_vals = spd_fn(
                        self.slambda,
                        params["sigma_dist"][j],
                        params["rho_dist"][j],
                    )
                    rgp = jnp.sqrt(spd_vals) * params["dist_z"][:, j]
                    dist_v = dist_v.at[:, j].set(self.Xgp @ rgp)
            else:
                cov_fn = GP_COV_FNS[self.dist_cov_type]
                for j in range(J):
                    K = cov_fn(
                        self.coords,
                        params["sigma_dist"][j],
                        params["rho_dist"][j],
                    )
                    K = K + 1e-12 * jnp.eye(self.N_tips)
                    L_K = jnp.linalg.cholesky(K)
                    dist_v = dist_v.at[:, j].set(
                        L_K @ params["dist_z"][:, j]
                    )

        return tdrift_trees, residual_v, dist_v

    # ------------------------------------------------------------------
    # Priors  (cf. Stan 06-model.stan, lines 1-51)
    # ------------------------------------------------------------------

    def _compute_priors(self, params, L_residual):
        """Evaluate all prior log-densities.

        Mirrors the prior section of Stan 06-model.stan.
        """
        J = self.J
        lp = 0.0

        # b ~ prior_b
        lp = lp + prior_logp(params["b"], self.prior_specs["b"])

        # for (t in 1:N_tree) {
        #   eta_anc[t] ~ prior_eta_anc;
        #   z_drift[t, i] ~ std_normal();
        #   terminal_drift[t] ~ std_normal();  (conditional)
        # }
        lp = lp + prior_logp(params["eta_anc"], self.prior_specs["eta_anc"])
        lp = lp + dist.Normal(0.0, 1.0).log_prob(params["z_drift"]).sum()

        if self.needs_terminal_drift:
            has_normal = len(self.normal_j0) > 0
            if has_normal and not self.repeated:
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

        # A_offdiag ~ prior_A_offdiag;  A_diag ~ prior_A_diag
        lp = lp + prior_logp(params["A_diag"], self.prior_specs["A_diag"])
        if self.n_offdiag > 0:
            lp = lp + prior_logp(
                params["A_offdiag"], self.prior_specs["A_offdiag"]
            )

        # Q_sigma ~ prior_Q_sigma
        lp = lp + prior_logp(params["Q_sigma"], self.prior_specs["Q_sigma"])

        # c ~ prior_c  (ordered cutpoints — batched)
        if self.ordered_j:
            all_c = jnp.concatenate(
                [params[f"c{j1}"] for j1 in self.ordered_j]
            )
            lp = lp + prior_logp(all_c, self.prior_specs["c"])

        # phi ~ normal(inv_overdisp, inv_overdisp)  or manual prior
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

        # shape ~ prior_shape  (gamma shape)
        for j0 in self.gamma_j0:
            lp = lp + prior_logp(
                params[f"shape{j0 + 1}"], self.prior_specs["shape"]
            )

        # GP priors: dist_z ~ std_normal(); sigma_dist, rho_dist ~ priors
        if self.has_lon_lat:
            lp = lp + dist.Normal(0.0, 1.0).log_prob(params["dist_z"]).sum()
            lp = lp + prior_logp(
                params["sigma_dist"], self.prior_specs["sigma_dist"]
            )
            lp = lp + prior_logp(
                params["rho_dist"], self.prior_specs["rho_dist"]
            )

        # Residual priors (repeated measures)
        # Matches Stan 06-model.stan lines 35-51:
        #   for normal vars: if (miss[i,j]==0) residual_z[j,i] ~ std_normal()
        #   for non-normal vars: residual_z[j,i] ~ std_normal()
        if self.repeated:
            rz = params["residual_z"]  # (J, N_obs)
            rz_lp = jnp.zeros(())
            for j0 in self.normal_j0:
                obs_mask = self.miss[:, j0] == 0
                rz_lp = rz_lp + jnp.where(
                    obs_mask,
                    dist.Normal(0.0, 1.0).log_prob(rz[j0]),
                    0.0,
                ).sum()
            for j0 in self.nonnormal_j0:
                rz_lp = rz_lp + (
                    dist.Normal(0.0, 1.0).log_prob(rz[j0]).sum()
                )
            lp = lp + rz_lp

            # sigma_residual ~ prior_sigma_residual
            lp = lp + prior_logp(
                params["sigma_residual"], self.prior_specs["sigma_residual"]
            )

            # L_residual ~ lkj_corr_cholesky(eta)
            n_corr_res = J * (J - 1) // 2
            if n_corr_res > 0:
                lp = lp + lkj_cholesky_logp(
                    L_residual, self.lkj_eta_residual, J
                )

        return lp

    # ------------------------------------------------------------------
    # Likelihood  (cf. Stan 06-model.stan, lines 53-104)
    # ------------------------------------------------------------------

    def _likelihood(self, params, eta_trees, tip_L_VCV_trees,
                    tdrift_trees, residual_v, L_residual, dist_v=None):
        """Compute the log-likelihood contribution.

        Mirrors Stan 06-model.stan:
          for (i in 1:N_obs) {
            vector[N_tree] lp = rep_vector(0.0, N_tree);
            for (t in 1:N_tree) { ... }
            target += log_sum_exp(lp);
          }
        """
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

            # Linear model base: eta[t, tip_id[i]][j] + dist_v[tip_id[i], j]
            def base_lmod(j0, _eta_obs=eta_obs, _dist_v=dist_v):
                lm = _eta_obs[:, j0]
                if _dist_v is not None:
                    lm = lm + _dist_v[tid, j0]
                return lm

            # --- Normal variables ---
            # Non-repeated: set_tdrift path (tdrift vector evaluated under
            #   MVN with L_VCV_tips covariance)
            # Repeated: set_residuals path (residuals evaluated under
            #   MVN with diag(sigma_residual) @ L_residual covariance)

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
                    L_cov_res = jnp.broadcast_to(
                        L_cov_res, (self.N_obs, J, J)
                    )
                obs_lp = obs_lp + mvn_chol_logp(residuals, L_cov_res)

            # --- Non-normal likelihoods ---
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
                    y_int = (y_safe - 1).astype(jnp.int32)
                    ll = _ordered_logistic_logp(lmod, cutpoints, y_int)
                elif d == "poisson_softplus":
                    y_safe = jnp.where(miss_j, 0.0, y_obs)
                    mu = self.obs_means[j0] * jax.nn.softplus(lmod)
                    ll = dist.Poisson(rate=mu).log_prob(
                        y_safe.astype(jnp.int32)
                    )
                elif d == "negative_binomial_softplus":
                    y_safe = jnp.where(miss_j, 0.0, y_obs)
                    mu = self.obs_means[j0] * jax.nn.softplus(lmod)
                    phi = params[f"phi{j1}"].squeeze()
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

    # ------------------------------------------------------------------
    # Main log-density  (orchestrates all blocks)
    # ------------------------------------------------------------------

    def log_density(self, x):
        """Compute log-density at unconstrained parameter vector x.

        This is the function that gets JIT-compiled and differentiated by JAX.
        Mirrors the Stan block flow:
          transformed parameters -> model { priors; likelihood }
        """
        params, logdet_jac = self.unpack_params(x)
        J = self.J
        lp = logdet_jac

        # --- Build A and Q matrices (05-transformed-parameters) ---
        A_mat = self._build_A_matrix(params)
        Q, L_R, lp_Q = self._build_Q_matrix(params)
        lp = lp + lp_Q

        # --- Build L_residual (needed by both priors and likelihood) ---
        L_residual = None
        if self.repeated:
            n_corr_res = J * (J - 1) // 2
            if n_corr_res > 0:
                L_residual, lkj_jac_res = raw_to_cholesky(
                    params["_L_residual_raw"], J
                )
                lp = lp + lkj_jac_res
            else:
                L_residual = jnp.eye(J)

        # --- Priors (06-model) ---
        lp = lp + self._compute_priors(params, L_residual)

        # --- Branch-length caches (05-transformed-parameters) ---
        Q_inf = ksolve(A_mat, Q, J)
        A_delta_cache, L_VCV_cache, A_solve_cache, b_delta_cache = (
            self._compute_caches(A_mat, Q_inf, params)
        )

        # --- Tree traversal (05-transformed-parameters) ---
        eta_trees, tip_L_VCV_trees = self._tree_traversal(
            params, A_delta_cache, L_VCV_cache, b_delta_cache
        )

        # --- Derived quantities: tdrift, residual_v, dist_v ---
        tdrift_trees, residual_v, dist_v = self._transformed_params(
            params, tip_L_VCV_trees, L_residual
        )

        # --- Likelihood (06-model) ---
        if not self.prior_only:
            lp = lp + self._likelihood(
                params,
                eta_trees,
                tip_L_VCV_trees,
                tdrift_trees,
                residual_v,
                L_residual if self.repeated else None,
                dist_v,
            )

        return lp

    # ------------------------------------------------------------------
    # Expand function for nutpie output
    # ------------------------------------------------------------------

    def make_expand_fn(self, seed1=0, seed2=0, chain=0):
        """Return a function that maps flat unconstrained -> named params.

        Keys are returned in the exact order of _expand_var_list().
        The JAX computation is JIT-compiled; only the final
        jax->numpy conversion runs in Python.
        """
        J = self.J
        var_order = [n for n, _ in self._expand_var_list()]

        # JIT-compile the pure JAX part
        @jax.jit
        def _expand_jax(x):
            params, _ = self.unpack_params(x)
            vals = []

            # A and Q — use shared builders
            A_mat = self._build_A_matrix(params)
            Q, L_R, _ = self._build_Q_matrix(params)
            vals.append(A_mat)  # A
            vals.append(Q)      # Q

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

            if self.has_lon_lat:
                vals.append(params["sigma_dist"])
                vals.append(params["rho_dist"])

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


# --------------------------------------------------------------------------
# Entry points (called from R via reticulate)
# --------------------------------------------------------------------------


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


_sampling_state = {}


def start_sampling(compiled_model, **kwargs):
    """Start nutpie sampling in a background thread.

    Returns immediately so R can poll via check_sampling().
    """
    import nutpie
    import threading
    import time

    state = {
        "result": None,
        "error": None,
        "done": False,
        "start_time": time.time(),
    }

    def run():
        try:
            state["result"] = nutpie.sample(
                compiled_model, progress_bar=False, **kwargs
            )
        except Exception as e:
            state["error"] = str(e)
        finally:
            state["done"] = True

    t = threading.Thread(target=run, daemon=True)
    t.start()
    state["_thread"] = t
    _sampling_state["current"] = state


def check_sampling():
    """Return sampling status (called from R polling loop)."""
    import time

    state = _sampling_state.get("current", {})
    return {
        "done": state.get("done", True),
        "elapsed": time.time() - state.get("start_time", time.time()),
        "error": state.get("error"),
    }


def collect_result():
    """Block until sampling finishes and return the ArviZ trace."""
    state = _sampling_state.get("current")
    if not state:
        raise RuntimeError("No active sampling session")
    state["_thread"].join(timeout=600)
    if state["error"]:
        raise RuntimeError(state["error"])
    return state["result"]
