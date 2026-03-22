"""Evaluate PyMC coevolve model log-posterior at Stan primary parameters.

Used by compare_stan_pymc_logprob() for Stan vs PyMC density checks.

Normalization convention
------------------------
Stan's log_prob() returns the UNNORMALIZED log density (kernel), dropping
additive constants such as -0.5*log(2π) per Normal parameter.  PyMC's
compile_logp() returns the FULLY NORMALIZED log density.

The caller (compare_stan_pymc_logprob) passes stan_ref_params -- Stan's
constrained parameter values at unconstrained=0 (u_0).  This function
evaluates PyMC at the SAME constrained values, giving logp_pymc_0.  The
offset
  norm_const = logp_pymc_0 - logp_stan_0
is the total normalization constant across all shared parameters, which
cancels when the caller computes abs(lp_stan_adj - (lp_pymc - norm_const)).

By evaluating at the same constrained θ_ref, the likelihood terms at the
reference cancel exactly and only distribution normalization constants remain.
"""

import os
import sys

import numpy as np
import pytensor.tensor as pt
from pymc.util import get_transformed_name

# Import CoevPymcModel from the same package directory
_pkg_py_dir = os.path.dirname(os.path.abspath(__file__))
if _pkg_py_dir not in sys.path:
    sys.path.insert(0, _pkg_py_dir)
from coev_pymc_model import CoevPymcModel  # noqa: E402


def _reshape_param(constrained, target_shape):
    """Flattened CmdStan arrays -> target shape when sizes match.

    Stan stores multi-dimensional parameters in column-major (Fortran) order:
    first index varies fastest.  Use order='F' when the input is a flat 1-D
    vector (or 0-D scalar) from constrain_variables(); leave already-shaped
    multi-dimensional arrays alone.
    """
    a = np.asarray(constrained, dtype=np.float64)
    ts = tuple(int(x) for x in target_shape)
    if a.size != int(np.prod(ts)):
        return a
    if a.ndim <= 1:
        # Scalar or flat vector from constrain_variables: Stan column-major.
        return a.ravel().reshape(ts, order="F")
    # Already multi-dimensional (e.g. from a previous remap step) — return as-is.
    return a


def _remap_z_drift_and_orphan_prior(stan_z, tip, pymc_shape):
    """Map Stan z_drift [N_tree, N_seg-1, J] to PyMC [N_tree, N_internal, J].

    Stan assigns independent std_normal priors to z on every edge step
    i=2..N_seg. PyMC stores only internal-node slots. Tip-edge z components are
    not PyMC parameters; return their UNNORMALIZED log kernel so it can be
    subtracted from Stan's (also unnormalized) log_prob for apples-to-apples
    comparison.
    """
    stan_z = np.asarray(stan_z, dtype=np.float64)
    tip_arr = np.asarray(tip, dtype=np.int32)
    N_tree, _, J = pymc_shape
    N_seg = tip_arr.shape[1]
    exp_size = N_tree * (N_seg - 1) * J
    if stan_z.size == exp_size and stan_z.ndim != 3:
        stan_z = stan_z.ravel().reshape((N_tree, N_seg - 1, J), order="F")
    n_internal = pymc_shape[1]
    internal_mask = (tip_arr[:, 1:] == 0).astype(np.int32)
    internal_slot = np.cumsum(internal_mask, axis=1) - 1
    internal_slot[internal_mask == 0] = -1

    pymc_z = np.zeros((N_tree, n_internal, J), dtype=np.float64)
    orphan_log_kernel = 0.0

    for t in range(N_tree):
        for j in range(N_seg - 1):
            zv = stan_z[t, j, :]
            if internal_mask[t, j]:
                k = internal_slot[t, j]
                if k < 0 or k >= n_internal:
                    raise ValueError(f"Bad internal_slot {k} for tree {t} step {j}")
                pymc_z[t, k, :] = zv
            else:
                # Unnormalized std_normal kernel: -0.5 * sum(z^2)
                orphan_log_kernel += float(-0.5 * np.sum(zv * zv))

    return pymc_z, orphan_log_kernel


def _invert_lkj_cholesky(L):
    """Stan L_R (Cholesky of correlation matrix) -> PyMC _*_raw (pre-tanh stick-breaking)."""
    L = np.asarray(L, dtype=np.float64)
    J = L.shape[0]
    raw = []
    for i in range(1, J):
        cum_sq = 0.0
        for j in range(i):
            z = L[i, j] / np.sqrt(max(1.0 - cum_sq, 1e-10))
            raw.append(np.arctanh(np.clip(z, -1 + 1e-10, 1 - 1e-10)))
            cum_sq += L[i, j] ** 2
    return np.array(raw, dtype=np.float64)


def _std_normal_log_kernel(arr):
    """Sum of UNNORMALIZED std_normal log kernels: sum(-0.5 * x^2).

    Matches Stan's convention of dropping -0.5*log(2π) per element.
    """
    a = np.asarray(arr, dtype=np.float64).ravel()
    return float(np.sum(-0.5 * a ** 2))


def _build_pymc_point(model, data_dict, stan_params_dict, ip):
    """Map Stan constrained params → PyMC unconstrained point dict.

    Returns (point_dict, stan_tip_z_log_prior).  The orphan z prior is the
    unnormalized std_normal kernel for Stan-only tip-edge z_drift elements.
    """
    stan_params = dict(stan_params_dict)
    point = {k: np.asarray(v, dtype=np.float64) for k, v in ip.items()}
    pymc_free_names = {rv.name for rv in model.free_RVs}
    stan_tip_z_log_prior = 0.0

    # Remap z_drift: Stan [N_tree, N_seg-1, J] → PyMC [N_tree, N_internal, J].
    if "z_drift" in stan_params and "z_drift" in pymc_free_names:
        zrv = next(rv for rv in model.free_RVs if rv.name == "z_drift")
        tr = model.rvs_to_transforms[zrv]
        pkey_z = get_transformed_name("z_drift", tr) if tr is not None else "z_drift"
        tgt = point[pkey_z].shape
        sz, stan_tip_z_log_prior = _remap_z_drift_and_orphan_prior(
            stan_params["z_drift"], data_dict["tip"], tgt,
        )
        stan_params["z_drift"] = sz

    # Invert stick-breaking for LKJ Cholesky parameters.
    # PyMC free RV is `_L_R_raw` / `_L_residual_raw`; Stan free RV is `L_R` / `L_residual`.
    for rv in model.free_RVs:
        rname = rv.name
        if rname.startswith("_") and rname.endswith("_raw"):
            stan_key = rname[1:-4]  # e.g. "_L_R_raw" -> "L_R"
            if stan_key in stan_params:
                L = np.asarray(stan_params[stan_key], dtype=np.float64)
                if L.ndim == 2:
                    raw_vec = _invert_lkj_cholesky(L)
                    tr = model.rvs_to_transforms[rv]
                    pkey = get_transformed_name(rname, tr) if tr is not None else rname
                    target_shape = point[pkey].shape
                    point[pkey] = raw_vec.reshape(target_shape)

    # Ordered cutpoints: PyMC free RV is `c{j}_raw` (Normal, unsorted);
    # Stan constrained param is `c{j}` (sorted ordered[]).  Since sort is
    # idempotent on already-sorted values, assign directly to the raw slot.
    for rv in model.free_RVs:
        rname = rv.name
        if rname.endswith("_raw") and not rname.startswith("_"):
            stan_key = rname[:-4]  # e.g. "c2_raw" -> "c2"
            if stan_key in stan_params:
                constrained = np.asarray(stan_params[stan_key], dtype=np.float64)
                tr = model.rvs_to_transforms[rv]
                pkey = get_transformed_name(rname, tr) if tr is not None else rname
                target_shape = point[pkey].shape
                point[pkey] = _reshape_param(constrained, target_shape).astype(np.float64)

    # Apply forward transform (constrained → unconstrained) for remaining params.
    for rv in model.free_RVs:
        rname = rv.name
        if rname not in stan_params:
            continue
        constrained = np.asarray(stan_params[rname], dtype=np.float64)
        tr = model.rvs_to_transforms[rv]
        pkey = get_transformed_name(rname, tr) if tr is not None else rname
        if pkey not in point:
            raise KeyError(
                f"PyMC point key {pkey!r} missing from initial_point; rv={rname!r}"
            )
        target_shape = point[pkey].shape
        if tr is None:
            point[pkey] = _reshape_param(constrained, target_shape).astype(np.float64)
        else:
            cflat = constrained.reshape(-1)
            uflat = np.empty(cflat.shape[0], dtype=np.float64)
            rv_inputs = rv.owner.inputs
            for i in range(cflat.shape[0]):
                ci = pt.as_tensor_variable(np.float64(cflat[i]))
                val = np.asarray(tr.forward(ci, *rv_inputs).eval())
                uflat[i] = float(val.reshape(-1)[0])
            point[pkey] = uflat.reshape(target_shape)

    return point, stan_tip_z_log_prior


def pymc_logprob_at_stan_primary_params(
    data_dict, stan_params, compile_mode="cpu", stan_ref_params=None
):
    """Return PyMC log density, Stan-only z adjustments, and reference logp.

    Parameters
    ----------
    stan_ref_params : dict or None
        Stan's constrained parameter values at unconstrained=0 (u_0).  When
        provided, logp_pymc_0 is evaluated at the SAME constrained values as
        Stan's reference point so that norm_const = logp_pymc_0 - logp_stan_0
        cancels only distribution normalization constants (not likelihood terms).
        When None, falls back to PyMC's initial_point (less accurate for
        prior_only=FALSE models).

    Returns
    -------
    dict with keys:
      logp_pymc           : float -- normalized PyMC logp at target point
      logp_pymc_0         : float -- normalized PyMC logp at reference point
      stan_tip_z_log_prior: float -- unnormalized std_normal kernel for Stan-only
                            orphan tip-edge z_drift elements (to subtract from
                            Stan's log_prob before comparing)
      stan_terminal_drift_log_prior: float -- unnormalized kernel for terminal_drift
                            when Stan has it but PyMC does not
    """
    model = CoevPymcModel().build(data_dict, compile_mode)

    # jacobian=False: eliminates transform-specific Jacobian terms that differ
    # between Stan and PyMC (e.g. Stan uses exp-transform for upper-bounded
    # params; PyMC uses softplus-based interval transform).  Remaining
    # difference is purely distribution normalization constants.
    logp_fn = model.compile_logp(jacobian=False)
    ip = model.initial_point(random_seed=None)

    if not isinstance(stan_params, dict):
        stan_params = dict(stan_params)

    pymc_free_names = {rv.name for rv in model.free_RVs}
    prior_only = int(data_dict.get("prior_only", 0))

    # terminal_drift: Stan evaluates terminal_drift ~ std_normal() inside
    # if (!prior_only).  Only subtract when prior_only=FALSE (the density is
    # actually included in Stan's log_prob) AND PyMC omits terminal_drift
    # (all-normal fully-observed model).
    stan_terminal_drift_log_prior = 0.0
    if (
        "terminal_drift" in stan_params
        and "terminal_drift" not in pymc_free_names
        and not prior_only
    ):
        stan_terminal_drift_log_prior = _std_normal_log_kernel(
            stan_params["terminal_drift"]
        )

    # Compute logp_pymc_0 at Stan's reference constrained values so that
    # norm_const = logp_pymc_0 - logp_stan_0 contains only normalization
    # constants and cancels across both reference and target evaluations.
    if stan_ref_params is not None:
        if not isinstance(stan_ref_params, dict):
            stan_ref_params = dict(stan_ref_params)
        ref_point, _ = _build_pymc_point(model, data_dict, stan_ref_params, ip)
        logp_pymc_0 = float(logp_fn(ref_point))
    else:
        # Fallback: PyMC initial_point (slightly inaccurate for bounded params)
        fallback_point = {k: np.asarray(v, dtype=np.float64) for k, v in ip.items()}
        logp_pymc_0 = float(logp_fn(fallback_point))

    target_point, stan_tip_z_log_prior = _build_pymc_point(
        model, data_dict, stan_params, ip
    )
    logp_pymc = float(logp_fn(target_point))

    return {
        "logp_pymc": logp_pymc,
        "logp_pymc_0": logp_pymc_0,
        "stan_tip_z_log_prior": float(stan_tip_z_log_prior),
        "stan_terminal_drift_log_prior": float(stan_terminal_drift_log_prior),
    }
