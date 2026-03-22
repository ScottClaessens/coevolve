"""Evaluate PyMC coevolve model log-posterior at Stan primary parameters.

Used by compare_stan_pymc_logprob() for Stan vs PyMC density checks.
Stan returns an unnormalized log density; PyMC returns a fully normalized one.
The normalization offset is cancelled via a shared reference point evaluated
in both backends, so only distribution normalization constants remain.
"""

import os
import sys

import numpy as np
import pytensor.tensor as pt
from pymc.util import get_transformed_name

_pkg_py_dir = os.path.dirname(os.path.abspath(__file__))
if _pkg_py_dir not in sys.path:
    sys.path.insert(0, _pkg_py_dir)
from coev_pymc_model import CoevPymcModel  # noqa: E402


def _reshape_param(constrained, target_shape):
    """Reshape a flat CmdStan array to target_shape using column-major order."""
    a = np.asarray(constrained, dtype=np.float64)
    ts = tuple(int(x) for x in target_shape)
    if a.size != int(np.prod(ts)):
        return a
    if a.ndim <= 1:
        return a.ravel().reshape(ts, order="F")
    return a


def _remap_z_drift_and_orphan_prior(stan_z, tip, pymc_shape):
    """Map Stan z_drift [N_tree, N_seg-1, J] to PyMC [N_tree, N_internal, J].

    Tip-edge z components exist in Stan but not PyMC; returns their unnormalized
    std_normal log kernel so it can be subtracted from Stan's log_prob.
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
                orphan_log_kernel += float(-0.5 * np.sum(zv * zv))

    return pymc_z, orphan_log_kernel


def _invert_lkj_cholesky(L):
    """Invert Stan's LKJ Cholesky factor to PyMC stick-breaking raw parameters."""
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
    """Unnormalized std_normal log kernel: sum(-0.5 * x^2)."""
    a = np.asarray(arr, dtype=np.float64).ravel()
    return float(np.sum(-0.5 * a ** 2))


def _build_pymc_point(model, data_dict, stan_params_dict, ip):
    """Map Stan constrained params to PyMC unconstrained point dict.

    Returns (point_dict, stan_tip_z_log_prior) where stan_tip_z_log_prior is
    the unnormalized kernel for Stan-only tip-edge z_drift elements.
    """
    stan_params = dict(stan_params_dict)
    point = {k: np.asarray(v, dtype=np.float64) for k, v in ip.items()}
    pymc_free_names = {rv.name for rv in model.free_RVs}
    stan_tip_z_log_prior = 0.0

    if "z_drift" in stan_params and "z_drift" in pymc_free_names:
        zrv = next(rv for rv in model.free_RVs if rv.name == "z_drift")
        tr = model.rvs_to_transforms[zrv]
        pkey_z = get_transformed_name("z_drift", tr) if tr is not None else "z_drift"
        tgt = point[pkey_z].shape
        sz, stan_tip_z_log_prior = _remap_z_drift_and_orphan_prior(
            stan_params["z_drift"], data_dict["tip"], tgt,
        )
        stan_params["z_drift"] = sz

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

    # Ordered cutpoints: assign Stan's sorted c{j} directly to PyMC's c{j}_raw.
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
    """Evaluate PyMC log density at Stan constrained parameters.

    Parameters
    ----------
    stan_ref_params : dict or None
        Stan's constrained values at unconstrained=0. When provided, the
        normalization offset is computed at the same point in both backends
        so likelihood terms cancel exactly.

    Returns
    -------
    dict with keys logp_pymc, logp_pymc_0, stan_tip_z_log_prior,
    and stan_terminal_drift_log_prior.
    """
    model = CoevPymcModel().build(data_dict, compile_mode)

    logp_fn = model.compile_logp(jacobian=False)
    ip = model.initial_point(random_seed=None)

    if not isinstance(stan_params, dict):
        stan_params = dict(stan_params)

    pymc_free_names = {rv.name for rv in model.free_RVs}
    prior_only = int(data_dict.get("prior_only", 0))

    # Subtract terminal_drift prior only when Stan includes it (prior_only=FALSE)
    # but PyMC omits it (all-normal fully-observed models).
    stan_terminal_drift_log_prior = 0.0
    if (
        "terminal_drift" in stan_params
        and "terminal_drift" not in pymc_free_names
        and not prior_only
    ):
        stan_terminal_drift_log_prior = _std_normal_log_kernel(
            stan_params["terminal_drift"]
        )

    if stan_ref_params is not None:
        if not isinstance(stan_ref_params, dict):
            stan_ref_params = dict(stan_ref_params)
        ref_point, _ = _build_pymc_point(model, data_dict, stan_ref_params, ip)
        logp_pymc_0 = float(logp_fn(ref_point))
    else:
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
