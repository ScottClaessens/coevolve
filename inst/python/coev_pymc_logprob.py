"""Evaluate PyMC coevolve model log-posterior at Stan primary parameters.

Used by compare_stan_pymc_logprob() for Stan vs PyMC density checks.
"""

import numpy as np
import pytensor.tensor as pt
from pymc.util import get_transformed_name


def _reshape_param(constrained, target_shape):
    """Flattened CmdStan arrays -> target shape when sizes match."""
    a = np.asarray(constrained, dtype=np.float64)
    ts = tuple(int(x) for x in target_shape)
    if a.size == int(np.prod(ts)):
        return a.reshape(ts)
    return a


def _remap_z_drift_and_orphan_prior(stan_z, tip, pymc_shape):
    """Map Stan z_drift [N_tree, N_seg-1, J] to PyMC [N_tree, N_internal, J].

    Stan assigns independent std_normal priors to z on every edge step
    i=2..N_seg. PyMC stores only internal-node slots. Tip-edge z components are
    not PyMC parameters; return their log N(0,1) density so it can be
    subtracted from Stan's log_prob for apples-to-apples comparison.
    """
    stan_z = np.asarray(stan_z, dtype=np.float64)
    tip_arr = np.asarray(tip, dtype=np.int32)
    N_tree, _, J = pymc_shape
    N_seg = tip_arr.shape[1]
    exp_size = N_tree * (N_seg - 1) * J
    if stan_z.size == exp_size and stan_z.ndim != 3:
        stan_z = stan_z.reshape((N_tree, N_seg - 1, J))
    n_internal = pymc_shape[1]
    internal_mask = (tip_arr[:, 1:] == 0).astype(np.int32)
    internal_slot = np.cumsum(internal_mask, axis=1) - 1
    internal_slot[internal_mask == 0] = -1

    pymc_z = np.zeros((N_tree, n_internal, J), dtype=np.float64)
    orphan_log_prior = 0.0
    log_norm_const = -0.5 * J * np.log(2.0 * np.pi)

    for t in range(N_tree):
        for j in range(N_seg - 1):
            zv = stan_z[t, j, :]
            if internal_mask[t, j]:
                k = internal_slot[t, j]
                if k < 0 or k >= n_internal:
                    raise ValueError(f"Bad internal_slot {k} for tree {t} step {j}")
                pymc_z[t, k, :] = zv
            else:
                orphan_log_prior += float(-0.5 * np.sum(zv * zv) + log_norm_const)

    return pymc_z, orphan_log_prior


def _std_normal_logp(arr):
    """Sum of std_normal log densities for all elements of arr."""
    a = np.asarray(arr, dtype=np.float64).ravel()
    return float(np.sum(-0.5 * a ** 2 - 0.5 * np.log(2.0 * np.pi)))


def pymc_logprob_at_stan_primary_params(pymc_code, data_dict, stan_params, compile_mode="cpu"):
    """Return PyMC log density and Stan-only z adjustments (see module doc)."""
    ns = {}
    exec(pymc_code, ns)
    build_model = ns["build_model"]
    model = build_model(data_dict, compile_mode)

    logp_fn = model.compile_logp(jacobian=True)
    ip = model.initial_point(random_seed=None)
    point = {k: np.asarray(v, dtype=np.float64) for k, v in ip.items()}

    if not isinstance(stan_params, dict):
        stan_params = dict(stan_params)

    stan_params = dict(stan_params)
    stan_tip_z_log_prior = 0.0
    stan_terminal_drift_log_prior = 0.0

    # terminal_drift: Stan always declares it as a parameter; PyMC omits it when
    # needs_terminal_drift = FALSE (all-normal, fully-observed, no repeated measures).
    # When prior_only=FALSE, Stan evaluates std_normal() per element in the likelihood;
    # subtract that log density so both densities are comparable.
    pymc_free_names = {rv.name for rv in model.free_RVs}
    prior_only = int(data_dict.get("prior_only", 0))
    if "terminal_drift" in stan_params and "terminal_drift" not in pymc_free_names and not prior_only:
        stan_terminal_drift_log_prior = _std_normal_logp(stan_params["terminal_drift"])

    if "z_drift" in stan_params and "z_drift" in pymc_free_names:
        zrv = next(rv for rv in model.free_RVs if rv.name == "z_drift")
        tr = model.rvs_to_transforms[zrv]
        pkey_z = get_transformed_name("z_drift", tr) if tr is not None else "z_drift"
        tgt = point[pkey_z].shape
        sz, stan_tip_z_log_prior = _remap_z_drift_and_orphan_prior(
            stan_params["z_drift"],
            data_dict["tip"],
            tgt,
        )
        stan_params["z_drift"] = sz

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
                val = np.asarray(tr.backward(ci, *rv_inputs).eval())
                uflat[i] = float(val.reshape(-1)[0])
            point[pkey] = uflat.reshape(target_shape)

    logp_pymc = float(logp_fn(point))

    return {
        "logp_pymc": logp_pymc,
        "stan_tip_z_log_prior": float(stan_tip_z_log_prior),
        "stan_terminal_drift_log_prior": float(stan_terminal_drift_log_prior),
    }
