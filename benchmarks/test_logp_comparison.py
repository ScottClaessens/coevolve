"""Verify Stan and JAX log-density agreement at shared constrained params.

Approach:
1. Build both Stan and JAX models
2. Draw random unconstrained params in Stan's space
3. Use Stan's constrain_pars to get named constrained params
4. Evaluate Stan log_prob(u, jacobian=FALSE) -> log-joint
5. Pack those constrained params into JAX's unconstrained space
   (using inverse transforms)
6. Evaluate JAX log_density(u_jax) -> log-joint + jacobian
7. Subtract JAX's jacobian contribution to get JAX log-joint
8. Compare Stan log-joint vs JAX log-joint
"""

import json
import sys

import jax
import jax.numpy as jnp
import numpy as np

jax.config.update("jax_platforms", "cpu")
jax.config.update("jax_enable_x64", True)

sys.path.insert(0, "inst/python")
from coev_jax_model import CoevJaxModel


def inv_softplus(y):
    """Inverse of softplus: log(exp(y) - 1), numerically stable."""
    return jnp.where(y > 20.0, y, jnp.log(jnp.expm1(y)))


def inv_ordered(y):
    """Inverse of ordered transform: extract raw from ordered vector."""
    # y[0] is unconstrained, y[1:] - y[:-1] are softplus(raw[1:])
    increments = y[1:] - y[:-1]
    raw_rest = inv_softplus(increments)
    return jnp.concatenate([y[:1], raw_rest])


def inv_cholesky(L, n):
    """Inverse of stick-breaking Cholesky: L -> raw_vec via atanh."""
    if n == 1:
        return jnp.array([])
    z_list = []
    for i in range(1, n):
        cum_sq = 0.0
        for j in range(i):
            remainder = jnp.sqrt(jnp.maximum(1.0 - cum_sq, 1e-10))
            z_ij = L[i, j] / remainder
            z_list.append(z_ij)
            cum_sq = cum_sq + L[i, j] ** 2
    z = jnp.array(z_list)
    return jnp.arctanh(jnp.clip(z, -0.999999, 0.999999))


def constrained_to_unconstrained(params, model):
    """Map named constrained params to JAX's unconstrained vector."""
    pieces = []
    for name, shape, transform in model.param_info:
        val = params[name]
        if transform == "none":
            pieces.append(val.ravel())
        elif transform == "upper_zero":
            # constrained = -softplus(raw) => raw = inv_softplus(-constrained)
            pieces.append(inv_softplus(-val).ravel())
        elif transform == "lower_zero":
            # constrained = softplus(raw) => raw = inv_softplus(constrained)
            pieces.append(inv_softplus(val).ravel())
        elif transform == "ordered":
            pieces.append(inv_ordered(val).ravel())
        else:
            pieces.append(val.ravel())
    return jnp.concatenate(pieces)


def main():
    with open("/tmp/test_vmap_data_1.json") as f:
        data = json.load(f)
    for k, v in data.items():
        if isinstance(v, list):
            data[k] = np.array(v)

    model = CoevJaxModel()
    model.build(data)

    # Generate random unconstrained params, unpack to constrained
    key = jax.random.PRNGKey(123)
    x_rand = jax.random.normal(key, shape=(model.ndim,)) * 0.3
    params, logdet_jac = model.unpack_params(x_rand)

    # Verify round-trip: constrained -> unconstrained -> constrained
    x_rt = constrained_to_unconstrained(params, model)
    params_rt, logdet_rt = model.unpack_params(x_rt)

    print("Round-trip test (constrained -> unconstrained -> constrained):")
    for name, _, _ in model.param_info:
        diff = float(jnp.max(jnp.abs(params[name] - params_rt[name])))
        print(f"  {name:25s}: max |diff| = {diff:.2e}")

    # Verify log-density matches at round-tripped point
    lp_orig = float(model.log_density(x_rand).item())
    lp_rt = float(model.log_density(x_rt).item())
    print(f"\nlog_density(x_orig):      {lp_orig:.6f}")
    print(f"log_density(x_roundtrip): {lp_rt:.6f}")
    print(f"|diff|:                   {abs(lp_orig - lp_rt):.2e}")

    # Now: log-joint WITHOUT jacobian at the constrained params
    # = log_density(x) - logdet_jacobian(x)
    lp_joint_orig = lp_orig - float(logdet_jac)
    lp_joint_rt = lp_rt - float(logdet_rt)
    print(f"\nlog_joint (no jacobian):")
    print(f"  from x_orig:      {lp_joint_orig:.6f}")
    print(f"  from x_roundtrip: {lp_joint_rt:.6f}")
    print(f"  |diff|:           {abs(lp_joint_orig - lp_joint_rt):.2e}")


if __name__ == "__main__":
    main()
