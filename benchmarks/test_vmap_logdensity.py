"""Test vmap at full log_density + grad level.

Patches the model's _tree_traversal method with a vmapped version
and compares full value_and_grad timing.
"""

import json
import sys
import time

import jax
import jax.numpy as jnp
import numpy as np

jax.config.update("jax_platforms", "cpu")
jax.config.update("jax_enable_x64", True)

sys.path.insert(0, "inst/python")
from coev_jax_model import CoevJaxModel


def patched_tree_traversal(self, params, A_delta_cache, L_VCV_cache,
                           b_delta_cache):
    """Vmapped tree traversal: all trees processed in parallel per level."""
    J = self.J
    N_tree = self.N_tree

    eta = jnp.zeros((N_tree, self.N_seg, J))
    for t in range(N_tree):
        eta = eta.at[t, int(self.root_ids[t])].set(params["eta_anc"][t])

    for _l in range(self.n_levels):
        nids = self.jax_level_node_ids[:, _l]
        pids = self.jax_level_parent_ids[:, _l]
        li = self.jax_level_length_idx[:, _l]
        is_int = self.jax_level_is_internal[:, _l]
        didx = self.jax_level_drift_idx[:, _l]

        def _one_tree_level(eta_t, nids_t, pids_t, li_t, is_int_t, didx_t,
                            z_drift_t):
            _base = (
                jnp.einsum('bij,bj->bi', A_delta_cache[li_t], eta_t[pids_t])
                + b_delta_cache[li_t]
            )
            _noise = jnp.einsum(
                'bij,bj->bi', L_VCV_cache[li_t], z_drift_t[didx_t]
            )
            return eta_t.at[nids_t].set(_base + is_int_t[:, None] * _noise)

        eta = jax.vmap(_one_tree_level)(
            eta, nids, pids, li, is_int, didx, params["z_drift"]
        )

    for t in range(N_tree):
        eta = eta.at[t, int(self.root_ids[t])].set(params["eta_anc"][t])

    # Gather tip L_VCV per tree
    tip_li = jnp.array([
        self.length_index[t][self.tip_to_seg[t]]
        for t in range(N_tree)
    ])
    tip_L_VCV = L_VCV_cache[tip_li]

    eta_trees = [eta[t] for t in range(N_tree)]
    tip_L_VCV_trees = [tip_L_VCV[t] for t in range(N_tree)]
    return eta_trees, tip_L_VCV_trees


def load_model(n_trees):
    with open(f"/tmp/test_vmap_data_{n_trees}.json") as f:
        data = json.load(f)
    for k, v in data.items():
        if isinstance(v, list):
            data[k] = np.array(v)
    return data


def bench_model(data, label, n_iters=300):
    model = CoevJaxModel()
    model.build(data)

    vg = jax.jit(jax.value_and_grad(model.log_density))
    key = jax.random.PRNGKey(42)
    x = jax.random.normal(key, shape=(model.ndim,), dtype=jnp.float64) * 0.5

    lp, grad = vg(x)
    lp.block_until_ready()

    t0 = time.perf_counter()
    for _ in range(n_iters):
        lp, grad = vg(x)
        lp.block_until_ready()
    elapsed = time.perf_counter() - t0
    us = elapsed / n_iters * 1e6
    print(f"  {label}: {us:.0f} us/call")
    return lp


def main():
    for n_trees in [1, 2, 4]:
        print(f"\n--- N_trees = {n_trees} ---")
        data = load_model(n_trees)

        # Original (loop)
        lp_loop = bench_model(data, "loop ")

        # Patched (vmap)
        orig_method = CoevJaxModel._tree_traversal
        CoevJaxModel._tree_traversal = patched_tree_traversal
        lp_vmap = bench_model(data, "vmap ")
        CoevJaxModel._tree_traversal = orig_method

        print(f"  logp diff: {abs(float(lp_loop) - float(lp_vmap)):.2e}")


if __name__ == "__main__":
    main()
