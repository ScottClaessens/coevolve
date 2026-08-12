"""Quick experiment: vmap tree traversal over the tree dimension."""

import json
import sys
import time

import jax
import jax.numpy as jnp
import numpy as np

jax.config.update("jax_platforms", "cpu")
jax.config.update("jax_enable_x64", True)

sys.path.insert(0, "inst/python")
from coev_jax_model import CoevJaxModel, ksolve


def tree_traversal_loop(model, params, A_delta_cache, L_VCV_cache,
                        b_delta_cache):
    """Original: Python loop over trees."""
    J = model.J
    eta_trees = []
    for t in range(model.N_tree):
        eta = jnp.zeros((model.N_seg, J))
        eta = eta.at[int(model.root_ids[t])].set(params["eta_anc"][t])
        for _l in range(model.n_levels):
            _nids = model.jax_level_node_ids[t, _l]
            _pids = model.jax_level_parent_ids[t, _l]
            _li = model.jax_level_length_idx[t, _l]
            _is_int = model.jax_level_is_internal[t, _l]
            _didx = model.jax_level_drift_idx[t, _l]
            _base = (
                jnp.einsum('bij,bj->bi', A_delta_cache[_li], eta[_pids])
                + b_delta_cache[_li]
            )
            _noise = jnp.einsum(
                'bij,bj->bi', L_VCV_cache[_li],
                params["z_drift"][t, _didx]
            )
            eta = eta.at[_nids].set(_base + _is_int[:, None] * _noise)
        eta = eta.at[int(model.root_ids[t])].set(params["eta_anc"][t])
        eta_trees.append(eta)
    return eta_trees


def tree_traversal_vmap(model, params, A_delta_cache, L_VCV_cache,
                        b_delta_cache):
    """Vmapped: process all trees in parallel per level."""
    J = model.J
    N_tree = model.N_tree

    eta = jnp.zeros((N_tree, model.N_seg, J))
    for t in range(N_tree):
        eta = eta.at[t, int(model.root_ids[t])].set(params["eta_anc"][t])

    for _l in range(model.n_levels):
        nids = model.jax_level_node_ids[:, _l]
        pids = model.jax_level_parent_ids[:, _l]
        li = model.jax_level_length_idx[:, _l]
        is_int = model.jax_level_is_internal[:, _l]
        didx = model.jax_level_drift_idx[:, _l]

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
        eta = eta.at[t, int(model.root_ids[t])].set(params["eta_anc"][t])

    return [eta[t] for t in range(N_tree)]


def load_model(n_trees):
    with open(f"/tmp/test_vmap_data_{n_trees}.json") as f:
        data = json.load(f)
    for k, v in data.items():
        if isinstance(v, list):
            data[k] = np.array(v)
    model = CoevJaxModel()
    model.build(data)
    return model


def benchmark(model, params, A_delta_cache, L_VCV_cache,
              b_delta_cache, fn, name, n_iters=300):
    # Warmup
    result = fn(model, params, A_delta_cache, L_VCV_cache, b_delta_cache)
    result[0][0].block_until_ready()

    t0 = time.perf_counter()
    for _ in range(n_iters):
        result = fn(model, params, A_delta_cache, L_VCV_cache, b_delta_cache)
        result[0][0].block_until_ready()
    elapsed = time.perf_counter() - t0
    us = elapsed / n_iters * 1e6
    print(f"  {name:5s}: {us:.0f} us/call")
    return result


def main():
    for n_trees in [1, 2, 4]:
        print(f"\n--- N_trees = {n_trees} ---")
        model = load_model(n_trees)
        key = jax.random.PRNGKey(42)
        x = jax.random.normal(key, shape=(model.ndim,), dtype=jnp.float64) * 0.5
        params, _ = model.unpack_params(x)

        A_mat = model._build_A_matrix(params)
        Q, _, _ = model._build_Q_matrix(params)
        Q_inf = ksolve(A_mat, Q, model.J)
        caches = model._compute_caches(A_mat, Q_inf, params)
        A_delta_cache, L_VCV_cache, _, b_delta_cache = caches

        r_loop = benchmark(model, params, A_delta_cache, L_VCV_cache,
                           b_delta_cache, tree_traversal_loop, "loop")
        r_vmap = benchmark(model, params, A_delta_cache, L_VCV_cache,
                           b_delta_cache, tree_traversal_vmap, "vmap")

        for t in range(model.N_tree):
            diff = float(jnp.max(jnp.abs(r_loop[0][t] - r_vmap[0][t])))
            print(f"  tree {t}: max |diff| = {diff:.2e}")


if __name__ == "__main__":
    main()
