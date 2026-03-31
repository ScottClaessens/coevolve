"""Lightweight JAX log-density benchmark using cached authority data.

Run: conda run -n nutenv python benchmarks/bench_lightweight.py
"""

import json
import sys
import time

import jax
import jax.numpy as jnp
import numpy as np

jax.config.update("jax_enable_x64", True)
jax.config.update("jax_platforms", "cpu")

sys.path.insert(0, "inst/python")
from coev_jax_model import CoevJaxModel  # noqa: E402

INTEGER_VARS = {
    "node_seq", "parent", "tip", "effects_mat",
    "tip_id", "N_tips", "N_tree", "N_obs", "J", "N_seg",
    "num_effects", "prior_only", "miss", "length_index",
    "tip_to_seg", "N_unique_lengths",
    "n_levels", "max_level_size", "root_ids",
    "level_node_ids", "level_parent_ids",
    "level_length_idx", "level_is_internal",
    "level_drift_idx", "level_sizes",
    "NBgp",
}


def load_cached_data(path="benchmarks/authority_data_cache.json"):
    with open(path) as f:
        raw = json.load(f)

    shapes = {}
    data_keys = []
    for k in list(raw.keys()):
        if k.endswith("__shape"):
            shapes[k.replace("__shape", "")] = raw[k]
        else:
            data_keys.append(k)

    data = {}
    for k in data_keys:
        v = raw[k]
        if isinstance(v, list) and not isinstance(v[0] if v else None, (dict, str)):
            is_int = k in INTEGER_VARS
            arr = np.array(v, dtype=np.int32 if is_int else np.float64)
            if k in shapes:
                s = shapes[k]
                if isinstance(s, list):
                    arr = arr.reshape(s)
            data[k] = arr
        elif isinstance(v, dict):
            # prior_specs or nested dict
            data[k] = v
        elif isinstance(v, (int, float)):
            data[k] = int(v) if k in INTEGER_VARS else float(v)
        else:
            data[k] = v
    return data


def main():
    data = load_cached_data()
    model = CoevJaxModel()
    model.build(data)

    logp_and_grad = jax.jit(jax.value_and_grad(model.log_density))
    x0 = jnp.zeros(model.ndim)

    # JIT compilation
    t_compile_start = time.perf_counter()
    val, grad = logp_and_grad(x0)
    val.block_until_ready()
    t_compile = time.perf_counter() - t_compile_start

    logp_val = float(val)

    # Warm up (2 more calls)
    for _ in range(2):
        val, grad = logp_and_grad(x0)
        val.block_until_ready()

    # Benchmark: 500 reps
    times = []
    for _ in range(500):
        t0 = time.perf_counter()
        val, grad = logp_and_grad(x0)
        val.block_until_ready()
        t1 = time.perf_counter()
        times.append((t1 - t0) * 1e6)

    times = np.array(times)
    med = np.median(times)
    p5 = np.percentile(times, 5)
    p95 = np.percentile(times, 95)

    print(
        f"BENCH: median_us={med:.1f} compile_s={t_compile:.2f} "
        f"logp={logp_val:.4f} p5={p5:.1f} p95={p95:.1f} ndim={model.ndim}"
    )


if __name__ == "__main__":
    main()
