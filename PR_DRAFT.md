# PR: Add pure JAX backend via `nuts_sampler = "nutpie"`

> Feedback welcome on design, scope, and open questions below.

## Motivation

JAX traces the entire log-density function into a single XLA computation
graph, which is then compiled holistically — enabling operation fusion,
memory optimization, and efficient vectorization that Stan's
per-operation C++ compilation doesn't support. This also makes it
natural to express batched operations that Stan's language requires loops
for:

- **Batched matrix operations.** Stan loops over unique branch lengths
  one at a time. JAX computes matrix exponentials, Cholesky
  decompositions, and drift caches for all branch lengths in single
  batched calls.
- **Level-batched tree traversal.** Stan propagates states node-by-node.
  JAX groups nodes by tree depth and processes all nodes at the same
  level in one `einsum` call. Levels are still sequential (parent-child
  dependency), but within-level work is batched.
- **Vmap across trees.** For multiphylo models, `jax.vmap` runs the
  tree traversal, derived quantities, and likelihood in parallel across
  trees — something Stan's language has no mechanism for.

On the authority dataset (97 taxa, 2 ordered-logistic variables), JAX
is **3-5x faster wall-clock** than Stan across 1-50 trees.
Population-parameter ESS/s (effective samples per second) is **modestly
better** with JAX at moderate tree counts (~1.5x median at 10 trees)
and roughly comparable at 50 trees. JAX produces fewer effective samples
per draw than Stan, but the per-draw speed advantage more than
compensates. These are local CPU benchmarks; GPU acceleration is
supported but not yet tested.

## How this replaces the existing nutpie-via-BridgeStan path

On `main`, `backend = "nutpie"` compiles the same Stan model to C++
via BridgeStan (`nutpie.compile_stan_model`), then samples with nutpie's
Rust NUTS. This uses nutpie's Rust NUTS sampler but the model itself is still the
compiled Stan C++ binary (via BridgeStan). The log-density and gradient
computations are the same as Stan — nutpie just replaces the sampler,
not the model evaluation.

This branch replaces that entire path. Instead of compiling Stan to
C++ and calling it through BridgeStan:

1. `coev_make_model_config()` generates a Python dict describing the
   model structure (variable types, tree structure, priors, GP config).
2. `inst/python/coev_jax_model.py` (`CoevJaxModel`) builds a pure JAX
   log-density function from that config — no Stan code is generated or
   compiled.
3. `jax.value_and_grad` provides gradients; nutpie samples using its
   Rust NUTS with those JAX-computed gradients.

The old BridgeStan helpers (`nutpie_compile_stan_model`,
`nutpie_sample`, `convert_nutpie_draws`) are replaced by the new `jax_*`
modules. The `backend` argument is deprecated in favor of
`nuts_sampler` (`"stan"` or `"nutpie"`).

## Correctness verification

Since the JAX backend reimplements the Stan log-density in a different
language, correctness is verified by evaluating both models at the
**same unconstrained parameter vector** and comparing log-density and
gradient. This is analogous to the approach used by https://github.com/pymc-labs/alchemize.

`compare_stan_jax_logprob()` does this:

1. Compile both models for the same data/config.
2. Draw `n_points` random unconstrained vectors with a fixed seed.
3. At each point, compute `log_prob()` and `grad_log_prob()` in both
   backends.
4. Check that the Stan/JAX log-density **difference is constant** across
   points (any constant offset comes from normalization terms Stan
   drops) and that **gradients agree** to a tolerance.

This comparison is wired into the test suite as
`tests/testthat/test-logp_stan_jax.R`, which runs 12 prior-only
configurations (all response distributions, GP/HSGP, repeated measures,
multiphylo, measurement error, effects matrices, correlated drift on/off)
plus 7 full-likelihood configurations.

**Prior-only tests pass at machine precision** (max grad diff ~1e-16)
**Likelihood tests pass at ~1e-3 tolerance**

## What changed

### New files

| File | Lines | Purpose |
|------|-------|---------|
| `inst/python/coev_jax_model.py` | 1297 | Core JAX log-density, mirrors Stan blocks (data transforms, parameters, model, generated quantities) |
| `R/jax_helpers.R` | 338 | R-side helpers: Python availability check, JAX model instantiation, expand-fn compilation |
| `R/jax_sample.R` | 159 | Orchestrates JAX sampling: build model, JIT-compile, call nutpie, collect draws |
| `R/jax_wrapper.R` | 167 | `jax_fit` S3 class with `draws()`, `summary()`, `metadata()` methods |
| `R/coev_make_model_config.R` | 146 | Generates Python-side model config dict (analogous to `coev_make_stancode` + `coev_make_standata`) |
| `R/compare_stan_jax_logprob.R` | ~170 | Evaluates Stan and JAX log-densities + gradients at the same unconstrained point |
| `tests/testthat/test-logp_stan_jax.R` | ~340 | 19 tests covering all supported configs (prior-only and full-likelihood) |
| `benchmarks/` (various) | — | Benchmark scripts, debug utilities, PR drafting |

### Modified files

| File | What changed |
|------|-------------|
| `R/coev_fit.R` | New `nuts_sampler` argument (`"stan"` default, `"nutpie"` for JAX). Legacy `backend` argument preserved with deprecation warning. |
| `R/summary.R` | Cutpoint variable parsing now handles both Stan (`c[i,j]`) and JAX (`c1[j]`) naming conventions. |
| `R/stancode.R` | Returns Stan code when available, or a message indicating JAX backend was used. |
| `tests/testthat/test-coev_fit_nutpie.R` | Rewritten for the new `nuts_sampler` API plus integration tests for JAX fits. |
| `README.Rmd` / `README.md` | Updated sampler docs: `nuts_sampler` API, setup instructions, reticulate config. |
| `DESCRIPTION` | Title/description updated, version bumped to 1.0.0.9001. |
| `NEWS.md` | Added JAX backend entry. |
| `NAMESPACE` | New exports: `coev_make_model_config`, `compare_stan_jax_logprob`, plus `jax_fit` S3 methods. |

### Feature coverage

The JAX backend supports all model features available in the Stan backend:

- All response distributions (`normal`, `bernoulli_logit`,
  `ordered_logistic`, `poisson_softplus`, `negative_binomial_softplus`,
  `gamma_log`)
- Repeated measures
- Effects matrices
- Exact GP spatial control (`lon_lat`)
- HSGP approximate GP (`lon_lat` + `dist_k`)
- Correlated drift estimation (on/off)
- Residual correlation estimation
- Measurement error
- Multiphylo (phylogenetic uncertainty)
- Prior-only sampling

All combinations are covered by the logp agreement test suite.

## What it does NOT change

- The Stan backend (`nuts_sampler = "stan"`) is untouched.
- All downstream functions (`coev_plot_*`, `coev_calculate_*`,
  `coev_pred_series`, `extract_samples`, `summary`, `plot`) work with
  JAX fits via the same `coevfit` S3 class.
- No changes to the Stan model code or Stan data generation.

## Open questions

### 1. Maintainability of dual Stan/JAX code paths

This is the main long-term concern. The JAX log-density
(`coev_jax_model.py`, 1297 lines) is a manual reimplementation of the
Stan model. Any future change to the Stan model must be mirrored in the
Python code, or the two backends will diverge.

Mitigations in place:

- `tests/testthat/test-logp_stan_jax.R` runs the logp-agreement check
  across all supported configs — will catch most divergence.
- The Python code is structured to mirror Stan block names
  (`_compute_priors` ~ `model{}` priors, `_likelihood` ~ likelihood
  block, etc.) for easier cross-reference.

Mitigations worth considering:

- Wire `test-logp_stan_jax.R` into CI (currently runs locally; would add
  significant CI runtime since it compiles Stan).
- Contributor checklist: "if you change the Stan model, add a matching
  case to `test-logp_stan_jax.R`."
- Longer term: auto-generate the JAX log-density from a shared model
  spec.

### 2. Python dependency management

Python deps are now pinned in `inst/python/requirements.txt`
(jax==0.5.3, numpyro==0.19.0, nutpie==0.16.8). `check_jax_available()`
passes these pinned specs to `reticulate::py_require()`, which resolves
them via uv into a managed environment.

This is a substantial improvement over unpinned deps (a JAX update
broke numpyro during development), but still has limitations:

- GPU support would need `jax[cuda]` which has different deps.
- Users with existing Python environments may still hit conflicts.
- A proper lockfile (e.g., pixi.toml or uv.lock) would be more
  robust than pinned versions in requirements.txt.

### 3. API naming

`nuts_sampler = "nutpie"` selects the JAX backend, which is slightly
misleading — nutpie is the sampler, but the key change is the JAX
log-density. Should this be `nuts_sampler = "jax"` or remain
`"nutpie"` for brms/etc. consistency?

### 4. Scope of this PR

Large diff (~4,000 lines added). Could be split into: R-side plumbing
first, then the Python model, then tests/benchmarks. The pieces are
tightly coupled, so splitting may not buy much.

### 5. CI runtime cost

The logp comparison test suite takes ~10-15 min locally because each
test compiles Stan. If wired into CI it becomes the dominant cost. Could
be gated on `COEVOLVE_EXTENDED_TESTS=true` or run nightly instead of on
every PR.

## CI status

- **lintr**: passing (no lints)
- **R CMD check**: 1 WARNING (CRAN incoming feasibility, expected for
  non-CRAN deps), 0 NOTEs on package code itself
- **Test suite**: 1091 tests pass locally (including 19 logp agreement
  tests at machine precision for priors, ~1e-3 for likelihood)

## How to test

```r
# Install Python deps
# pip install jax numpyro nutpie

# Fit with JAX backend
fit <- coev_fit(
  data = authority$data,
  variables = list(
    political_authority = "ordered_logistic",
    religious_authority = "ordered_logistic"
  ),
  id = "language",
  tree = authority$phylogeny,
  nuts_sampler = "nutpie"
)

# Verify log-density agreement with Stan
compare_stan_jax_logprob(
  data = authority$data,
  variables = list(
    political_authority = "ordered_logistic",
    religious_authority = "ordered_logistic"
  ),
  id = "language",
  tree = authority$phylogeny
)
```
