# PR: Add pure JAX backend via `nuts_sampler = "nutpie"`

> **Status: DRAFT / RFC** -- This is a proposal for discussion, not ready to
> merge. Feedback welcome on the design, scope, and open questions below.

## Motivation

The coevolve Stan model involves heavy linear algebra at every leapfrog
step: matrix exponentials for each branch length, Cholesky decompositions,
multivariate normal evaluations, and a tree traversal that propagates
expectations node-by-node.

A pure JAX implementation of the same log-density can restructure these
computations in ways that Stan's language doesn't express:

- **Batched matrix exponentials and caches**: Stan loops over unique
  branch lengths one at a time (`for (u in 1:N_unique_lengths)`),
  calling `matrix_exp()`, `cholesky_decompose()`, and solving
  `A \ (A_delta - I)` sequentially per length. The JAX version computes
  all of these in single batched calls -- one `matrix_exp_batch()` over
  shape `(n_unique_lengths, J, J)`, one `jnp.linalg.cholesky()` over
  the same batch, etc.
- **Level-batched tree traversal**: Stan propagates states node-by-node
  (`for i in 2:N_seg`). The JAX version groups nodes by tree depth and
  processes all nodes at the same level in one `jnp.einsum()` call.
  Still loops over levels (inherently sequential due to parent-child
  dependencies), but the within-level work is batched.
- **Lightweight dependency footprint**: the Python side requires only
  `jax`, `numpyro` (used solely for `numpyro.distributions`), and
  `nutpie`. An earlier prototype used PyMC/PyTensor, but that stack is
  substantially heavier and harder to install.

In preliminary local benchmarks across 7 model configurations, JAX shows
**2-4x improvement in bottleneck ESS/s** on models where both backends
converge well, with larger gains on models where Stan's NUTS struggles to
mix. One model (repeated measures) shows a slight regression. More
rigorous benchmarking is needed before making strong claims here.

The multiphylo case is worth highlighting since phylogenetic uncertainty
is a distinctive feature of this package and scales the model
substantially. A scaling benchmark varying the number of trees (1-4)
on the authority dataset shows:

| Trees | Wall speedup | Bottleneck ESS/s speedup | Stan max Rhat | JAX max Rhat |
|-------|-------------|-------------------------|---------------|-------------|
| 1 | 7.7x | 2.4x | 1.01 | 1.01 |
| 2 | 7.8x | 25x | 1.53 | 1.12 |
| 3 | 7.1x | 18x | 1.08 | 1.02 |
| 4 | 4.0x | 5.3x | 1.01 | 1.01 |

Wall-clock speedup is consistently 4-8x. When both backends converge
(1 and 4 trees), JAX holds a 2-5x ESS/s advantage. At 2-3 trees Stan's
NUTS struggles to mix (Rhat >> 1), inflating the apparent ESS/s gap.
These are single runs, so exact numbers should be taken directionally.

The goal is to offer a fast, optional alternative without touching the
Stan path that is already tested and trusted.

## How this replaces the existing nutpie-via-BridgeStan path

On `main`, `backend = "nutpie"` compiles the **same Stan model** to C++
via BridgeStan (`nutpie.compile_stan_model`), then samples with nutpie's
Rust NUTS. This gives nutpie's sampler benefits (parallel chains,
low-rank mass matrix adaptation) but still pays the full Stan C++
compilation cost. The gradient computation runs through BridgeStan's C FFI
into the compiled Stan binary.

This branch **replaces that entire path**. Instead of compiling Stan to
C++ and calling it through BridgeStan:

1. `coev_make_model_config()` generates a Python dict describing the model
   structure (variable types, tree topology, priors, GP config, etc.).
2. `inst/python/coev_jax_model.py` (`CoevJaxModel`) builds a pure JAX
   log-density function from that config -- no Stan code is generated or
   compiled.
3. `jax.value_and_grad` provides gradients. nutpie samples using its Rust
   NUTS with those JAX-computed gradients.

The old BridgeStan helpers (`nutpie_compile_stan_model`,
`nutpie_sample`, `convert_nutpie_draws`, etc.) are replaced by the new
`jax_*` modules. The `backend` argument is deprecated in favor of
`nuts_sampler` (`"stan"` or `"nutpie"`).

## Log-density comparison: how we verify correctness

Since the JAX backend reimplements the Stan model's log-density in a
different language, we need a way to confirm they agree numerically.
`compare_stan_jax_logprob()` does this without any MCMC sampling:

1. Compile both the Stan model (via cmdstanr) and the JAX model for the
   same data/config.
2. Run Stan for 1 iteration just to discover the unconstrained parameter
   dimension and get access to `fit$log_prob()`.
3. Evaluate both models at **two** unconstrained points: the origin
   (`u = 0`) and a seeded random draw (`u ~ N(0, 0.5)`).
4. Compare the **difference** in log-density between the two points
   (`logp(u_1) - logp(u_0)`), rather than the raw values. This cancels
   normalization constants that may differ between backends (e.g., Stan
   drops constants that JAX includes, or vice versa).
5. Warn if `|delta_stan - delta_jax| > tol` (default 0.1 on log scale).

Note: Stan and JAX use different unconstrained parameterizations (e.g.,
different LKJ decompositions, Stan includes tip-edge `z_drift`), so the
parameter vectors have different dimensions. Each backend is evaluated at
its own random point -- what must agree is the log-density *surface shape*,
not the raw value at identical coordinates.

## What changed

### New files

| File | Lines | Purpose |
|------|-------|---------|
| `inst/python/coev_jax_model.py` | 1297 | Core JAX log-density, mirrors Stan model blocks (data transforms, parameters, model, generated quantities) |
| `R/jax_helpers.R` | 338 | R-side helpers: Python availability check, JAX model instantiation, expand-fn compilation |
| `R/jax_sample.R` | 159 | Orchestrates JAX sampling: build model, compile, call nutpie, collect draws |
| `R/jax_wrapper.R` | 167 | `jax_fit` S3 class with `draws()`, `summary()`, `metadata()` methods |
| `R/coev_make_model_config.R` | 146 | Generates the Python-side model config dict from R arguments (analogous to `coev_make_stancode` + `coev_make_standata`) |
| `R/compare_stan_jax_logprob.R` | 166 | Diagnostic: evaluates Stan and JAX log-densities at the same unconstrained point to verify numerical agreement |
| `benchmarks/` (various) | -- | Benchmark scripts and a Quarto report comparing Stan vs JAX |

### Modified files

| File | What changed |
|------|-------------|
| `R/coev_fit.R` | New `nuts_sampler` argument (`"stan"` default, `"nutpie"` for JAX). Legacy `backend` argument preserved with deprecation warning. JAX path builds model config instead of Stan code, calls `jax_sample()` instead of cmdstanr. |
| `R/summary.R` | Cutpoint variable parsing now handles both Stan (`c[i,j]`) and JAX (`c1[j]`) naming conventions. |
| `R/stancode.R` | `stancode()` returns the Stan code string when available, or a message indicating JAX backend was used. |
| `tests/testthat/test-coev_fit_nutpie.R` | Rewritten: tests for `normalize_nuts_sampler()`, `check_jax_available()`, `stop_if_jax_not_available()`, plus integration tests for JAX fits (normal, ordered logistic, exact GP, HSGP). |
| `README.Rmd` / `README.md` | Updated sampler docs: `nuts_sampler` API, setup instructions, reticulate config. |
| `DESCRIPTION` | Title/description updated, version bumped to 1.0.0.9001. |
| `NEWS.md` | Added JAX backend entry under New Features. |
| `NAMESPACE` | New exports: `coev_make_model_config`, `compare_stan_jax_logprob`, plus `jax_fit` S3 methods. |

### Feature coverage

The JAX backend supports all model features available in the Stan backend:

- All response distributions (normal, bernoulli_logit, ordered_logistic, poisson_softplus, negative_binomial_softplus, lognormal)
- Repeated measures
- Effects matrices
- Exact GP spatial control (`lon_lat`)
- HSGP approximate GP (`lon_lat` + `dist_k`)
- Correlated drift estimation
- Residual correlation estimation
- Prior-only sampling

## What it does NOT change

- The Stan backend (`nuts_sampler = "stan"`) is completely untouched.
- All existing downstream functions (`coev_plot_*`, `coev_calculate_*`,
  `coev_pred_series`, `extract_samples`, `summary`, `plot`) work with JAX
  fits via the same `coevfit` S3 class.
- No changes to the Stan model code or Stan data generation.

## Open questions

### 1. Maintainability of dual Stan/JAX code paths

This is the biggest concern. The JAX log-density in `coev_jax_model.py`
(1297 lines) is a manual reimplementation of the Stan model. Any future
change to the Stan model (new response distribution, changed prior
structure, modified generated quantities) must be mirrored in the Python
code, or the two backends will silently diverge.

Concretely, the maintenance surface is:

- **Log-density**: Stan's `model {}` block vs. `CoevJaxModel.log_density()`
  and its submethods (`_lp_priors`, `_lp_drift`, `_lp_likelihood`, etc.)
- **Data transforms**: Stan's `transformed data {}` vs. `CoevJaxModel.build()`
- **Parameter transforms**: Stan's unconstrained-to-constrained mapping vs.
  `CoevJaxModel._unpack_params()`
- **Generated quantities**: Stan's `generated quantities {}` vs.
  `CoevJaxModel.generated_quantities()`
- **Cutpoint naming**: `summary.R` now handles both `c[i,j]` (Stan) and
  `c1[j]` (JAX) formats

Mitigations currently in place:
- `compare_stan_jax_logprob()` can verify numerical agreement for any
  model configuration
- The Python code is structured to mirror Stan block names, making it
  easier to cross-reference

Mitigations worth considering:
- CI job that runs `compare_stan_jax_logprob()` across a matrix of model
  configs (all distributions, with/without GP, with/without effects)
- A contributor checklist: "if you change the Stan model, run
  `compare_stan_jax_logprob()` on the affected config"
- Longer term: auto-generate the JAX log-density from a shared model spec

### 2. Python dependency management

Currently Python deps (jax, numpyro, nutpie) are installed via bare `pip`
and reticulate's `"managed"` environment. This is fragile:

- No version pinning -- a JAX or NumPyro update could silently break the
  backend.
- No lockfile or reproducible environment specification.
- Users with existing Python environments may hit conflicts.

**Options to consider:**
- **pixi** (conda-forge based): lockfile, cross-platform, handles CUDA/GPU
  deps natively. Would add a `pixi.toml` to the repo.
- **uv** with a `requirements.txt` or `pyproject.toml`: lightweight, fast,
  pip-compatible. Reticulate already has some uv support.
- **renv + reticulate**: use renv's Python integration to snapshot the venv.
- Pin exact versions in a `requirements.txt` at minimum, regardless of tool.

### 3. Numerical validation strategy

`compare_stan_jax_logprob()` checks that Stan and JAX agree at a single
random unconstrained point. Should this be:

- Run automatically as part of CI (currently it's a manual diagnostic)?
- Extended to check gradients as well as log-density values?
- Part of the test suite with a tolerance threshold?

### 4. API naming

The current API uses `nuts_sampler = "nutpie"` to select the JAX backend,
which is slightly misleading since nutpie is only the sampler -- the key
change is the JAX log-density. Should this be `nuts_sampler = "jax"` or
remain `"nutpie"` for consistency with brms/other packages?

### 5. Generated quantities parity

The JAX backend computes generated quantities (eta, cutpoints) in Python
and maps them to Stan-compatible parameter names. This duplicates logic
from the Stan `generated quantities {}` block. Worth auditing systematically
for parity -- any mismatch here would cause subtle downstream issues in
`coev_plot_*` and `coev_calculate_*` functions that consume these parameters.

### 6. Scope of this PR

This is a large diff (+4,053 / -648 across 38 files). It might be worth
splitting into smaller PRs (e.g., R-side plumbing first, then the Python
model, then tests/benchmarks), though the pieces are tightly coupled.

## CI status

- **lintr**: passing (no lints)
- **R CMD check**: 1 WARNING (CRAN incoming feasibility, expected for
  non-CRAN deps), 0 NOTEs on package code itself
- **Tests**: all passing locally (251s)

## How to test

```r
# Install Python deps
pip install jax numpyro nutpie

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

# Verify numerical agreement with Stan
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
