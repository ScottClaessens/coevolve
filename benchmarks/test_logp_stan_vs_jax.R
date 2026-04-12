## Compare Stan and JAX log-joint at shared constrained parameters.
##
## Strategy:
## 1. Compile Stan, run 1 iter for log_prob() + constrain_pars()
## 2. Draw random unconstrained points in Stan's space
## 3. constrain_pars() -> named constrained params
## 4. Stan: log_prob(u, jacobian=FALSE)
## 5. Pass constrained params to JAX, invert transforms, evaluate
## 6. Compare log-joints (JAX log_density minus its jacobian)

devtools::load_all()
reticulate::py_require(c("jax", "numpyro", "nutpie"))

compare_logp <- function(data, variables, id, tree,
                         effects_mat = NULL,
                         lon_lat = NULL, dist_k = NA,
                         prior = NULL,
                         seed = 1L, n_points = 5) {
  distributions <- as.character(variables)

  # --- Stan side ---
  sc <- coev_make_stancode(
    data = data, variables = variables, id = id, tree = tree,
    effects_mat = effects_mat, lon_lat = lon_lat, dist_k = dist_k,
    prior = prior, log_lik = FALSE, prior_only = TRUE
  )
  sd <- coev_make_standata(
    data = data, variables = variables, id = id, tree = tree,
    effects_mat = effects_mat, lon_lat = lon_lat, dist_k = dist_k,
    prior = prior, log_lik = FALSE, prior_only = TRUE
  )
  mod <- cmdstanr::cmdstan_model(
    cmdstanr::write_stan_file(sc), force_recompile = TRUE
  )
  fit <- suppressWarnings(mod$sample(
    data = sd, chains = 1L, iter_warmup = 1L, iter_sampling = 1L,
    seed = 1L, refresh = 0, show_messages = FALSE, init = 0
  ))
  n_upars <- ncol(posterior::as_draws_matrix(
    fit$unconstrain_draws(draws = fit$draws())
  ))

  # --- JAX side ---
  cfg <- coev_make_model_config(
    data = data, variables = variables, id = id, tree = tree,
    effects_mat = effects_mat, lon_lat = lon_lat, dist_k = dist_k,
    prior = prior, prior_only = TRUE
  )
  sd_jax <- embed_model_config(
    standata_to_jax(sd, distributions), cfg
  )
  py_data <- convert_r_to_python_data_jax(sd_jax)

  jax_mod <- load_jax_model_module(convert = FALSE)
  jax <- reticulate::import("jax", convert = FALSE)
  jax$config$update("jax_platforms", "cpu")
  jax$config$update("jax_enable_x64", TRUE)

  model_obj <- jax_mod$CoevJaxModel()
  model_obj$build(py_data)
  # Expose to Python's __main__ so py_run_string can see it
  py <- reticulate::import("builtins", convert = FALSE)
  py_main <- reticulate::import("__main__", convert = FALSE)
  py_main$model <- model_obj

  # Evaluate at multiple random points
  set.seed(seed)
  results <- data.frame(
    point = integer(), stan_logp = numeric(),
    jax_logp = numeric(), diff = numeric()
  )

  for (i in seq_len(n_points)) {
    u_stan <- rnorm(n_upars, 0, 0.3)
    logp_stan <- as.numeric(
      fit$log_prob(u_stan, jacobian = FALSE)
    )
    cpars <- fit$constrain_variables(u_stan)
    cpars_json <- jsonlite::toJSON(
      lapply(cpars, as.numeric),
      auto_unbox = TRUE, digits = 17
    )

    # Evaluate JAX log-joint from same constrained params
    reticulate::py_run_string(sprintf("
import json, jax, jax.numpy as jnp, numpy as np
from coev_jax_model import (
    transform_upper_zero_logdet, transform_lower_zero_logdet,
    transform_ordered_logdet
)

cpars = json.loads('''%s''')
model = model

# Map Stan constrained params to JAX param_info layout
params = {}
for name, shape, transform in model.param_info:
    if name.startswith('_L_R_raw') or name.startswith('_L_residual_raw'):
        params[name] = jnp.zeros(shape)
        continue
    if name in cpars:
        params[name] = jnp.array(
            np.array(cpars[name], dtype=np.float64)
        ).reshape(shape)
    else:
        raise ValueError(
            f'JAX param {name} not in Stan: {list(cpars.keys())}'
        )

# Invert transforms to get unconstrained vector
pieces = []
logdet_jac = 0.0
for name, shape, transform in model.param_info:
    val = params[name]
    if transform == 'none':
        pieces.append(val.ravel())
    elif transform == 'upper_zero':
        raw = jnp.log(jnp.expm1(jnp.clip(-val, 1e-10, None)))
        pieces.append(raw.ravel())
        logdet_jac += float(transform_upper_zero_logdet(raw))
    elif transform == 'lower_zero':
        raw = jnp.log(jnp.expm1(jnp.clip(val, 1e-10, None)))
        pieces.append(raw.ravel())
        logdet_jac += float(transform_lower_zero_logdet(raw))
    elif transform == 'ordered':
        increments = val[1:] - val[:-1]
        raw_rest = jnp.log(jnp.expm1(jnp.clip(increments, 1e-10, None)))
        raw = jnp.concatenate([val[:1], raw_rest])
        pieces.append(raw.ravel())
        logdet_jac += float(transform_ordered_logdet(raw))
    else:
        pieces.append(val.ravel())

u_jax = jnp.concatenate(pieces)
lp_full = float(model.log_density(u_jax).item())
lp_jax_joint = lp_full - logdet_jac
", cpars_json))

    logp_jax <- reticulate::py$lp_jax_joint
    d <- logp_stan - logp_jax
    results <- rbind(results, data.frame(
      point = i, stan_logp = logp_stan,
      jax_logp = logp_jax, diff = d
    ))
    cat(sprintf(
      "  point %d: Stan=%.4f  JAX=%.4f  diff=%.4f\n",
      i, logp_stan, logp_jax, d
    ))
  }

  cat(sprintf(
    "\n  mean |diff| = %.6f, max |diff| = %.6f\n",
    mean(abs(results$diff)), max(abs(results$diff))
  ))

  # Constant offset check: if all diffs are ~same, the models
  # agree up to a normalization constant
  if (nrow(results) > 1) {
    diff_sd <- sd(results$diff)
    cat(sprintf("  sd(diff) = %.6f (should be ~0 if constant offset)\n",
                diff_sd))
  }
  invisible(results)
}

cat("=== Authority (ordered logistic, single tree) ===\n")
compare_logp(
  data = authority$data,
  variables = list(
    political_authority = "ordered_logistic",
    religious_authority = "ordered_logistic"
  ),
  id = "language",
  tree = authority$phylogeny,
  prior = list(A_offdiag = "normal(0, 2)")
)
