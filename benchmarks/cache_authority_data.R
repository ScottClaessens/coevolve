#!/usr/bin/env Rscript
# Cache the authority model data dict as JSON for Python-only benchmarking.
# Run: Rscript benchmarks/cache_authority_data.R

devtools::load_all(quiet = TRUE)

variables <- list(
  political_authority = "ordered_logistic",
  religious_authority = "ordered_logistic"
)
distributions <- as.character(variables)

# Generate stan data (effects_mat defaults to all-to-all internally)
sd <- coevolve:::coev_make_standata(
  data = authority$data,
  variables = variables,
  id = "language",
  tree = authority$phylogeny,
  prior = list(A_offdiag = "normal(0, 2)")
)

# Generate model config
model_cfg <- coevolve:::coev_make_model_config(
  data = authority$data,
  variables = variables,
  id = "language",
  tree = authority$phylogeny,
  prior = list(A_offdiag = "normal(0, 2)")
)

# Convert to JAX format and embed config
sd_jax <- coevolve:::embed_model_config(
  coevolve:::standata_to_jax(sd, distributions),
  model_cfg
)

# Attach shape metadata so Python can reconstruct arrays
shapes <- list()
for (nm in names(sd_jax)) {
  val <- sd_jax[[nm]]
  if (is.matrix(val) || is.array(val)) {
    shapes[[paste0(nm, "__shape")]] <- dim(val)
  } else if (is.numeric(val) && length(val) > 1) {
    shapes[[paste0(nm, "__shape")]] <- length(val)
  }
}

out <- c(sd_jax, shapes)
jsonlite::write_json(out, "benchmarks/authority_data_cache.json",
                     auto_unbox = TRUE, digits = 17)
cat("Cached to benchmarks/authority_data_cache.json\n")
cat(sprintf("  Keys: %d\n", length(names(sd_jax))))
