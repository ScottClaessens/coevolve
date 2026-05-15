## Benchmark: multiphylo scaling (1-4 trees)
## Compares Stan vs JAX bottleneck ESS/s as number of trees increases.

devtools::load_all()
library(ape)
library(posterior)
library(dplyr)

n_chains  <- 4L
n_warmup  <- 500L
n_samples <- 1000L
seed      <- 42L

fit_and_summarise <- function(model_name, ..., nuts_sampler) {
  t0 <- proc.time()[3]
  fit <- coev_fit(
    ...,
    nuts_sampler    = nuts_sampler,
    chains          = n_chains,
    parallel_chains = n_chains,
    iter_warmup     = n_warmup,
    iter_sampling   = n_samples,
    seed            = seed
  )
  elapsed <- proc.time()[3] - t0
  summ <- fit$fit$summary()
  data.frame(
    model    = model_name,
    backend  = nuts_sampler,
    variable = summ$variable,
    ess_bulk = summ$ess_bulk,
    ess_tail = summ$ess_tail,
    rhat     = summ$rhat,
    elapsed  = elapsed
  )
}

results <- list()

for (n_trees in 1:4) {
  cat("\n========== N_trees =", n_trees, "==========\n")
  tree_list <- rep(list(authority$phylogeny), n_trees)
  class(tree_list) <- "multiPhylo"

  label <- paste0("multiphylo_", n_trees)

  cat("  Fitting Stan...\n")
  r_stan <- fit_and_summarise(
    label,
    data = authority$data,
    variables = list(
      political_authority = "ordered_logistic",
      religious_authority = "ordered_logistic"
    ),
    id = "language",
    tree = tree_list,
    prior = list(A_offdiag = "normal(0, 2)"),
    nuts_sampler = "stan"
  )
  results[[paste0(label, "_stan")]] <- r_stan

  cat("  Fitting JAX...\n")
  r_jax <- fit_and_summarise(
    label,
    data = authority$data,
    variables = list(
      political_authority = "ordered_logistic",
      religious_authority = "ordered_logistic"
    ),
    id = "language",
    tree = tree_list,
    prior = list(A_offdiag = "normal(0, 2)"),
    nuts_sampler = "nutpie"
  )
  results[[paste0(label, "_jax")]] <- r_jax
}

all_results <- bind_rows(results)

# Compute bottleneck ESS/s
summary_table <- all_results %>%
  mutate(ess_per_s = ess_bulk / elapsed) %>%
  group_by(model, backend) %>%
  summarise(
    min_ess_bulk_per_s = min(ess_per_s, na.rm = TRUE),
    median_ess_bulk_per_s = median(ess_per_s, na.rm = TRUE),
    elapsed = first(elapsed),
    max_rhat = max(rhat, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n\n========== RESULTS ==========\n\n")
print(as.data.frame(summary_table), digits = 3)

write.csv(all_results, "benchmarks/multiphylo_scaling_raw.csv",
          row.names = FALSE)
write.csv(as.data.frame(summary_table),
          "benchmarks/multiphylo_scaling_summary.csv",
          row.names = FALSE)
cat("\nSaved to benchmarks/multiphylo_scaling_*.csv\n")
