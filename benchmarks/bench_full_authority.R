#!/usr/bin/env Rscript
# Full authority model benchmark: 200 warmup / 200 draws, 4 chains, nutpie+JAX.
# Run: Rscript benchmarks/bench_full_authority.R

devtools::load_all(quiet = TRUE)

set.seed(42)

# Capture compile/sample timing from messages
timing <- list(compile = NA_real_, sample = NA_real_)
msg_handler <- function(m) {
  txt <- conditionMessage(m)
  pat <- "Compile:\\s*([0-9.]+)s\\s*\\|\\s*Sample:\\s*([0-9.]+)s"
  match <- regmatches(txt, regexec(pat, txt))[[1]]
  if (length(match) == 3) {
    timing$compile <<- as.numeric(match[2])
    timing$sample <<- as.numeric(match[3])
  }
  invokeRestart("muffleMessage")
}

t_wall <- system.time({
  fit <- withCallingHandlers(
    coev_fit(
      data = authority$data,
      variables = list(
        political_authority = "ordered_logistic",
        religious_authority = "ordered_logistic"
      ),
      id = "language",
      tree = authority$phylogeny,
      prior = list(A_offdiag = "normal(0, 2)"),
      nuts_sampler = "nutpie",
      chains = 4,
      iter_sampling = 200,
      iter_warmup = 200,
      seed = 1
    ),
    message = msg_handler
  )
})[["elapsed"]]

# ESS diagnostics
draws <- fit$fit$draws()
summ <- posterior::summarize_draws(draws, "ess_bulk", "ess_tail")
min_ess_bulk <- min(summ$ess_bulk, na.rm = TRUE)
min_ess_tail <- min(summ$ess_tail, na.rm = TRUE)
ess_per_sec <- if (!is.na(timing$sample) && timing$sample > 0) {
  min_ess_bulk / timing$sample
} else {
  min_ess_bulk / t_wall
}

cat(sprintf(
  paste0(
    "FULL: wall_s=%.1f compile_s=%.1f sample_s=%.1f ",
    "min_ess_bulk=%.0f min_ess_tail=%.0f ess_per_sec=%.1f\n"
  ),
  t_wall, timing$compile, timing$sample,
  min_ess_bulk, min_ess_tail, ess_per_sec
))
