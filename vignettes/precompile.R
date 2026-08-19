# Pre-compile vignettes that take a long time to run
# see: https://ropensci.org/blog/2019/12/08/precompute-vignettes/
knitr::knit("vignettes/coevolve.Rmd.orig",   "vignettes/coevolve.Rmd")
knitr::knit("vignettes/priors.Rmd.orig",     "vignettes/priors.Rmd")
knitr::knit("vignettes/missing.Rmd.orig",    "vignettes/missing.Rmd")
knitr::knit("vignettes/repeated.Rmd.orig",   "vignettes/repeated.Rmd")
knitr::knit("vignettes/multiphylo.Rmd.orig", "vignettes/multiphylo.Rmd")
knitr::knit("vignettes/spatial.Rmd.orig",    "vignettes/spatial.Rmd")
knitr::knit("vignettes/compare.Rmd.orig",    "vignettes/compare.Rmd")
knitr::knit("vignettes/jax.Rmd.orig",        "vignettes/jax.Rmd")
knitr::knit("vignettes/ancestral_states.Rmd.orig",
            "vignettes/ancestral_states.Rmd")

# remove "vignettes/" to ensure correct figure paths
edit_figure_paths <- function(file) {
  lines <- readLines(file, warn = FALSE)
  lines <- gsub("vignettes/", "", lines, fixed = TRUE)
  writeLines(lines, file)
}
edit_figure_paths("vignettes/coevolve.Rmd")
edit_figure_paths("vignettes/priors.Rmd")
edit_figure_paths("vignettes/missing.Rmd")
edit_figure_paths("vignettes/repeated.Rmd")
edit_figure_paths("vignettes/multiphylo.Rmd")
edit_figure_paths("vignettes/spatial.Rmd")
edit_figure_paths("vignettes/compare.Rmd")
edit_figure_paths("vignettes/jax.Rmd")
edit_figure_paths("vignettes/ancestral_states.Rmd")
