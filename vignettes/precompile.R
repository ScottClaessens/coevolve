# Pre-compiled vignettes that take a long time to run
# Need to change figure paths in Rmd files afterwards (remove "vignettes/")
# see: https://ropensci.org/blog/2019/12/08/precompute-vignettes/
knitr::knit("vignettes/coevolve.Rmd.orig", "vignettes/coevolve.Rmd")
