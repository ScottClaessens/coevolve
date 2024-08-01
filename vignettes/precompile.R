# Pre-compiled vignettes that take a long time to run
# Must manually move image files from coevolve/ to coevolve/vignettes/ after
# knit, see: https://ropensci.org/blog/2019/12/08/precompute-vignettes/
knitr::knit("vignettes/coevolve.Rmd.orig", "vignettes/coevolve.Rmd")
