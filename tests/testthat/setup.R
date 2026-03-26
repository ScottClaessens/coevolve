# global teardown: runs after all tests complete
withr::defer({
  try({
    # kill any leftover cmdstan processes
    if (.Platform$OS.type == "unix") {
      system("pkill -f cmdstan", ignore.stdout = TRUE, ignore.stderr = TRUE)
      system("pkill -f stan", ignore.stdout = TRUE, ignore.stderr = TRUE)
    }
    # force garbage collection
    gc()
    # finalise python (if used via reticulate)
    if (requireNamespace("reticulate", quietly = TRUE)) {
      if (reticulate::py_available(initialize = FALSE)) {
        reticulate::py_finalize()
      }
    }
    # extra garbage collection after python shutdown
    gc()
  }, silent = TRUE)
}, envir = testthat::teardown_env())
