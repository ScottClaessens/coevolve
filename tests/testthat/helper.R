# reload fit from cmdstan csv in fixtures folder
reload_fit <- function(coevfit, filename) {
  coevfit$fit <-
    cmdstanr::as_cmdstan_fit(
      testthat::test_path("fixtures", filename)
    )
  coevfit
}
