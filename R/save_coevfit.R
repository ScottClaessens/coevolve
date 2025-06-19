#' Save a fitted \code{coevfit} object to a file
#'
#' A wrapper around \code{\link[base]{saveRDS}} that saves a fitted
#' \code{coevfit} model object while ensuring that all posterior draws and
#' diagnostics are correctly saved by \pkg{cmdstanr}.
#'
#' @param object An object of class \code{coevfit}
#' @param file A string declaring the path where the file should be saved
#' @param ... Other arguments to pass to \code{\link[base]{saveRDS}} besides
#'   \code{object} and \code{file}
#'
#' @returns An .RDS file containing the fitted \code{coevfit} model object
#'
#' @export
save_coevfit <- function(object, file, ...) {
  # stop if object is not of class coevfit
  if (!methods::is(object, "coevfit")) {
    stop2(
      paste0(
        "Argument 'object' must be a fitted coevolutionary model ",
        "of class coevfit."
      )
    )
  }
  # stop if file is not a string of length one
  if (!is.character(file) | length(file) != 1) {
    stop2("Argument 'file' must be a string of length one.")
  }
  # save cmdstanr model object to temporary rds file
  temp_rds_file <- tempfile(fileext = ".rds")
  object$fit$save_object(file = temp_rds_file)
  # replace fit in coevfit object with read-in rds file
  object$fit <- readRDS(temp_rds_file)
  # remove temporary file
  file.remove(temp_rds_file)
  # save coevfit as rds
  saveRDS(object, file = file, ...)
}
