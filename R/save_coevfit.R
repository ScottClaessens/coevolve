#' Save a fitted \code{coevfit} object to a file
#'
#' A wrapper around \code{\link[base]{saveRDS}} that saves a fitted
#' \code{coevfit} model object while ensuring that all posterior draws and
#' diagnostics are correctly saved by \pkg{cmdstanr}.
#'
#' @srrstats {G1.4, G2.1a} Function is documented with expected data types
#'
#' @param object An object of class \code{coevfit}
#' @param file A string declaring the path where the file should be saved
#' @param ... Other arguments to pass to \code{\link[base]{saveRDS}} besides
#'   \code{object} and \code{file}
#'
#' @returns An .RDS file containing the fitted \code{coevfit} model object
#'
#' @examples
#' \dontrun{
#' # fit dynamic coevolutionary model
#' fit <- coev_fit(
#'   data = authority$data,
#'   variables = list(
#'     political_authority = "ordered_logistic",
#'     religious_authority = "ordered_logistic"
#'   ),
#'   id = "language",
#'   tree = authority$phylogeny,
#'   # additional arguments for cmdstanr::sample()
#'   chains = 4,
#'   parallel_chains = 4,
#'   seed = 1
#'   )
#'
#' # save model to .rds file
#' save_coevfit(fit, file = "model.rds")
#' }
#'
#' @export
save_coevfit <- function(object, file, ...) {
  # stop if object is not of class coevfit
  #' @srrstats {G2.1} Assertion on type of input
  if (!methods::is(object, "coevfit")) {
    stop2(
      paste0(
        "Argument 'object' must be a fitted coevolutionary model ",
        "of class coevfit."
      )
    )
  }
  # stop if file is not a string of length one
  #' @srrstats {G2.0, G2.1} Assertion on length and type of input
  if (!is.character(file) || length(file) != 1) {
    stop2("Argument 'file' must be a string of length one.")
  }
  # if file does not end in .rds, add suffix
  #' @srrstats {G4.0} Provide .rds suffix where not provided
  if (!grepl("\\.rds$", file)) {
    file <- paste0(file, ".rds")
  }
  # save fit object - handle both cmdstanr and nutpie
  if (inherits(object$fit, "nutpie_fit")) {
    # For nutpie, we can save the fit directly (it's already an R object)
    # No need for intermediate save/load
  } else {
    # For cmdstanr, save to temporary rds file first
    temp_rds_file <- tempfile(fileext = ".rds")
    object$fit$save_object(file = temp_rds_file)
    # replace fit in coevfit object with read-in rds file
    object$fit <- readRDS(temp_rds_file)
    # remove temporary file
    file.remove(temp_rds_file)
  }
  # save coevfit as rds
  saveRDS(object, file = file, ...)
}
