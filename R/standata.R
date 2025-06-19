#' Expose the Stan data list for a fitted model represented by a \code{coevfit}
#'   object
#'
#' @param object An object of class \code{coevfit}.
#'
#' @return Named list of data for Stan
#'
#' @export
standata <- function(object){
  UseMethod("standata")
}

#' @export
standata.coevfit <- function(object){
  object$stan_data
}
