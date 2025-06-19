#' Expose the Stan code of a fitted model represented by a \code{coevfit} object
#'
#' @param object An object of class \code{coevfit}.
#'
#' @returns Printed Stan code
#'
#' @export
stancode <- function(object){
  UseMethod("stancode")
}

#' @export
stancode.coevfit <- function(object){
  cat(object$stan_code)
}
