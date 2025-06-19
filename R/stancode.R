#' Expose the Stan code of a fitted model represented by a \code{coevfit} object
#'
#' @param object An object of class \code{coevfit}.
#'
#' @returns Printed Stan code
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
#' # expose stan code
#' stancode(fit)
#' }
#'
#' @export
stancode <- function(object){
  UseMethod("stancode")
}

#' @export
stancode.coevfit <- function(object){
  cat(object$stan_code)
}
