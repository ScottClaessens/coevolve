#' Expose the Stan data list for a fitted model represented by a \code{coevfit}
#'   object
#'
#' @param object An object of class \code{coevfit}.
#'
#' @return Named list of data for Stan
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
#' # expose stan data
#' standata(fit)
#' }
#'
#' @export
standata <- function(object){
  UseMethod("standata")
}

#' @export
standata.coevfit <- function(object){
  object$stan_data
}
