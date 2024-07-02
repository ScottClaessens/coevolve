#' Extract delta theta from a fitted \code{coevfit} object
#'
#' @param object An object of class \code{coevfit}
#' @param response A character string equal to one of the coevolving variables in the model
#' @param predictor A character string equal to one of the coevolving variables in the model
#'
#' @return Posterior samples in the draws_array format
#' @export
#'
#' @examples
#' \dontrun{
#' # simulate data
#' n <- 20
#' tree <- ape::rcoal(n)
#' d <- data.frame(
#'   id = tree$tip.label,
#'   x = rbinom(n, size = 1, prob = 0.5),
#'   y = ordered(sample(1:4, size = n, replace = TRUE))
#' )
#' # fit dynamic coevolutionary model
#' m <- coev_fit(
#'   data = d,
#'   variables = list(
#'     x = "bernoulli_logit",
#'     y = "ordered_logistic"
#'   ),
#'   id = "id",
#'   tree = tree,
#'   # additional arguments for cmdstanr::sample()
#'   chains = 4,
#'   parallel_chains = 4,
#'   iter_warmup = 500,
#'   seed = 1
#' )
#' # get delta theta
#' coev_get_delta_theta(m, response = "y", predictor = "x")
#' }
coev_get_delta_theta <- function(object, response, predictor) {
  # stop if object is not of class coevfit
  if (!methods::is(object, "coevfit")) {
    stop2("Argument 'object' must be a fitted coevolutionary model of class coevfit.")
  }
  if (!is.character(response) | length(response) != 1) {
    # stop if response not character string of length one
    stop2("Argument 'response' must be a character string of length one.")
  } else if (!(response %in% names(object$variables))) {
    # stop if response not included in model
    stop2("Argument 'response' must be a variable included in the fitted model.")
  }
  if (!is.character(predictor) | length(predictor) != 1) {
    # stop if predictor not character string of length one
    stop2("Argument 'predictor' must be a character string of length one.")
  } else if (!(predictor %in% names(object$variables))) {
    # stop if predictor not included in model
    stop2("Argument 'predictor' must be a variable included in the fitted model.")
  }
  # stop if response and predictor are the same variable
  if (response == predictor) {
    stop2("Argument 'response' and 'predictor' must refer to different variables.")
  }
  # extract posterior draws
  draws <- as_draws_rvars(object$fit$draws())
  # get variables
  A <- draws$A
  b <- draws$b
  # ids for response and predictor
  id_resp <- which(response  == names(object$variables))
  id_pred <- which(predictor == names(object$variables))
  # medians and median absolute deviations for both variables
  eta_resp <- draws$eta[1:nrow(object$data), id_resp]
  eta_pred <- draws$eta[1:nrow(object$data), id_pred]
  med_resp <- rvar_median(eta_resp)
  med_pred <- rvar_median(eta_pred)
  diff_resp <- rvar_mad(eta_resp)
  diff_pred <- rvar_mad(eta_pred)
  # calculate theta for response given value of predictor
  calculate_theta <- function(pred_value) {
    -( pred_value * A[id_resp,id_pred] + b[id_resp] ) / A[id_resp,id_resp]
  }
  # calculate delta theta (change in theta with +1 mad increase in predictor)
  med_theta  <- calculate_theta(med_pred)
  diff_theta <- calculate_theta(med_pred + diff_pred)
  delta_theta <- (diff_theta - med_theta) / diff_resp
  # save as draws array for output
  delta_theta <- as_draws_array(delta_theta)
  dimnames(delta_theta)$variable <- "delta_theta"
  return( delta_theta )
}
