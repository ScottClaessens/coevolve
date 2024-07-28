#' Calculate delta theta from a fitted \code{coevfit} object
#'
#' @param object An object of class \code{coevfit}
#' @param response A character string equal to one of the coevolving variables
#'   in the model
#' @param predictor A character string equal to one of the coevolving variables
#'   in the model
#'
#' @return Posterior samples in the draws_array format
#' @export
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
#' # calculate delta theta
#' coev_calculate_delta_theta(
#'   object = fit,
#'   response = "political_authority",
#'   predictor = "religious_authority"
#'   )
#' }
coev_calculate_delta_theta <- function(object, response, predictor) {
  # stop if object is not of class coevfit
  if (!methods::is(object, "coevfit")) {
    stop2(
      paste0(
        "Argument 'object' must be a fitted coevolutionary model ",
        "of class coevfit."
        )
      )
  }
  if (!is.character(response) | length(response) != 1) {
    # stop if response not character string of length one
    stop2("Argument 'response' must be a character string of length one.")
  } else if (!(response %in% names(object$variables))) {
    # stop if response not included in model
    stop2(
      "Argument 'response' must be a variable included in the fitted model."
      )
  }
  if (!is.character(predictor) | length(predictor) != 1) {
    # stop if predictor not character string of length one
    stop2("Argument 'predictor' must be a character string of length one.")
  } else if (!(predictor %in% names(object$variables))) {
    # stop if predictor not included in model
    stop2(
      "Argument 'predictor' must be a variable included in the fitted model."
      )
  }
  # stop if response and predictor are the same variable
  if (response == predictor) {
    stop2(
      "Argument 'response' and 'predictor' must refer to different variables."
      )
  }
  # extract posterior draws
  draws <- posterior::as_draws_rvars(object$fit$draws())
  # ids for response and predictor
  id_resp <- which(response  == names(object$variables))
  id_pred <- which(predictor == names(object$variables))
  # medians and median absolute deviations for all variables
  eta  <- draws$eta[1:nrow(object$data),]
  med  <- apply(eta, 2, posterior::rvar_median)
  diff <- apply(eta, 2, posterior::rvar_mad)
  # construct intervention values list
  values1 <- list()
  values2 <- list()
  for (j in 1:length(names(object$variables))) {
    variables <- names(object$variables)
    if (j == id_resp) {
      # if variable is response, set to NA
      values1[[variables[j]]] <- NA
      values2[[variables[j]]] <- NA
    } else if (j == id_pred) {
      # else if variable is predictor, set to median and median+mad
      values1[[variables[j]]] <- stats::median(med[[j]])
      values2[[variables[j]]] <- stats::median(med[[j]] + diff[[j]])
    } else {
      # else for all other variables, set to median
      values1[[variables[j]]] <- stats::median(med[[j]])
      values2[[variables[j]]] <- stats::median(med[[j]])
    }
  }
  # calculate delta theta (change in theta with +1 mad increase in predictor)
  # all other variables are held at their median values
  theta1 <- coev_calculate_theta(object, intervention_values = values1)
  theta2 <- coev_calculate_theta(object, intervention_values = values2)
  delta_theta <-
    (theta2[,id_resp] - theta1[,id_resp]) / stats::median(diff[[id_resp]])
  # save as draws array for output
  delta_theta <- posterior::as_draws_array(as.matrix(delta_theta))
  dimnames(delta_theta)$variable <- "delta_theta"
  return( delta_theta )
}
