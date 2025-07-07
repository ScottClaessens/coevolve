#' Calculate delta theta from a fitted \code{coevfit} object
#'
#' Calculate \eqn{\Delta\theta} from a fitted \code{coevfit} object.
#' \eqn{\Delta\theta} is defined as the difference in the equilibrium value for
#' a "response" trait that results from a standardised increase in the value of
#' another "predictor" trait. This value can be used to assess contingencies
#' and directionality between variables in the coevolutionary process.
#'
#' @srrstats {G1.3, G1.4, G2.1a} Function documentation begins here, with
#'   expected data types and definitions of statistical terminology and inputs
#' @srrstats {G2.0a} Secondary documentation on expected argument length
#' @srrstats {G2.3, G2.3b} Documenting that character parameters are
#'   strictly case-sensitive
#'
#' @param object An object of class \code{coevfit}
#' @param response A character string of length one equal to one of the
#'   coevolving variables in the model (strictly case-sensitive)
#' @param predictor A character string of length one equal to one of the
#'   coevolving variables in the model (strictly case-sensitive)
#'
#' @returns Posterior samples in the draws_array format
#'
#' @author Scott Claessens \email{scott.claessens@@gmail.com}, Erik Ringen
#'   \email{erikjacob.ringen@@uzh.ch}
#'
#' @details This function calculates \eqn{\Delta\theta}, which is defined as the
#'   difference in the equilibrium value for a "response" trait that results
#'   from a standardised increase in the value of another "predictor" trait. The
#'   function first calculates the equilibrium trait value for the response
#'   trait when the predictor trait is held at its empirical median value (the
#'   \code{\link{coev_calculate_theta}} function is used for this purpose, see
#'   \code{help(coev_calculate_theta)} for further details). The function then
#'   calculates the equilibrium trait value for the response variable after
#'   increasing the predictor trait from its median by one median absolute
#'   deviation. The function then returns the posterior difference between these
#'   values. The resulting \eqn{\Delta\theta} samples can be used to infer
#'   whether increases in one trait have a positive or negative selective effect
#'   on another trait in the model.
#'
#' @references
#' Ringen, E., Martin, J. S., & Jaeggi, A. (2021). Novel phylogenetic methods
#' reveal that resource-use intensification drives the evolution of "complex"
#' societies. \emph{EcoEvoRXiv}. \code{doi:10.32942/osf.io/wfp95}
#'
#' Sheehan, O., Watts, J., Gray, R. D., Bulbulia, J., Claessens, S., Ringen,
#' E. J., & Atkinson, Q. D. (2023). Coevolution of religious and political
#' authority in Austronesian societies. \emph{Nature Human Behaviour},
#' \emph{7}(1), 38-45. \code{10.1038/s41562-022-01471-y}
#'
#' @seealso \code{\link{coev_calculate_delta_theta}},
#'   \code{\link{coev_plot_delta_theta}}
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
#'
#' @export
coev_calculate_delta_theta <- function(object, response, predictor) {
  #' @srrstats {G5.2, G5.2a} Unique error messages for each input
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
  if (!is.character(response) || length(response) != 1) {
    # stop if response not character string of length one
    #' @srrstats {G2.0, G2.1, G2.2} Assertion on type and length of input
    stop2("Argument 'response' must be a character string of length one.")
  } else if (!(response %in% names(object$variables))) {
    # stop if response not included in model
    #' @srrstats {G2.3, G2.3a} Permit only expected character input
    stop2(
      "Argument 'response' must be a variable included in the fitted model."
    )
  }
  if (!is.character(predictor) || length(predictor) != 1) {
    # stop if predictor not character string of length one
    #' @srrstats {G2.0, G2.1, G2.2} Assertion on type and length of input
    stop2("Argument 'predictor' must be a character string of length one.")
  } else if (!(predictor %in% names(object$variables))) {
    # stop if predictor not included in model
    #' @srrstats {G2.3, G2.3a} Permit only expected character input
    stop2(
      "Argument 'predictor' must be a variable included in the fitted model."
    )
  }
  # stop if response and predictor are the same variable
  #' @srrstats {G2.3, G2.3a} Permit only expected character input
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
  eta  <- draws$eta[, 1:object$stan_data$N_tips, ]
  med  <- apply(eta, 3, posterior::rvar_median)
  diff <- apply(eta, 3, posterior::rvar_mad)
  # construct intervention values list
  values1 <- list()
  values2 <- list()
  for (j in seq_along(names(object$variables))) {
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
    (theta2[, id_resp] - theta1[, id_resp]) / stats::median(diff[[id_resp]])
  # save as draws array for output
  delta_theta <- posterior::as_draws_array(as.matrix(delta_theta))
  dimnames(delta_theta)$variable <- "delta_theta"
  delta_theta
}
