#' Trace and density plots for \code{coevfit} objects
#'
#' @srrstats {G1.3, G1.4, G2.1a} Function documentation begins here, with
#'   expected data types and definitions of statistical terminology and inputs
#' @srrstats {G2.3, G2.3b} Documenting that character parameters are
#'   strictly case-sensitive
#' @srrstats {BS6.1, BS6.2, BS6.3, BS6.5} This plot method displays
#'   distributional estimates (densities) and sequences of samples (trace plots)
#'   from the model object
#'
#' @param object An object of class \code{coevfit}.
#' @param parameters If NULL (default), the function returns a list of plots for
#'   the main parameters from the model (the selection matrix "A", the drift
#'   matrix "Q", and the continuous time intercepts "b"). Otherwise, a character
#'   vector declaring the parameters to plot (strictly case-sensitive).
#' @param combo A character vector with at least two elements. Each element of
#'   \code{combo} corresponds to a column in the resulting plot and should be
#'   the name of one of the available
#'   \code{\link[bayesplot:MCMC-overview]{MCMC}} functions (omitting the
#'   \code{mcmc_} prefix).
#' @param npars Number of parameters displayed per page (default = 5).
#' @param plot Logical. If TRUE (default), plots are plotted in the active
#'   graphic device. If FALSE, plots are not plotted.
#'
#' @return A list of \code{ggplot} objects
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
#' # print trace and density plots
#' plot(fit)
#' }
#'
#' @method plot coevfit
#' @export
plot.coevfit <- function(object, parameters = NULL, combo = c("dens", "trace"),
                         npars = 5, plot = TRUE) {
  if (!is.null(parameters)) {
    if (!all(parameters %in% object$fit$metadata()$model_params)) {
      stop2("Argument 'parameters' contains invalid parameter names.")
    }
  } else {
    # keep only main parameters: A, Q, and b
    parameters <- object$fit$metadata()$model_params
    parameters <- parameters[
      stringr::str_starts(parameters, stringr::fixed("A[")) |
        stringr::str_starts(parameters, stringr::fixed("Q[")) |
        stringr::str_starts(parameters, stringr::fixed("b["))
    ]
  }
  if (!is.numeric(npars) || npars <= 0) {
    stop2("Argument 'npars' is not a positive number.")
  }
  if (!is.logical(plot)) {
    stop2("Argument 'plot' is not logical.")
  }
  n_plots <- ceiling(length(parameters) / npars)
  plots <- vector(mode = "list", length = n_plots)
  for (i in seq_len(n_plots)) {
    sub <- ((i - 1) * npars + 1):min(i * npars, length(parameters))
    sub_pars <- parameters[sub]
    plots[[i]] <- bayesplot::mcmc_combo(
      posterior::as_draws(
        object$fit,
        variable = sub_pars
      ),
      combo = combo
    )
    if (plot) {
      plot(plots[[i]], newpage = TRUE)
    }
  }
  invisible(plots)
}
