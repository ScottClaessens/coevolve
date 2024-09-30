#' Plot delta theta values from a fitted \code{coevfit} object
#'
#' Plot delta theta values for all trait pairs from a fitted \code{coevfit}
#' object. This plot can be used to visually assess contingencies and
#' directionality between different variables in the coevolutionary process.
#'
#' @param object An object of class \code{coevfit}
#' @param variables If NULL (default), the plot includes all coevolving
#'   variables from the model. Otherwise, a character vector of length >= 2
#'   declaring the variables to be included in the plot.
#' @param prob Probability mass to include in the inner interval. Default is
#'   0.66 (66% interval).
#' @param prob_outer Probability mass to include in the outer interval. Default
#'   is 0.95 (95% interval).
#' @param limits If NULL (default), limits are scaled automatically to include
#'   all posterior samples. Otherwise, a numeric vector of length 2 specifying
#'   the lower and upper limits for the x-axis.
#'
#' @return A \code{ggplot} object
#'
#' @author Scott Claessens \email{scott.claessens@@gmail.com}, Erik Ringen
#'   \email{erikjacob.ringen@@uzh.ch}
#'
#' @details This function repeatedly uses the
#'   \code{\link{coev_calculate_delta_theta}} function under the hood to
#'   generate a pairs plot of \eqn{\Delta\theta} for all variables in the model.
#'   For more details on the definition and calculation of \eqn{\Delta\theta},
#'   see \code{help(coev_calculate_delta_theta)}. Note that often the posterior
#'   distribution for \eqn{\Delta\theta} has long tails, meaning that the
#'   distribution for different traits can be difficult to visualise in a single
#'   pairs plot. If this plot does not produce satisfactory visualisations, the
#'   user should either specify narrower limits for the x-axis or use the
#'   \code{\link{coev_calculate_delta_theta}} function to create plots manually.
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
#'   \code{\link{coev_calculate_theta}}
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
#' # plot delta theta values for all effects
#' coev_plot_delta_theta(fit)
#' }
#'
#' @export
coev_plot_delta_theta <- function(object, variables = NULL, prob = 0.66,
                                  prob_outer = 0.95, limits = NULL) {
  # stop if object is not of class coevfit
  if (!methods::is(object, "coevfit")) {
    stop2(
      paste0(
        "Argument 'object' must be a fitted coevolutionary model ",
        "of class coevfit."
        )
      )
  }
  # if user specifies variables argument:
  if (!is.null(variables)) {
    if (!methods::is(variables, "character")) {
      # stop if variables not character vector
      stop2("Argument 'variables' must be a character vector.")
    } else if (!(length(variables) > 1)) {
      # stop if variables not > length 1
      stop2("Argument 'variables' must be of length > 1.")
    } else if (!all(variables %in% names(object$variables))) {
      # stop if variables not included in model
      stop2(
        "Some variables in 'variables' are not included in the fitted model."
        )
    } else if (any(duplicated(variables))) {
      # stop if variables contains duplicates
      stop2("Argument 'variables' contains duplicates.")
    }
  } else {
    # otherwise, default is to include all variables in order
    variables <- names(object$variables)
  }
  if (!methods::is(prob, "numeric")) {
    # stop if prob not numeric
    stop2("Argument 'prob' must be numeric.")
  } else if (length(prob) != 1) {
    # stop if prob not == length 1
    stop2("Argument 'prob' must be of length 1.")
  } else if (prob <= 0 | prob >= 1) {
    # stop if prob not between 0 and 1
    stop2("Argument 'prob' must be between 0 and 1.")
  }
  if (!methods::is(prob_outer, "numeric")) {
    # stop if prob_outer not numeric
    stop2("Argument 'prob_outer' must be numeric.")
  } else if (length(prob_outer) != 1) {
    # stop if prob_outer not == length 1
    stop2("Argument 'prob_outer' must be of length 1.")
  } else if (prob_outer <= 0 | prob_outer >= 1) {
    # stop if prob_outer not between 0 and 1
    stop2("Argument 'prob_outer' must be between 0 and 1.")
  } else if (prob_outer <= prob) {
    # stop if prob is greater than or equal to prob_outer
    stop2("Argument 'prob_outer' must be greater than argument 'prob'.")
  }
  # if user specifies limits argument
  if (!is.null(limits)) {
    if (!methods::is(limits, "numeric")) {
      # stop if limits not numeric vector
      stop2("Argument 'limits' must be a numeric vector.")
    } else if (!(length(limits) == 2)) {
      # stop if limits not == length 2
      stop2("Argument 'limits' must be of length 2.")
    }
  }
  # prepare data for plot
  d <- tidyr::expand_grid(
    response = variables,
    predictor = variables
    )
  # remove autocorrelation effects
  d <- d[d$response != d$predictor,]
  # for each combination, calculate delta theta
  d <- dplyr::mutate(
    d,
    delta_theta = purrr::map2(
      .x = .data$response,
      .y = .data$predictor,
      .f = function(x, y) as.numeric(
        coev_calculate_delta_theta(
          object,
          response = x,
          predictor = y
        )
      )
    )
  )
  # unlist result
  d <- tidyr::unnest(d, "delta_theta")
  # response and predictor as factors
  d$response  <- factor(d$response, levels = variables)
  d$predictor <- factor(d$predictor, levels = variables)
  # get median and lower/upper
  dd <- dplyr::group_by(d, .data$response, .data$predictor)
  dd <- dplyr::summarise(
    dd,
    med = stats::median(.data$delta_theta),
    lower = stats::quantile(.data$delta_theta, 0.5 - (prob / 2)),
    upper = stats::quantile(.data$delta_theta, 0.5 + (prob / 2)),
    lower_outer = stats::quantile(.data$delta_theta, 0.5 - (prob_outer / 2)),
    upper_outer = stats::quantile(.data$delta_theta, 0.5 + (prob_outer / 2)),
    .groups = "drop"
  )
  # plot
  p <-
    ggplot2::ggplot() +
    ggplot2::geom_density(
      data = d,
      mapping = ggplot2::aes(x = .data$delta_theta),
      colour = NA,
      fill = "darkgrey"
    ) +
    ggplot2::geom_linerange(
      data = dd,
      mapping = ggplot2::aes(
        y = 0,
        xmin = .data$lower_outer,
        xmax = .data$upper_outer
      )
    ) +
    ggplot2::geom_pointrange(
      data = dd,
      mapping = ggplot2::aes(
        x = .data$med,
        y = 0,
        xmin = .data$lower,
        xmax = .data$upper
        ),
      linewidth = 1
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed"
      ) +
    ggplot2::geom_rect(
      data = data.frame(
        response  = factor(variables, levels = variables),
        predictor = factor(variables, levels = variables),
        delta_theta = 0
        ),
      fill = "grey95",
      xmin = -Inf, xmax = Inf,
      ymin = -Inf, ymax = Inf
    ) +
    ggplot2::facet_grid(
      .data$predictor ~ .data$response,
      switch = "y"
      ) +
    ggplot2::labs(
      x = expression(paste(Delta, theta[z])),
      y = "From this variable...",
      title = "... to this variable."
      ) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect("grey"),
      panel.background = ggplot2::element_rect("white"),
      axis.line.x = ggplot2::element_line("black"),
      axis.title.x = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 11)
    )
  # if limits specified by user
  if (!is.null(limits)) p <- p + ggplot2::xlim(limits)
  return(p)
}
