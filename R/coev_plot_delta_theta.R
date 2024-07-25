#' Plot delta theta values from a fitted \code{coevfit} object
#'
#' @param object An object of class \code{coevfit}
#' @param variables If NULL (default), the plot includes all coevolving
#'   variables from the model. Otherwise, a character vector of length >= 2
#'   declaring the variables to be included in the plot.
#' @param ... Additional arguments passed to \code{ggdist::stat_slabinterval}
#'
#' @return A \code{ggplot} object
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
#' # plot delta theta values for all effects
#' coev_plot_delta_theta(m)
#' }
coev_plot_delta_theta <- function(object, variables = NULL, ...) {
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
  # get range for plotting
  dd <- dplyr::group_by(d, .data$response, .data$predictor)
  dd <- dplyr::summarise(
    dd,
    lower = stats::quantile(.data$delta_theta, 0.05),
    upper = stats::quantile(.data$delta_theta, 0.95),
    .groups = "drop"
    )
  # plot
  ggplot2::ggplot(
    data = d,
    mapping = ggplot2::aes(x = .data$delta_theta)
    ) +
    ggdist::stat_slabinterval(
      .width = c(0.5, 0.89), # 50% and 89% credible intervals
      n = 1e4                # increased resolution
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
    ggplot2::coord_cartesian(
      xlim = c(min(dd$lower), max(dd$upper))
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
}
