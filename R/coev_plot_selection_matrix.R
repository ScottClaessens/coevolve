#' Plot selection matrix from a fitted \code{coevfit} object
#'
#' Plot the selection matrix (i.e., the A matrix) for all trait pairs from a
#' fitted \code{coevfit} object. Values can be standardised for easier
#' interpretation.
#'
#' @param object An object of class \code{coevfit}
#' @param variables If NULL (default), the plot includes all coevolving
#'   variables from the model. Otherwise, a character vector of length >= 2
#'   declaring the variables to be included in the plot.
#' @param std Logical. If \code{FALSE} (default), values are not standardised.
#'   If \code{TRUE}, values are scaled by the standard deviation of each
#'   variable among taxa.
#' @param prob Probability mass to include in the inner interval. Default is
#'   0.66 (66% interval).
#' @param prob_outer Probability mass to include in the outer interval. Default
#'   is 0.95 (95% interval).
#'
#' @return A \code{ggplot} object
#'
#' @author Scott Claessens \email{scott.claessens@@gmail.com}, Erik Ringen
#'   \email{erikjacob.ringen@@uzh.ch}
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
#' # plot selection matrix
#' coev_plot_selection_matrix(fit)
#'
#' # standardise the cross-selection effects
#' coev_plot_selection_matrix(fit, std = TRUE)
#' }
#'
#' @export
coev_plot_selection_matrix <- function(object, variables = NULL, std = FALSE,
                                       prob = 0.66, prob_outer = 0.95) {
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
  # extract samples for A matrix
  post <- extract_samples(object)
  A <- post$A
  # if std = TRUE, scale by the standard deviation of each variable among taxa
  if (std) {
    N_tree <- object$stan_data$N_tree
    N_tips <- object$stan_data$N_tips
    eta_tips <- post$eta[, 1, 1:N_tips, ] # uses only first tree
    eta_sd <- apply(eta_tips, c(1, 3), stats::sd) # creates [samples, traits]
    A_temp <- A
    for (i in 1:object$stan_data$J) {
      for (j in 1:object$stan_data$J) {
        A_temp[, i, j] <- A[, i, j] * (eta_sd[, j] / eta_sd[, i])
      }
    }
    A <- A_temp
  }
  # get data for plot
  d <-
    tidyr::expand_grid(
      from = variables,
      to = variables
    ) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      A = list(A[, which(names(object$variables) == .data$to),
                 which(names(object$variables) == .data$from)])
    ) |>
    tidyr::unnest("A")
  # get quantiles
  dd <-
    d |>
    dplyr::group_by(.data$from, .data$to) |>
    dplyr::summarise(
      med = stats::median(.data$A),
      lower = stats::quantile(.data$A, 0.5 - (prob / 2)),
      upper = stats::quantile(.data$A, 0.5 + (prob / 2)),
      lower_outer = stats::quantile(.data$A, 0.5 - (prob_outer / 2)),
      upper_outer = stats::quantile(.data$A, 0.5 + (prob_outer / 2)),
      .groups = "drop"
    )
  # plot A matrix
  ggplot2::ggplot() +
    tidybayes::stat_halfeye(
      data = d,
      mapping = ggplot2::aes(x = .data$A),
      colour = NA,
      fill = "darkgrey",
      trim = TRUE,
      normalize = "panels"
    ) +
    ggplot2::facet_grid(
      as.formula("to ~ from"),
      scales = "free",
      switch = "y",
      axes = "all"
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
    ggplot2::labs(
      x = NULL,
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
}
