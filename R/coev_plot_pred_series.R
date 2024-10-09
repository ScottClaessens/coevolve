#' Plot a predicted coevolutionary time series from a fitted \code{coevfit}
#' object
#'
#' This function plots a predicted coevolutionary time series using the
#' estimated parameters from a fitted \code{coevfit} model. By default, the
#' plot uses the posterior ancestral states estimated by the model as the
#' starting values, but users can also set their own starting values for
#' traits. Plots can be generated with or without stochastic drift. For more
#' details on the underlying predictive function, see
#' \code{help(coev_pred_series)}.
#'
#' @param object An object of class \code{coevfit}
#' @param prob A value between 0 and 1 indicating the desired probability
#'   to be covered by the uncertainty intervals. The default is 0.95.
#' @param ... Additional arguments passed to \code{\link{coev_pred_series}}
#'
#' @return A ggplot object.
#'
#' @author Scott Claessens \email{scott.claessens@@gmail.com}, Erik Ringen
#'   \email{erikjacob.ringen@@uzh.ch}
#'
#' @seealso \code{\link{coev_pred_series}}
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
#' # simulated trait co-evolution
#' coev_plot_pred_series(
#'   object = fit,
#'   stochastic = TRUE
#'   )
#'
#' # expected trait co-evolution, no drift
#' coev_plot_pred_series(
#'   object = fit,
#'   stochastic = FALSE
#'   )
#' }
#'
#' @export
coev_plot_pred_series <- function(object, prob = 0.95, ...) {
  # stop if object not of class coevfit
  if (!methods::is(object, "coevfit")) {
    stop2(
      paste0(
        "Argument 'object' must be a fitted coevolutionary model ",
        "of class coevfit."
        )
      )
  }
  # stop if prob is outside of range 0 - 1
  if (prob <= 0 | prob >= 1) {
    stop2("Argument 'prob' is not between 0 and 1.")
  }
  # user-supplied arguments
  user_args <- list(...)
  # default arguments from coev_pred_series
  default_args <- as.list(formals(coev_pred_series))
  combined_args <- utils::modifyList(default_args, user_args)
  combined_args$object <- object
  preds <- do.call(coev_pred_series, combined_args)
  # create plot
  if (combined_args$stochastic == FALSE) {
      # get CI quantiles
      probs <- c(((1 - prob) / 2), 1 - ((1 - prob) / 2))
      # get expected values in cube format
      epreds_long <-
        preds |>
        cubelyr::as.tbl_cube(
          met_name = "est"
          ) |>
        as.data.frame()
      # summarise expected values
      epreds_summary <-
        epreds_long |>
        dplyr::group_by(.data$response, .data$time) |>
        dplyr::summarise(
          mean = mean(.data$est),
          lower_CI = stats::quantile(.data$est, probs[1]),
          upper_CI = stats::quantile(.data$est, probs[2])
      )
      # plot
      p <-
        ggplot2::ggplot(
          epreds_summary,
          ggplot2::aes(
            x = .data$time,
            y = .data$mean,
            fill = .data$response,
            color = .data$response,
            linetype = .data$response
            )
          ) +
      ggplot2::geom_line(lwd = 1) +
      ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin = .data$lower_CI,
          ymax = .data$upper_CI
          ),
        alpha = 0.25,
        lwd = 0.25
        ) +
      ggplot2::theme_classic(base_size = 14) +
      ggplot2::scale_x_continuous(
        breaks = c(0, max(epreds_summary$time)),
        labels = c("LCA", "present")
        ) +
      ggplot2::ylab("trait value (latent scale)") +
      ggplot2::xlab("time") +
      ggplot2::theme(
        legend.title = ggplot2::element_blank(),
        strip.background = ggplot2::element_blank(),
        strip.text = ggplot2::element_blank()
        ) +
      ggplot2::ggtitle("Expected trait coevolution")
  } else if (combined_args$stochastic == TRUE) {
    # get predictions
    sims_long <-
      preds |>
      cubelyr::as.tbl_cube(met_name = "est") |>
      as.data.frame() |>
      dplyr::mutate(
        sim = factor(
          paste("sim", .data$samps),
          levels = paste("sim", sort(unique(.data$samps)))
          )
        )
    # plot
    p <-
      sims_long |>
      dplyr::filter(.data$samps <= 15) |>
      ggplot2::ggplot(
        ggplot2::aes(
          x = .data$time,
          y = .data$est,
          color = .data$response
          )
        ) +
    ggplot2::facet_wrap( ~ .data$sim, ncol = 5) +
    ggplot2::geom_line(lwd = 1) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::scale_x_continuous(
      breaks = c(0, max(sims_long$time)),
      labels = c("LCA", "present")
      ) +
    ggplot2::ylab("trait value (latent scale)") +
    ggplot2::xlab("time") +
    ggplot2::theme(
      legend.title = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      panel.spacing.x = ggplot2::unit(1, "lines"),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
      ) +
    ggplot2::ggtitle("Predicted trait coevolution")
  }
  return(p)
}
