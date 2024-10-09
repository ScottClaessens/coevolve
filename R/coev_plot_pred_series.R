#' Plot a predicted co-evolutionary time series from a fitted \code{coevfit}
#' object
#'
#' @param object An object of class \code{coevfit}
#' @param prob A value between 0 and 1 indicating the desired probability
#'   to be covered by the uncertainty intervals. The default is 0.95.
#' @param ... Additional arguments passed to \code{\link{coev_pred_series}}
#' @return A ggplot object.
#'
#' @author Scott Claessens \email{scott.claessens@@gmail.com}, Erik Ringen
#'   \email{erikjacob.ringen@@uzh.ch}
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
#'   stochastic = T
#'   )
#' 
#' # expected trait co-evolution, no drift
#' coev_plot_pred_series(
#'   object = fit,
#'   stochastic = F
#'   )
#' }
#'
#' @export
coev_plot_pred_series <- function(object, prob = 0.95, ...){
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
  combined_args <- modifyList(default_args, user_args)
  combined_args$object <- object
  preds <- do.call(coev_pred_series, combined_args)
  
  if (combined_args$stochastic == F) {
      # for CI quantiles
      probs <- c(((1 - prob) / 2), 1 - ((1 - prob) / 2))

      epreds_long <- preds |>
      cubelyr::as.tbl_cube(met_name = "est") |> 
      as.data.frame()

      epreds_summary <- epreds_long |> 
      dplyr::group_by(response, time) |> 
      dplyr::summarise(mean = mean(est),
      lower_CI = quantile(est, probs[1]),
      upper_CI = quantile(est, probs[2]),
      )
    
      p <- ggplot2::ggplot(epreds_summary, aes(x = time, y = mean, fill = response, color = response, linetype = response)) +
      geom_line(lwd = 1) + 
        geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), alpha = 0.25, lwd = 0.25) +
      theme_classic(base_size = 14) + 
      scale_x_continuous(breaks = c(0, max(epreds_summary$time)), labels = c("LCA", "present")) +
      ylab("trait value (latent scale)") + 
      xlab("time") +
      theme(legend.title = element_blank(),  strip.background = element_blank(),strip.text = element_blank()) + 
      ggtitle("Expected trait coevolution")
  }
  else if (combined_args$stochastic == T) {
    sims_long <- preds |>
      cubelyr::as.tbl_cube(met_name = "est") |> 
      as.data.frame() |> 
      dplyr::mutate(sim = factor(paste("sim", samps), levels = paste("sim", sort(unique(samps)))))

    p <- ggplot2::ggplot(sims_long |> dplyr::filter(samps <= 15), aes(x = time, y = est, color = response)) +
    facet_wrap(~sim, ncol = 5) +
    geom_line(lwd = 1) + 
    theme_minimal(base_size = 14) + 
    scale_x_continuous(breaks = c(0, max(sims_long$time)), labels = c("LCA", "present")) +
    ylab("trait value (latent scale)") + 
    xlab("time") +
    theme(legend.title = element_blank(),  strip.background = element_blank(), panel.spacing.x = unit(1, "lines"), axis.text.x = element_blank(),  axis.ticks.x = element_blank()) + 
    ggtitle("Predicted trait coevolution")
  }
  return(p)
}