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
#' @author Scott Claessens \email{scott.claessens@@gmail.com}, 
#'   Erik Ringen \email{erikjacob.ringen@@uzh.ch}
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
#'   stochastic = FALSE,
#'   eta_anc = list(political_authority = -2, religious_authority = 1.5)
#'   )
#' }
#'
#' @export
coev_plot_pred_series <- function(object, prob = 0.95, ...) {
  # Check if object is of class coevfit
  if (!methods::is(object, "coevfit")) {
    stop2("Argument 'object' must be a fitted coevolutionary model of class 'coevfit'.")
  }
  # Check if prob is between 0 and 1
  if (!is.numeric(prob) || length(prob) != 1 || prob <= 0 || prob >= 1) {
    stop2("Argument 'prob' must be a single numeric value between 0 and 1.")
  }
  # Collect additional arguments
  user_args <- list(...)
  # Define default arguments from coev_pred_series, excluding 'object'
  default_args <- formals(coev_pred_series)
  default_args <- default_args[names(default_args) != "object"]
  # Convert default_args to a list with their default values
  default_args <- as.list(default_args)
  # Combine user arguments with defaults, giving precedence to user arguments
  combined_args <- modifyList(default_args, user_args)
  # Ensure 'object' is set
  combined_args$object <- object
  # Call coev_pred_series with combined arguments
  preds <- do.call(coev_pred_series, combined_args)
  # Find held indices, if any
  held_indices <- which(!is.na(unlist(combined_args$intervention_values[names(combined_args$intervention_values)])))
  # Check if 'preds' has the necessary dimension names
  if (is.null(dimnames(preds)) || 
      !all(c("samps", "time", "response") %in% names(dimnames(preds)))) {
    stop2("The 'preds' array must have dimension names: 'samps', 'time', 'response'.")
  }
  # Convert the 3D array 'preds' to a long-format data frame using tidyr::pivot_longer
  preds_long <- preds |>
    as.data.frame.table(responseName = "est") |>
    tidyr::pivot_longer(
      cols = -c(samps, time, response),
      names_to = NULL,
      values_to = "est"
    ) |>
    dplyr::mutate(
      samps = as.numeric(as.character(samps)),
      time = as.numeric(as.character(time))
    )
  # For deterministic predictions
  if (combined_args$stochastic == FALSE) {
    # Calculate confidence interval bounds
    probs <- c((1 - prob) / 2, 1 - (1 - prob) / 2)
    # Summarize mean and confidence intervals for each response and time
    epreds_summary <- preds_long |>
      dplyr::group_by(response, time) |>
      dplyr::summarise(
        mean = mean(est, na.rm = TRUE),
        lower_CI = stats::quantile(est, probs[1], na.rm = TRUE),
        upper_CI = stats::quantile(est, probs[2], na.rm = TRUE),
        .groups = "drop"
      )
    title <- "Expected trait coevolution"
    # label which traits are held
    if (length(held_indices) > 0) {
      epreds_summary$response <- ifelse(epreds_summary$response %in% names(held_indices), paste(epreds_summary$response, "(held)"), paste(epreds_summary$response, "(free)"))
      title <- paste(title, "| intervention")
    }
    # Create the deterministic plot
    p <- ggplot2::ggplot(epreds_summary, ggplot2::aes(x = time, y = mean, color = response, fill = response, linetype = response)) +
      ggplot2::geom_line(size = 1) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower_CI, ymax = upper_CI), alpha = 0.25, color = NA) +
      ggplot2::theme_classic(base_size = 14) +
      ggplot2::scale_x_continuous(
        breaks = c(1, max(epreds_summary$time)),
        labels = c("LCA", "Present")
      ) +
      ggplot2::ylab("Trait Value (Latent Scale)") +
      ggplot2::xlab("Time") +
      ggplot2::theme(
        legend.title = ggplot2::element_blank(),
        strip.background = ggplot2::element_blank(),
        strip.text = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(size = 14)
      ) +
      ggplot2::ggtitle(title)
  } 
  # For stochastic predictions
  else if (combined_args$stochastic == TRUE) {
    # Limit to first 15 simulations for clarity
    max_sims <- min(15, max(preds_long$samps))
    sampled_sim <- sample(1:max(preds_long$samps), size = max_sims, replace = F)

    sims_long <- preds_long |>
      dplyr::filter(samps %in% sampled_sim) |>
      dplyr::mutate(
        sim = factor(paste("Sim", match(samps, unique(samps))))
      )
    sim_num <- as.numeric(gsub("Sim ", "", unique(sims_long$sim)))
    sims_long$sim <- factor(sims_long$sim , levels = unique(sims_long$sim)[order(sim_num)])
    # Create the stochastic plot
    p <- ggplot2::ggplot(sims_long, ggplot2::aes(x = time, y = est, color = response, linetype = response)) +
      ggplot2::geom_line(size = 1) +
      ggplot2::facet_wrap(~ sim, ncol = 5) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::scale_x_continuous(
        breaks = c(1, max(sims_long$time)),
        labels = c("LCA", "Present")
      ) +
      ggplot2::ylab("Trait Value (Latent Scale)") +
      ggplot2::xlab("Time") +
      ggplot2::theme(
        legend.title = ggplot2::element_blank(),
        strip.background = ggplot2::element_blank(),
        panel.spacing = ggplot2::unit(1, "lines"),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(size = 14)
      ) +
      ggplot2::ggtitle("Predicted trait coevolution (stochastic)")
  } else {
    stop2("Invalid value for 'stochastic'. Must be TRUE or FALSE.")
  }
  return(p)
}