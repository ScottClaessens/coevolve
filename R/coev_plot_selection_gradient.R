#' Plot selection gradient heatmap from a fitted \code{coevfit} object
#'
#' @param object An object of class \code{coevfit}
#' @param var1 A character string equal to one of the coevolving variables in the model
#' @param var2 A character string equal to one of the coevolving variables in the model
#' @param contour Logical (defaults to FALSE); whether to show white contour lines
#'   to indicate where selection is stronger than drift
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
#' # plot selection gradient
#' coev_plot_selection_gradient(m)
#' }
coev_plot_selection_gradient <- function(object, var1, var2, contour = FALSE) {
  # stop if object is not of class coevfit
  if (!methods::is(object, "coevfit")) {
    stop2("Argument 'object' must be a fitted coevolutionary model of class coevfit.")
  }
  if (!is.character(var1) | length(var1) != 1) {
    # stop if var1 not character string of length one
    stop2("Argument 'var1' must be a character string of length one.")
  } else if (!(var1 %in% names(object$variables))) {
    # stop if var1 not included in model
    stop2("Argument 'var1' must be a variable included in the fitted model.")
  }
  if (!is.character(var2) | length(var2) != 1) {
    # stop if var2 not character string of length one
    stop2("Argument 'var2' must be a character string of length one.")
  } else if (!(var2 %in% names(object$variables))) {
    # stop if var2 not included in model
    stop2("Argument 'var2' must be a variable included in the fitted model.")
  }
  # stop if var1 and var2 are the same variable
  if (var1 == var2) {
    stop2("Argument 'var1' and 'var2' must refer to different variables.")
  }
  # stop if contour not logical
  if (!is.logical(contour)) {
    stop2("Argument 'contour' must be logical.")
  }
  # get IDs for variables
  id_var1 <- which(names(object$variables) == var1)
  id_var2 <- which(names(object$variables) == var2)
  # get posterior draws
  draws <- posterior::as_draws_rvars(object$fit)
  # medians and median absolute deviations for both variables
  eta_var1 <- posterior::rvar_median(draws$eta[1:nrow(object$data), id_var1])
  eta_var2 <- posterior::rvar_median(draws$eta[1:nrow(object$data), id_var2])
  med_var1 <- stats::median(eta_var1)
  med_var2 <- stats::median(eta_var2)
  mad_var1 <- stats::mad(eta_var1)
  mad_var2 <- stats::mad(eta_var2)
  lower_var1 <- med_var1 - 2.5*mad_var1
  lower_var2 <- med_var2 - 2.5*mad_var2
  upper_var1 <- med_var1 + 2.5*mad_var1
  upper_var2 <- med_var2 + 2.5*mad_var2
  # get median parameter values for A, b, and Q_diag
  A <- stats::median(draws$A)
  b <- stats::median(draws$b)
  Q_diag <- stats::median(draws$Q_diag)
  # ornstein uhlenbeck sde function for response and predictor variable
  # this currently assumes that values for all other traits are set to zero
  OU_sde <- function(resp_value, pred_value, resp_id, pred_id) {
    # get median absolute deviation to scale by
    mad_scale <- ifelse(resp_id == id_var1, mad_var1, mad_var2)
    # autoregressive selection effect
    ((A[resp_id, resp_id] * resp_value +
        # cross-lagged selection effect
        A[resp_id, pred_id] * pred_value +
        # sde intercept
        b[resp_id]
      # scaled by mad
      ) / mad_scale) /
    # divided by sigma, scaled by mad
    (Q_diag[resp_id] / mad_scale)
  }
  # get predictions for different levels of traits
  preds <-
    expand.grid(
      var1_value = seq(
        from = lower_var1,
        to = upper_var1,
        length.out = 20
        ),
      var2_value = seq(
        from = lower_var2,
        to = upper_var2,
        length.out = 20
        ),
      var1_delta_sigma = NA,
      var2_delta_sigma = NA
    )
  for (i in 1:nrow(preds)) {
    preds$var1_delta_sigma[i] <-
      OU_sde(
        resp_value = preds$var1_value[i],
        pred_value = preds$var2_value[i],
        resp_id = id_var1,
        pred_id = id_var2
      )
    preds$var2_delta_sigma[i] <-
      OU_sde(
        resp_value = preds$var2_value[i],
        pred_value = preds$var1_value[i],
        resp_id = id_var2,
        pred_id = id_var1
      )
  }
  # plotting
  out <-
    tidyr::pivot_longer(
      preds,
      -c("var1_value", "var2_value")
      )
  out <-
    dplyr::mutate(
      out,
      name = factor(
        .data$name,
        labels = c(
          expression(paste(Delta, !!var1)),
          expression(paste(Delta, !!var2))
          )
        )
      )
  # plot as z-score
  out <-
    ggplot2::ggplot(
      data = out,
      mapping = ggplot2::aes(
        x = (.data$var1_value - med_var1) / mad_var1,
        y = (.data$var2_value - med_var2) / mad_var2
      )
    ) +
    ggplot2::facet_wrap(
      ~ .data$name,
      labeller = ggplot2::label_parsed
      ) +
    ggplot2::geom_raster(
      mapping = ggplot2::aes(
        fill = .data$value
      )
    )
  # add contour
  if (contour) {
    out <-
      out +
      ggplot2::stat_contour(
        mapping = ggplot2::aes(z = .data$value),
        colour = "white",
        breaks = c(-1, 1)
        )
  }
  out +
    ggplot2::scale_x_continuous(
      expand = c(0, 0)
      ) +
    ggplot2::scale_y_continuous(
      expand = c(0, 0)
      ) +
    colorspace::scale_fill_continuous_divergingx(
      palette = "Geyser",
      trans = "reverse",
      guide = ggplot2::guide_colourbar(reverse = TRUE)
      ) +
    ggplot2::theme_bw(
      base_size = 14
      ) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.spacing = ggplot2::unit(1.5, "lines"),
      strip.background = ggplot2::element_blank()
      ) +
    ggplot2::labs(
      fill = expression(frac(paste(Delta, alpha), sigma))
      ) +
    ggplot2::xlab(
      paste(var1, "(z-score)")
      ) +
    ggplot2::ylab(
      paste(var2, "(z-score)")
      ) +
    ggplot2::coord_fixed()
}
