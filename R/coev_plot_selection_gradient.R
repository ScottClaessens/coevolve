#' Plot selection gradient heatmap from a fitted \code{coevfit} object
#'
#' Plot a heatmap of the selection gradient for two variables from a fitted
#' \code{coevfit} object.
#'
#' @param object An object of class \code{coevfit}
#' @param var1 A character string equal to one of the coevolving variables in
#'   the model
#' @param var2 A character string equal to one of the coevolving variables in
#'   the model
#' @param contour Logical (defaults to FALSE); whether to show white contour
#'   lines to indicate where selection is stronger than drift
#' @param limits A numeric vector of length 2 (defaults to \code{c(-2.5, 2.5)});
#'   specifying the lower limit and the upper limit of the x and y axes.
#'
#' @return A \code{ggplot} object
#'
#' @author Scott Claessens \email{scott.claessens@@gmail.com}, Erik Ringen
#'   \email{erikjacob.ringen@@uzh.ch}
#'
#' @details The selection gradient is operationalised as the ratio of the
#'   change in a trait due to deterministic selection \eqn{\Delta\alpha} to
#'   the change in a trait due to stochastic drift \eqn{\sigma}. Values
#'   between -1 and 1 indicate parameter space where the change due to
#'   drift is greater than change due to selection on the trait. Conversely,
#'   values greater than 1 (or less than -1) indicate parameter space where
#'   positive (or negative) selection is stronger than drift. If three or more
#'   traits were included in the model, other traits are held at their median
#'   values during these computations.
#'
#'   If the plot does not look right, the user might try zooming out from the
#'   default parameter space by setting wider limits. For some variables (e.g.,
#'   continuous and count variables), the default limits may not be suitable.
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
#' @seealso \code{\link{coev_plot_flowfield}}
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
#' # plot selection gradient
#' coev_plot_selection_gradient(
#'   object = fit,
#'   var1 = "political_authority",
#'   var2 = "religious_authority"
#'   )
#' }
#'
#' @export
coev_plot_selection_gradient <- function(object, var1, var2,
                                         contour = FALSE,
                                         limits = c(-2.5, 2.5)) {
  # stop if object is not of class coevfit
  if (!methods::is(object, "coevfit")) {
    stop2(
      paste0(
        "Argument 'object' must be a fitted coevolutionary model of class ",
        "coevfit."
        )
      )
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
  # stop if limits is not a numeric vector of length 2
  if (!(is.numeric(limits) & is.vector(limits) & length(limits) == 2)) {
    stop2("Argument 'limits' must be a numeric vector of length 2.")
  }
  # get IDs for variables
  id_var1 <- which(names(object$variables) == var1)
  id_var2 <- which(names(object$variables) == var2)
  # get posterior draws
  draws <- posterior::as_draws_rvars(object$fit)
  # medians and median absolute deviations for all variables
  eta  <- apply(
    draws$eta[1:object$stan_data$N_tips,], 2, posterior::rvar_median
    )
  meds <- unlist(lapply(eta, stats::median))
  mads <- unlist(lapply(eta, stats::mad))
  lowers <- meds + limits[1]*mads
  uppers <- meds + limits[2]*mads
  # get median parameter values for A, b, and Q_diag
  A <- stats::median(draws$A)
  b <- stats::median(draws$b)
  Q_diag <- stats::median(draws$Q_diag)
  # ornstein uhlenbeck sde function for response and predictor variable
  OU_sde <- function(resp_value, pred_value, resp_id, pred_id) {
    # sde intercept
    out <- b[resp_id]
    # selection effects
    for (j in 1:length(names(object$variables))) {
      if (j == resp_id) {
        # autoregressive selection effect (response -> response)
        out <- out + A[resp_id,j] * resp_value
      } else if (j == pred_id) {
        # cross-lagged selection effect for predictor -> response
        out <- out + A[resp_id,j] * pred_value
      } else {
        # cross-lagged selection effects for any remaining variables
        # held at their median trait values
        out <- out + A[resp_id,j] * meds[j]
      }
    }
    # scale by mad for response variable
    out <- out / mads[resp_id]
    # divide by sigma scaled by mad for response variable
    sigma <- Q_diag[resp_id] / mads[resp_id]
    out <- out / sigma
    return(out)
  }
  # get predictions for different levels of traits
  preds <-
    expand.grid(
      var1_value = seq(
        from = lowers[id_var1],
        to = uppers[id_var1],
        length.out = 20
        ),
      var2_value = seq(
        from = lowers[id_var2],
        to = uppers[id_var2],
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
        x = (.data$var1_value - meds[id_var1]) / mads[id_var1],
        y = (.data$var2_value - meds[id_var2]) / mads[id_var2]
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
