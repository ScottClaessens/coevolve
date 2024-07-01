#' Plot selection gradient heatmap from a fitted \code{coevfit} object
#'
#' @param object An object of class \code{coevfit}
#' @param var1 A character string equal to one of the coevolving variables in the model
#' @param var2 A character string equal to one of the coevolving variables in the model
#' @param contour Logical (defaults to FALSE); whether to show white contour lines
#' to indicate where selection is stronger than drift
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
  # get posterior draws for eta
  suppressWarnings({eta <- tidybayes::gather_draws(object$fit, eta[node,variable])})
  # restrict to tips only (nodes 1-n)
  eta <- eta[eta$node %in% 1:nrow(object$data),]
  # calculate posterior median eta for each tip
  eta <- dplyr::summarise(eta, value = stats::median(.value), .groups = "drop")
  # calculate median, mad, high, and low across all tips
  eta <- dplyr::group_by(eta, variable)
  eta <- dplyr::summarise(
    eta,
    median = stats::median(value),
    mad    = stats::mad(value),
    high   = median + 2.5*mad,
    low    = median - 2.5*mad
  )
  # get median parameter values for A, b, and Q_diag
  A <- apply(object$fit$draws("A"), 3, stats::median)
  dim(A) <- rep(length(names(object$variables)), 2)
  b <- as.vector(apply(object$fit$draws("b"), 3, stats::median))
  Q_diag <- as.vector(apply(object$fit$draws("Q_diag"), 3, stats::median))
  # ornstein uhlenbeck sde function for response and predictor variable
  # this currently assumes that values for all other traits are set to zero
  OU_sde <- function(respValue, predValue, respVarNum, predVarNum) {
    delta_sigma <-
      # autoregressive selection effect
      ((A[respVarNum,respVarNum] * respValue +
          # cross-lagged selection effect
          A[respVarNum,predVarNum] * predValue +
          # sde intercept
          b[respVarNum]
        # scaled by mad
        ) / eta$mad[respVarNum]) /
      # divided by sigma, scaled by mad
      (Q_diag[respVarNum] / eta$mad[respVarNum])
    return(delta_sigma)
  }
  # get predictions for different levels of traits
  var1Num <- which(var1 == names(object$variables))
  var2Num <- which(var2 == names(object$variables))
  preds <-
    expand.grid(
      var1Value = seq(
        from = eta$low[eta$variable == var1Num],
        to = eta$high[eta$variable == var1Num],
        length.out = 20
        ),
      var2Value = seq(
        from = eta$low[eta$variable == var2Num],
        to = eta$high[eta$variable == var2Num],
        length.out = 20
        ),
      var1DeltaSigma = NA,
      var2DeltaSigma = NA
    )
  for (i in 1:nrow(preds)) {
    preds$var1DeltaSigma[i] <-
      OU_sde(
        resp = preds$var1Value[i],
        pred = preds$var2Value[i],
        respVarNum = var1Num,
        predVarNum = var2Num
      )
    preds$var2DeltaSigma[i] <-
      OU_sde(
        resp = preds$var2Value[i],
        pred = preds$var1Value[i],
        respVarNum = var2Num,
        predVarNum = var1Num
      )
  }
  # plotting
  out <- tidyr::pivot_longer(preds, -c(var1Value, var2Value))
  out <- dplyr::mutate(
    out,
    name = factor(
      name,
      labels = c(
        expression(paste(Delta, !!var1)),
        expression(paste(Delta, !!var2))
        )
      )
    )
  # plot as z-score
  out <- ggplot2::ggplot(
    data = out,
    ggplot2::aes(x = (var1Value - eta$median[var1Num]) / eta$mad[var1Num],
                 y = (var2Value - eta$median[var2Num]) / eta$mad[var2Num],
                 z = value, fill = value)
    ) +
    ggplot2::facet_wrap(~ name, labeller = ggplot2::label_parsed) +
    ggplot2::geom_raster()
  # add contour
  if (contour) {
    out <- out + ggplot2::geom_contour(colour = "white", breaks = c(-1, 1))
  }
  out +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    colorspace::scale_fill_continuous_divergingx(
      palette = "Geyser", trans = "reverse",
      guide = ggplot2::guide_colourbar(reverse = TRUE)
      ) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.spacing = ggplot2::unit(1.5, "lines"),
                   strip.background = ggplot2::element_blank()) +
    ggplot2::labs(fill = expression(frac(paste(Delta, alpha), sigma))) +
    ggplot2::xlab(paste(var1, "(z-score)")) +
    ggplot2::ylab(paste(var2, "(z-score)")) +
    ggplot2::coord_fixed()
}
