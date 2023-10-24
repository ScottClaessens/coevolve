#' Plot cross-selection effects from a fitted \code{coevfit} object
#'
#' @param object An object of class \code{coevfit}
#' @param ... Additional arguments passed to \code{tidybayes::stat_slabinterval}
#'
#' @return A \code{ggplot} object
#' @export
#'
#' @examples
#' \dontrun{
#' # simulate data
#' n <- 20
#' tree <- ape::rtree(n)
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
#' # plot cross selection effects
#' coev_plot_cross(m)
#' }
coev_plot_cross <- function(object, ...) {
  # stop if object is not of class coevfit
  if (!methods::is(object, "coevfit")) {
    stop2("Argument 'object' must be a fitted coevolutionary model of class coevfit.")
  }
  # extract posterior draws for alpha matrix
  suppressWarnings({
    post <- tidyr::pivot_longer(
      object$fit$draws("alpha", format = "df"),
      cols = tidyselect::starts_with("alpha"),
      names_to = "parameter"
    )
  })
  # add from and to variables
  post$from <- names(object$variables)[as.numeric(substr(post$parameter, 9, 9))]
  post$to   <- names(object$variables)[as.numeric(substr(post$parameter, 7, 7))]
  # identify auto and cross effects
  post$type <- ifelse(post$from == post$to, "auto", "cross")
  # plot
  ggplot2::ggplot(data = post, mapping = ggplot2::aes(x = .data$value)) +
    ggdist::stat_slabinterval() +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed"
      ) +
    ggplot2::geom_rect(
      data = post[post$type == "auto",],
      fill = "grey95",
      xmin = -Inf, xmax = Inf,
      ymin = -Inf, ymax = Inf
      ) +
    ggplot2::facet_grid(
      from ~ to,
      switch = "y"
      ) +
    ggplot2::labs(
      x = "Cross selection effect",
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
