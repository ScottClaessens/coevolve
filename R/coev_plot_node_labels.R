#' Plot phylogenetic tree with node labels from a fitted coevolutionary model
#'
#' Plots the phylogenetic tree from a fitted \code{coevfit} object with
#' internal node IDs displayed, making it easy to identify node numbers
#' for use with \code{\link{coev_ancestral_states}}. The tree is returned
#' invisibly so it can be passed to \code{\link[ape]{identify.phylo}} for
#' interactive node selection.
#'
#' @param object An object of class \code{coevfit}
#' @param tip_labels Logical. If \code{FALSE} (default), tip labels are
#'   hidden to reduce clutter. Set to \code{TRUE} to display tip labels.
#' @param node_cex Numeric scaling factor for node label size. If \code{NULL}
#'   (default), the size is chosen automatically based on the number of tips.
#' @param node_col Color for node labels. Default is \code{"#B2182B"} (red).
#' @param ... Additional arguments passed to \code{\link[ape]{plot.phylo}}
#'   (e.g. \code{cex} for tip label size, \code{edge.width}, \code{type}).
#'
#' @returns The \code{phylo} tree object, returned invisibly.
#'
#' @author Erik Ringen \email{erikjacob.ringen@@uzh.ch}
#'
#' @details This function is a convenience wrapper around
#'   \code{\link[ape]{plot.phylo}} and \code{\link[ape]{nodelabels}} that
#'   displays internal node IDs using ape's numbering convention: tips are
#'   numbered \code{1:N_tips} and internal nodes are
#'   \code{(N_tips + 1):(N_tips + N_internal)}.
#'
#'   Node label size scales automatically with tree size when
#'   \code{node_cex = NULL}. For large trees, consider hiding tip labels
#'   (the default) and using the interactive workflow below to identify
#'   nodes of interest.
#'
#'   In an interactive R session, the returned tree can be passed to
#'   \code{identify()} to click on nodes and retrieve their IDs:
#'
#'   \preformatted{
#'   tree <- coev_plot_node_labels(fit)
#'   selected <- identify(tree)
#'   coev_ancestral_states(fit, nodes = selected$nodes)
#'   }
#'
#' @seealso \code{\link{coev_ancestral_states}},
#'   \code{\link[ape]{plot.phylo}}, \code{\link[ape]{identify.phylo}}
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
#'   chains = 4,
#'   parallel_chains = 4,
#'   seed = 1
#' )
#'
#' # plot tree with node labels (no tip labels by default)
#' coev_plot_node_labels(fit)
#'
#' # show tip labels and customise node appearance
#' coev_plot_node_labels(fit, tip_labels = TRUE, node_cex = 0.4,
#'                       cex = 0.4)
#'
#' # interactive node selection (in interactive R sessions)
#' tree <- coev_plot_node_labels(fit)
#' selected <- identify(tree)
#' coev_ancestral_states(fit, nodes = selected$nodes)
#' }
#'
#' @export
coev_plot_node_labels <- function(object, tip_labels = FALSE,
                                  node_cex = NULL, node_col = "#B2182B",
                                  ...) {
  if (!methods::is(object, "coevfit")) {
    stop2(
      paste0(
        "Argument 'object' must be a fitted coevolutionary model ",
        "of class coevfit."
      )
    )
  }
  if (!is.logical(tip_labels) || length(tip_labels) != 1) {
    stop2("Argument 'tip_labels' must be TRUE or FALSE.")
  }
  if (!is.null(node_cex)) {
    if (!is.numeric(node_cex) || length(node_cex) != 1 || node_cex <= 0) {
      stop2("Argument 'node_cex' must be a single positive number.")
    }
  }
  tree <- object$tree
  if (inherits(tree, "multiPhylo")) tree <- tree[[1]]
  n_tips <- length(tree$tip.label)
  internal_ids <- (n_tips + 1):(n_tips + tree$Nnode)
  # auto-scale node label size based on number of tips
  if (is.null(node_cex)) {
    node_cex <- min(0.8, max(0.25, 5 / sqrt(n_tips)))
  }
  ape::plot.phylo(
    tree,
    show.tip.label = tip_labels,
    edge.width = 1.5,
    edge.color = "grey50",
    no.margin = !tip_labels,
    ...
  )
  ape::nodelabels(
    text = internal_ids,
    cex = node_cex,
    frame = "none",
    col = node_col,
    font = 2,
    adj = c(-0.3, 0.5)
  )
  invisible(tree)
}
