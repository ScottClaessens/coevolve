#' Interactively select multiple nodes from a plotted phylogeny
#'
#' A convenience wrapper around \code{\link[ape]{identify.phylo}} that
#' allows the user to click on \emph{multiple} internal nodes in succession
#' and returns the collected node IDs. By default,
#' \code{\link[ape]{identify.phylo}} returns after a single click; this
#' function calls it in a loop until the user signals completion (by
#' clicking outside the tree, pressing Esc, or right-clicking, depending
#' on graphics device) or until \code{n} clicks have been collected.
#'
#' Intended workflow:
#'
#' \preformatted{
#' tree <- coev_plot_node_labels(fit)
#' nodes <- coev_identify_nodes(tree)
#' coev_ancestral_states(fit, nodes = nodes)
#' }
#'
#' @param tree A \code{phylo} object, typically the value returned (invisibly)
#'   by \code{\link{coev_plot_node_labels}}.
#' @param n (optional) Maximum number of nodes to collect. If \code{NULL}
#'   (default), the loop continues until the user cancels.
#'
#' @returns An integer vector of unique node IDs, in the order they were
#'   clicked.
#'
#' @author Erik Ringen \email{erikjacob.ringen@@uzh.ch}
#'
#' @seealso \code{\link{coev_plot_node_labels}},
#'   \code{\link{coev_ancestral_states}},
#'   \code{\link[ape]{identify.phylo}}
#'
#' @examples
#' \dontrun{
#' # plot tree with internal node IDs, then click on nodes of interest
#' tree <- coev_plot_node_labels(fit)
#' nodes <- coev_identify_nodes(tree)
#' coev_ancestral_states(fit, nodes = nodes)
#' }
#'
#' @export
coev_identify_nodes <- function(tree, n = NULL) {
  if (!inherits(tree, "phylo")) {
    stop2("Argument 'tree' must be a 'phylo' object.")
  }
  if (!is.null(n)) {
    if (!is.numeric(n) || length(n) != 1 || n < 1 || as.integer(n) != n) {
      stop2("Argument 'n' must be a single positive integer.")
    }
  }
  if (!interactive()) {
    stop2("coev_identify_nodes() requires an interactive R session.") # nocov
  }
  # nocov start: interactive code path; cannot be exercised in non-interactive
  # test runs, but the loop logic is straightforward.
  collected <- integer()
  i <- 0L
  repeat {
    i <- i + 1L
    if (!is.null(n) && i > n) break
    sel <- tryCatch(
      ape::identify.phylo(tree, quiet = TRUE),
      error = function(e) NULL
    )
    if (is.null(sel) || length(sel$nodes) == 0L) break
    collected <- c(collected, sel$nodes)
  }
  unique(collected)
  # nocov end
}
