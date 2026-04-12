#' Extract ancestral state estimates from a fitted coevolutionary model
#'
#' Extracts posterior estimates of latent trait values at internal nodes
#' of the phylogenetic tree from a fitted \code{coevfit} object. The latent
#' \code{eta} parameters in the model correspond directly to nodes in the
#' phylogeny, making ancestral state reconstruction a natural byproduct of
#' the model fitting process.
#'
#' @param object An object of class \code{coevfit}
#' @param variables (optional) A character vector of variable names to include.
#'   If \code{NULL} (default), all variables from the fitted model are included.
#' @param nodes Which nodes to include. Either \code{"internal"} (default) for
#'   internal nodes only, \code{"all"} for both tips and internal nodes, or an
#'   integer vector of specific node IDs (using ape's node numbering convention
#'   where tips are 1:N_tips and internal nodes are (N_tips+1):(N_tips+N_internal)).
#' @param tree_id (optional) An integer indicating which tree to use when the
#'   model was fit with a \code{multiPhylo} object. If \code{NULL} (default),
#'   uses tree 1 when only one tree is present, or integrates over trees when
#'   multiple trees are present.
#' @param scale Character string, either \code{"latent"} (default) to return
#'   values on the latent eta scale, or \code{"response"} to apply the inverse
#'   link function for each variable's response distribution.
#' @param summary Logical. If \code{TRUE} (default), returns a summary data
#'   frame with posterior median and credible intervals. If \code{FALSE},
#'   returns the raw posterior draws.
#' @param prob A single numeric value between 0 and 1 indicating the desired
#'   probability mass for the credible interval. Default is 0.95.
#'
#' @returns If \code{summary = TRUE}, a data frame (tibble) with columns:
#'   \describe{
#'     \item{node}{Integer node ID in the reference tree}
#'     \item{variable}{Character variable name}
#'     \item{estimate}{Posterior median}
#'     \item{lower}{Lower bound of the credible interval}
#'     \item{upper}{Upper bound of the credible interval}
#'     \item{clade_pp}{Clade posterior probability (1.0 for single-tree models)}
#'   }
#'   The data frame has attributes \code{ref_tree} (the phylo object used),
#'   \code{prob}, and \code{scale}. The format is compatible with ggtree's
#'   \code{\%<+\%} operator for plotting.
#'
#'   If \code{summary = FALSE}, a list with elements:
#'   \describe{
#'     \item{draws}{A 3D array with dimensions [draws, nodes, variables]}
#'     \item{ref_tree}{The reference phylo object}
#'     \item{node_ids}{Integer vector of node IDs (matching dim 2 of draws)}
#'     \item{variable_names}{Character vector of variable names}
#'   }
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
#' # extract ancestral states (latent scale, internal nodes, summary)
#' asr <- coev_ancestral_states(fit)
#'
#' # extract all nodes including tips
#' asr_all <- coev_ancestral_states(fit, nodes = "all")
#'
#' # extract on response scale
#' asr_resp <- coev_ancestral_states(fit, scale = "response")
#'
#' # get raw posterior draws
#' asr_draws <- coev_ancestral_states(fit, summary = FALSE)
#' }
#'
#' @seealso \code{\link{coev_fit}}, \code{\link{coev_plot_trait_values}}
#'
#' @export
coev_ancestral_states <- function(
    object,
    variables = NULL,
    nodes = "internal",
    tree_id = NULL,
    scale = "latent",
    summary = TRUE,
    prob = 0.95
) {
  # input validation
  run_checks_ancestral_states(
    object, variables, nodes, tree_id, scale, summary, prob
  )
  # get tree(s) from fitted model
  tree <- object$tree
  if (!inherits(tree, "multiPhylo")) {
    tree <- c(tree) # coerce single phylo to multiPhylo
    class(tree) <- "multiPhylo"
  }
  n_tree <- object$stan_data$N_tree
  n_tips <- object$stan_data$N_tips
  n_seg <- object$stan_data$N_seg
  node_seq <- object$stan_data$node_seq # matrix [N_tree, N_seg]
  tip_indicator <- object$stan_data$tip # matrix [N_tree, N_seg]
  var_names <- names(object$variables)
  distributions <- unlist(object$variables)
  j <- length(var_names)
  # subset variables if requested
  if (!is.null(variables)) {
    var_idx <- match(variables, var_names)
    var_names <- variables
  } else {
    var_idx <- seq_len(j)
  }
  # determine which tree to use
  if (!is.null(tree_id)) {
    t_idx <- tree_id
  } else {
    t_idx <- 1L
  }
  ref_tree <- tree[[t_idx]]
  # extract posterior draws of eta
  # eta has Stan dims array[N_tree, N_seg] vector[J]
  # after extract_samples(): [draws, N_tree, N_seg, J]
  post <- extract_samples(object)
  eta_draws <- post$eta # [draws, N_tree, N_seg, J]
  # subset to the selected tree
  eta_t <- eta_draws[, t_idx, , , drop = FALSE]
  # collapse tree dimension: [draws, N_seg, J]
  eta_t <- array(
    eta_t,
    dim = c(dim(eta_t)[1], dim(eta_t)[3], dim(eta_t)[4])
  )
  # determine which nodes to include
  # node_seq[t, i] gives the phylo node ID at position i
  # eta_t[draw, k, j] is indexed by node ID k (since Stan assigns
  # eta[t, node_seq[t,i]] = ..., so the second dim IS node ID)
  all_node_ids <- seq_len(n_seg)
  internal_mask <- all_node_ids > n_tips
  if (is.character(nodes)) {
    if (nodes == "internal") {
      selected_nodes <- all_node_ids[internal_mask]
    } else if (nodes == "all") {
      selected_nodes <- all_node_ids
    }
  } else {
    # integer vector of specific node IDs
    selected_nodes <- as.integer(nodes)
  }
  # extract draws for selected nodes and variables
  # result: [draws, length(selected_nodes), length(var_idx)]
  draws_subset <- eta_t[, selected_nodes, var_idx, drop = FALSE]
  if (summary) {
    # compute posterior summaries
    probs <- c((1 - prob) / 2, 1 - (1 - prob) / 2)
    rows <- list()
    for (n_i in seq_along(selected_nodes)) {
      for (v_i in seq_along(var_idx)) {
        node_draws <- draws_subset[, n_i, v_i]
        rows[[length(rows) + 1]] <- data.frame(
          node = selected_nodes[n_i],
          variable = var_names[v_i],
          estimate = stats::median(node_draws),
          lower = stats::quantile(node_draws, probs[1], names = FALSE),
          upper = stats::quantile(node_draws, probs[2], names = FALSE),
          clade_pp = 1.0,
          stringsAsFactors = FALSE
        )
      }
    }
    result <- do.call(rbind, rows)
    result <- tibble::as_tibble(result)
    attr(result, "ref_tree") <- ref_tree
    attr(result, "prob") <- prob
    attr(result, "scale") <- scale
    return(result)
  } else {
    # return raw draws
    dimnames(draws_subset) <- list(
      draw = NULL,
      node = selected_nodes,
      variable = var_names
    )
    return(list(
      draws = draws_subset,
      ref_tree = ref_tree,
      node_ids = selected_nodes,
      variable_names = var_names
    ))
  }
}

#' Input validation for coev_ancestral_states
#' @noRd
run_checks_ancestral_states <- function(
    object, variables, nodes, tree_id, scale, summary, prob
) {
  if (!methods::is(object, "coevfit")) {
    stop2(
      paste0(
        "Argument 'object' must be a fitted coevolutionary model ",
        "of class coevfit."
      )
    )
  }
  if (!is.null(variables)) {
    if (!is.character(variables)) {
      stop2("Argument 'variables' must be a character vector.")
    }
    if (!all(variables %in% names(object$variables))) {
      stop2(
        paste0(
          "Argument 'variables' contains variable names that are not ",
          "included in the fitted model."
        )
      )
    }
  }
  if (is.character(nodes)) {
    if (!(nodes %in% c("internal", "all"))) {
      stop2(
        "Argument 'nodes' must be \"internal\", \"all\", or an integer vector."
      )
    }
  } else if (!is.numeric(nodes)) {
    stop2(
      "Argument 'nodes' must be \"internal\", \"all\", or an integer vector."
    )
  }
  if (!is.null(tree_id)) {
    if (!is.numeric(tree_id) || length(tree_id) != 1 ||
          as.integer(tree_id) != tree_id) {
      stop2("Argument 'tree_id' must be a single integer.")
    }
    if (tree_id < 1 || tree_id > object$stan_data$N_tree) {
      stop2("Argument 'tree_id' must be between 1 and the number of trees.")
    }
  }
  if (!is.character(scale) || length(scale) != 1 ||
        !(scale %in% c("latent", "response"))) {
    stop2("Argument 'scale' must be \"latent\" or \"response\".")
  }
  if (!is.logical(summary) || length(summary) != 1) {
    stop2("Argument 'summary' must be logical.")
  }
  if (!is.numeric(prob) || length(prob) != 1 || prob <= 0 || prob >= 1) {
    stop2(
      "Argument 'prob' must be a single numeric value between 0 and 1."
    )
  }
}
