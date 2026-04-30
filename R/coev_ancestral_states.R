#' Extract ancestral state estimates from a fitted coevolutionary model
#'
#' Extracts posterior estimates of trait values at internal nodes of the
#' phylogenetic tree from a fitted \code{coevfit} object. By default,
#' estimates are returned on the latent (eta) scale; use
#' \code{scale = "response"} to apply the inverse link for each variable's
#' response distribution. The latent \code{eta} parameters in the model
#' correspond directly to nodes in the phylogeny, making ancestral state
#' reconstruction a natural byproduct of model fitting.
#'
#' @param object An object of class \code{coevfit}
#' @param variables (optional) A character vector of variable names to include.
#'   If \code{NULL} (default), all variables from the fitted model are included.
#' @param nodes Which nodes to include. Either \code{"internal"} (default) for
#'   internal nodes only, \code{"all"} for both tips and internal nodes, or an
#'   integer vector of specific node IDs (using ape's node
#'   numbering: tips = 1:N_tips, internal = (N_tips+1):total).
#' @param tree_id (optional) An integer indicating which tree to use when the
#'   model was fit with a \code{multiPhylo} object. If \code{NULL} (default),
#'   uses tree 1. (Integration over trees with topological uncertainty is
#'   not yet implemented.)
#' @param scale Character string, either \code{"latent"} (default) to return
#'   values on the latent eta scale, or \code{"response"} to apply the inverse
#'   link function for each variable's response distribution.
#' @param summary Logical. If \code{TRUE} (default), returns a long-format
#'   data frame with posterior point estimates and credible intervals. If
#'   \code{FALSE}, returns the raw posterior draws.
#' @param prob A single numeric value between 0 and 1 indicating the desired
#'   probability mass for the credible interval. Default is 0.95.
#'
#' @returns If \code{summary = TRUE}, a long-format data frame (tibble) with
#'   columns:
#'   \describe{
#'     \item{node}{Integer node ID in the reference tree}
#'     \item{variable}{Character variable name}
#'     \item{category}{For ordinal variables on the response scale, the
#'       category label (\code{"cat_1"}, \code{"cat_2"}, ...). \code{NA}
#'       otherwise.}
#'     \item{estimate}{Posterior point estimate. For ordinal-response
#'       category probabilities this is the posterior \emph{mean} (so the
#'       per-node category estimates sum to 1, matching the convention of
#'       \code{predict.brmsfit}); for all other cases it is the posterior
#'       median.}
#'     \item{lower}{Lower bound of the credible interval (per-category
#'       quantile for ordinal-response rows)}
#'     \item{upper}{Upper bound of the credible interval (per-category
#'       quantile for ordinal-response rows)}
#'   }
#'   The data frame has attributes \code{ref_tree} (the phylo object used),
#'   \code{prob}, and \code{scale}.
#'
#'   If \code{summary = FALSE}, a list with elements:
#'   \describe{
#'     \item{draws}{A 3D array (draws x nodes x variables)}
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
coev_ancestral_states <- function(object, variables = NULL,
                                  nodes = "internal",
                                  tree_id = NULL,
                                  scale = "latent",
                                  summary = TRUE,
                                  prob = 0.95) {
  # input validation
  run_checks_ancestral_states(
    object, variables, nodes, tree_id, scale, summary, prob
  )
  # get tree(s) from fitted model
  tree <- object$tree
  if (!inherits(tree, "multiPhylo")) {
    tree <- list(tree)
    class(tree) <- "multiPhylo"
  }
  n_tips <- object$stan_data$N_tips
  n_seg <- object$stan_data$N_seg
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
  # apply response-scale transform if requested
  if (scale == "response") {
    # compute mean of observed (non-missing) values per variable
    # needed for poisson_softplus and negative_binomial_softplus
    obs_means <- compute_obs_means(object, var_idx)
    # get cutpoints for ordered_logistic variables
    cutpoints <- extract_cutpoints(post, distributions, var_idx)
    draws_subset <- apply_inverse_link(
      draws_subset, distributions[var_idx], obs_means, cutpoints
    )
  }
  if (summary) {
    # compute posterior summaries (long format: one row per
    # node x variable x category)
    probs <- c((1 - prob) / 2, 1 - (1 - prob) / 2)
    rows <- list()
    for (n_i in seq_along(selected_nodes)) {
      for (v_i in seq_along(var_idx)) {
        dist <- distributions[var_idx[v_i]]
        if (scale == "response" && dist == "ordered_logistic") {
          cat_draws <- draws_subset[[v_i]][, n_i, , drop = FALSE]
          n_cats <- dim(cat_draws)[3]
          for (k in seq_len(n_cats)) {
            cat_k <- cat_draws[, 1, k]
            # use mean so per-node category estimates sum to 1
            rows[[length(rows) + 1]] <- data.frame(
              node = selected_nodes[n_i],
              variable = var_names[v_i],
              category = paste0("cat_", k),
              estimate = mean(cat_k),
              lower = stats::quantile(cat_k, probs[1], names = FALSE),
              upper = stats::quantile(cat_k, probs[2], names = FALSE),
              stringsAsFactors = FALSE
            )
          }
        } else {
          if (is.list(draws_subset)) {
            node_draws <- draws_subset[[v_i]][, n_i, 1]
          } else {
            node_draws <- draws_subset[, n_i, v_i]
          }
          rows[[length(rows) + 1]] <- data.frame(
            node = selected_nodes[n_i],
            variable = var_names[v_i],
            category = NA_character_,
            estimate = stats::median(node_draws),
            lower = stats::quantile(node_draws, probs[1], names = FALSE),
            upper = stats::quantile(node_draws, probs[2], names = FALSE),
            stringsAsFactors = FALSE
          )
        }
      }
    }
    result <- do.call(dplyr::bind_rows, rows)
    result <- tibble::as_tibble(result)
    attr(result, "ref_tree") <- ref_tree
    attr(result, "prob") <- prob
    attr(result, "scale") <- scale
    result
  } else {
    # return raw draws
    if (!is.list(draws_subset)) {
      dimnames(draws_subset) <- list(
        draw = NULL,
        node = selected_nodes,
        variable = var_names
      )
    }
    list(
      draws = draws_subset,
      ref_tree = ref_tree,
      node_ids = selected_nodes,
      variable_names = var_names
    )
  }
}

#' Compute mean of observed (non-missing) values per variable
#' @noRd
compute_obs_means <- function(object, var_idx) {
  var_names <- names(object$variables)
  data <- object$data
  obs_means <- numeric(length(var_idx))
  for (i in seq_along(var_idx)) {
    vname <- var_names[var_idx[i]]
    vals <- data[[vname]]
    # remove NAs and convert to numeric
    vals <- as.numeric(vals[!is.na(vals)])
    obs_means[i] <- mean(vals)
  }
  obs_means
}

#' Extract cutpoints for ordered_logistic variables from posterior
#' @noRd
extract_cutpoints <- function(post, distributions, var_idx) {
  cutpoints <- list()
  all_dists <- distributions
  for (i in seq_along(var_idx)) {
    j <- var_idx[i]
    if (all_dists[j] == "ordered_logistic") {
      # cutpoints are named c1, c2, ... in posterior (by original position)
      cp_name <- paste0("c", j)
      if (cp_name %in% names(post)) {
        cutpoints[[i]] <- post[[cp_name]] # [draws, n_cutpoints]
      }
    } else {
      cutpoints[[i]] <- NULL
    }
  }
  cutpoints
}

#' Apply inverse link function to draws
#' @noRd
apply_inverse_link <- function(draws, distributions, obs_means, cutpoints) {
  n_draws <- dim(draws)[1]
  n_nodes <- dim(draws)[2]
  n_vars <- dim(draws)[3]
  has_ordinal <- any(distributions == "ordered_logistic")
  if (has_ordinal) {
    # need to return a list since ordinal variables have different dimensions
    result <- vector("list", n_vars)
    for (v in seq_len(n_vars)) {
      dist <- distributions[v]
      eta_v <- draws[, , v] # [draws, nodes]
      if (is.null(dim(eta_v))) {
        eta_v <- matrix(eta_v, nrow = n_draws, ncol = n_nodes)
      }
      if (dist == "ordered_logistic" && !is.null(cutpoints[[v]])) {
        cp <- cutpoints[[v]] # [draws, n_cutpoints]
        if (is.null(dim(cp))) cp <- matrix(cp, ncol = 1) # nocov
        n_cats <- ncol(cp) + 1
        # result: [draws, nodes, categories]
        cat_probs <- array(NA, dim = c(n_draws, n_nodes, n_cats))
        for (d in seq_len(n_draws)) {
          for (n in seq_len(n_nodes)) {
            cat_probs[d, n, ] <- ordered_logistic_probs(
              eta_v[d, n], cp[d, ]
            )
          }
        }
        result[[v]] <- cat_probs
      } else {
        # wrap scalar transform in 3D array [draws, nodes, 1]
        transformed <- transform_eta(eta_v, dist, obs_means[v])
        result[[v]] <- array(transformed, dim = c(n_draws, n_nodes, 1))
      }
    }
    result
  } else {
    # all variables are non-ordinal: apply element-wise transforms
    for (v in seq_len(n_vars)) {
      eta_v <- draws[, , v]
      if (is.null(dim(eta_v))) {
        eta_v <- matrix(eta_v, nrow = n_draws, ncol = n_nodes)
      }
      draws[, , v] <- transform_eta(eta_v, distributions[v], obs_means[v])
    }
    draws
  }
}

#' Apply scalar inverse link for a single distribution
#' @noRd
transform_eta <- function(eta, distribution, obs_mean = NULL) {
  switch(
    distribution,
    "normal" = eta,
    "bernoulli_logit" = stats::plogis(eta),
    "poisson_softplus" = obs_mean * log1p(exp(eta)),
    "negative_binomial_softplus" = obs_mean * log1p(exp(eta)),
    "gamma_log" = exp(eta),
    eta # nocov
  )
}

#' Compute ordered logistic probabilities from eta and cutpoints
#' @noRd
ordered_logistic_probs <- function(eta, cutpoints) {
  n_cats <- length(cutpoints) + 1
  probs <- numeric(n_cats)
  for (k in seq_len(n_cats)) {
    if (k == 1) {
      probs[k] <- stats::plogis(cutpoints[1] - eta)
    } else if (k == n_cats) {
      probs[k] <- 1 - stats::plogis(cutpoints[k - 1] - eta)
    } else {
      probs[k] <- stats::plogis(cutpoints[k] - eta) -
        stats::plogis(cutpoints[k - 1] - eta)
    }
  }
  probs
}

#' Input validation for coev_ancestral_states
#' @noRd
run_checks_ancestral_states <- function(object, variables, nodes,
                                        tree_id, scale,
                                        summary, prob) {
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
