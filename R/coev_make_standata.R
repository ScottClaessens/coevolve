#' Make Stan data for dynamic coevolutionary model
#'
#' @param data An object of class \code{data.frame} (or one that can be coerced
#'   to that class) containing data of all variables used in the model.
#' @param variables A named list identifying variables that should coevolve in
#'   the model and their associated response distributions as character strings (e.g.
#'   \code{list(var1 = "bernoulli_logit", var2 = "ordered_logistic")}). Must identify
#'   at least two variables. Variable names must refer to valid column names in data.
#'   Currently, the only supported response distributions are \code{bernoulli_logit},
#'   \code{ordered_logistic}, \code{poisson_softmax}, \code{normal}, and \code{lognormal}.
#' @param id A character of length one identifying the variable in the data that
#'   links rows to tips on the phylogeny. Must refer to a valid column name in
#'   the data. The id column must exactly match the tip labels in the phylogeny.
#' @param tree A phylogenetic tree object of class \code{phylo}.
#' @param dist_mat (optional) A distance matrix with row and column names exactly
#'   matching the tip labels in the phylogeny. If specified, the model will
#'   additionally control for spatial location by including a separate Gaussian
#'   Process over locations for every coevolving variable in the model.
#' @param prior (optional) A named list of priors for the model. If not specified,
#'   the model uses default priors (see Stan code). Alternatively, the user can
#'   specify a named list of priors. The list must contain non-duplicated entries
#'   for any of the following variables: the autoregressive and cross-effects
#'   (\code{alpha}), the drift scale parameters (\code{sigma}), the continuous
#'   time intercepts (\code{b}), the ancestral states for the traits (\code{eta_anc}),
#'   the cutpoints for ordinal variables (\code{c}), the sigma parameter(s) for
#'   Gaussian Processes over locations (\code{sigma_dist}), and the rho parameter(s)
#'   for Gaussian Processes over locations (\code{rho_dist}). These must be
#'   entered with valid prior strings, e.g. \code{list(alpha = "normal(0, 2)")}.
#'   Invalid prior strings will throw an error when the function internally checks
#'   the syntax of resulting Stan code.
#' @param prior_only Logical. If \code{FALSE} (default), the model is fitted to
#'   the data and returns a posterior distribution. If \code{TRUE}, the model
#'   samples from the prior only, ignoring the likelihood.
#'
#' @return A list containing the data for fitting the dynamic coevolutionary model in \pkg{Stan}.
#' @export
#'
#' @examples
#' # simulate data
#' n <- 20
#' tree <- ape::rcoal(n)
#' d <- data.frame(
#'    id = tree$tip.label,
#'    x = rbinom(n, size = 1, prob = 0.5),
#'    y = ordered(sample(1:4, size = n, replace = TRUE))
#' )
#' # make stan code
#' coev_make_standata(
#'    data = d,
#'    variables = list(
#'        x = "bernoulli_logit",
#'        y = "ordered_logistic"
#'    ),
#'    id = "id",
#'    tree = tree
#' )
coev_make_standata <- function(data, variables, id, tree, dist_mat = NULL,
                               prior = NULL, prior_only = FALSE) {
  # check arguments
  run_checks(data, variables, id, tree, dist_mat, prior, prior_only)
  # match data to tree tip label ordering
  data <- data[match(tree$tip.label, data[,id]),]
  # stop if data and tips do not match
  if (!identical(tree$tip.label, data[,id])) {
    stop2("Data and phylogeny tips do not match.")
  }
  # cut up tree into segments
  times <- ape::node.depth.edgelength(tree)
  # line up date of each node with the split points in the tree
  split_points <- sort(unique(times))
  node_time <- match(times, split_points)
  # create a sequence of nodes, respecting temporal order
  node_seq <- seq(from = 1, to = length(node_time))
  node_seq <- node_seq[order(node_time)]
  # find the "parent" node for each node and tip of the tree
  parent <- phangorn::Ancestors(tree, node_seq, type = "parent")
  # parent time indicates amount of time since the parent node
  # scaled by the total depth of the tree
  parent_time <- rep(NA, length(node_seq))
  parent_time[1] <- -99 # placeholder for ancestral state
  for (i in 2:length(parent_time)) {
    parent_time[i] <-
      (ape::node.depth.edgelength(tree)[node_seq[i]] - ape::node.depth.edgelength(tree)[parent[i]]) /
      max(ape::node.depth.edgelength(tree))
  }
  # get total num segments in the tree
  N_seg <- length(node_seq)
  # indicate whether a node in the seq is a tip
  tip <- ifelse(node_seq > length(tree$tip.label), 0, 1)
  # normalise distance matrix so that maximum distance = 1
  if (!is.null(dist_mat)) dist_mat <- dist_mat / max(dist_mat)
  # data list for stan
  sd <- list(
    N = length(tree$tip.label), # number of taxa
    J = length(variables),      # number of variables
    N_seg = N_seg,              # number of segments in the tree
    node_seq = node_seq,        # sequence of nodes
    parent = parent,            # parent node for each node
    ts = parent_time,           # amount of time since parent node
    tip = tip                   # is tip?
  )
  # add observed data variables one-by-one
  for (i in 1:length(variables)) {
    if (variables[[i]] %in% c("normal", "lognormal")) {
      sd[[paste0("y", i)]] <- as.numeric(data[,names(variables)[i]])
    } else if (variables[[i]] %in% c("bernoulli_logit", "ordered_logistic", "poisson_softmax")) {
      sd[[paste0("y", i)]] <- as.integer(data[,names(variables)[i]])
    }
  }
  # add distance matrix if specified
  if (!is.null(dist_mat)) sd[["dist_mat"]] <- dist_mat
  # add prior_only
  sd[["prior_only"]] <- as.numeric(prior_only)
  return(sd)
}
