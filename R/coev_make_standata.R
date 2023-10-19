#' Make Stan data for dynamic coevolutionary model
#'
#' @param data An object of class \code{data.frame} (or one that can be coerced
#'   to that class) containing data of all variables used in the model.
#' @param variables A named list identifying variables that should coevolve in
#'   the model and their associated response distributions as character strings (e.g.
#'   \code{list(var1 = "bernoulli_logit", var2 = "ordered_logistic")}). Must identify
#'   at least two variables. Variable names must refer to valid column names in data.
#'   Currently, the only supported response distributions are \code{bernoulli_logit}
#'   and \code{ordered_logistic}.
#' @param id A character of length one identifying the variable in the data that links rows to tips
#'   on the phylogeny. Must refer to a valid column name in the data. The id column
#'   must exactly match the tip labels in the phylogeny.
#' @param tree A phylogenetic tree object of class \code{phylo}.
#' @param prior A list of priors for the model.
#'
#' @return A list containing the data for fitting the dynamic coevolutionary model in \pkg{Stan}.
#' @export
#'
#' @examples
#' # simulate data
#' set.seed(1)
#' n <- 20
#' tree <- ape::rtree(n)
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
coev_make_standata <- function(data, variables, id, tree, prior = NULL) {
  # check arguments
  run_checks(data, variables, id, tree)
  # match data to tree tip label ordering
  data <- data[match(tree$tip.label, data[,id]),]
  # stop if data and tips do not match
  if (!identical(tree$tip.label, data[,id])) {
    stop("Data and phylogeny tips do not match.")
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
  # get observed data
  obs <- list()
  for (j in 1:length(variables)) {
    obs[[names(variables)[j]]] <- as.integer(data[,names(variables)[j]])
  }
  obs <- as.matrix(as.data.frame(obs))
  # data list for stan
  sd <- list(
    N = length(tree$tip.label), # number of taxa
    J = length(variables),      # number of variables
    N_seg = N_seg,              # number of segments in the tree
    node_seq = node_seq,        # sequence of nodes
    parent = parent,            # parent node for each node
    ts = parent_time,           # amount of time since parent node
    y = obs                     # observed data
  )
  return(sd)
}
