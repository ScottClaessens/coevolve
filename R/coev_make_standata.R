#' Make Stan data for dynamic coevolutionary model
#'
#' Make the data list for the \pkg{Stan} model. The function takes a dataset,
#' phylogeny, and a set of variables and performs the necessary computations
#' (e.g., cutting up the tree into segments, computing branch lengths,
#' determining parent and child nodes) for the \pkg{Stan} model.
#'
#' @param data An object of class \code{data.frame} (or one that can be coerced
#'   to that class) containing data of all variables used in the model.
#' @param variables A named list identifying variables that should coevolve in
#'   the model and their associated response distributions as character strings
#'   (e.g. \code{list(var1 = "bernoulli_logit", var2 = "ordered_logistic")}).
#'   Must identify at least two variables. Variable names must refer to valid
#'   column names in data. Currently, the only supported response distributions
#'   are \code{bernoulli_logit}, \code{ordered_logistic},
#'   \code{poisson_softplus}, \code{negative_binomial_softplus}, \code{normal},
#'   and \code{gamma_log}.
#' @param id A character of length one identifying the variable in the data that
#'   links rows to tips on the phylogeny. Must refer to a valid column name in
#'   the data. The id column must exactly match the tip labels in the phylogeny.
#' @param tree A phylogenetic tree object of class \code{phylo} or
#'   \code{multiPhylo}. The tree(s) must be rooted and must include positive
#'   non-zero branch lengths. All trees in \code{multiPhylo} objects must have
#'   the same number of internal nodes and branches.
#' @param effects_mat (optional) A boolean matrix with row and column names
#'   exactly matching the variables declared for the model. If not specified,
#'   all cross-lagged effects will be estimated in the model. If specified, the
#'   model will only estimate cross-lagged effects where cells in the matrix are
#'   TRUE and will ignore cross-lagged effects where cells in the matrix are
#'   FALSE. In the matrix, columns represent predictor variables and rows
#'   represent outcome variables. All autoregressive effects (e.g., X -> X) must
#'   be TRUE in the matrix.
#' @param complete_cases (optional) Logical. If \code{FALSE} (default), all
#'   missing values are imputed by the model. If \code{TRUE}, taxa with missing
#'   data are excluded.
#' @param dist_mat (optional) A distance matrix with row and column names
#'   exactly matching the tip labels in the phylogeny. If specified, the model
#'   will additionally control for spatial location by including a separate
#'   Gaussian Process over locations for every coevolving variable in the model.
#' @param dist_cov A string specifying the covariance kernel used for Gaussian
#'   Processes over locations. Currently supported are \code{"exp_quad"}
#'   (exponentiated-quadratic kernel; default), \code{"exponential"}
#'   (exponential kernel), and \code{"matern32"} (Matern 3/2 kernel).
#' @param measurement_error (optional) A named list of coevolving variables and
#'   their associated columns in the dataset containing standard errors. Only
#'   valid for normally-distributed variables. For example, if we declare
#'   \code{variables = list(x = "normal", y = "normal")}, then we could set
#'   \code{measurement_error = list(x = "x_std_err")} to tell the function to
#'   include measurement error on \code{x} using standard errors from the
#'   \code{x_std_err} column of the dataset.
#' @param prior (optional) A named list of priors for the model. If not
#'   specified, the model uses default priors (see \code{help(coev_fit)}).
#'   Alternatively, the user can specify a named list of priors. The list must
#'   contain non-duplicated entries for any of the following parameters: the
#'   autoregressive effects (\code{A_diag}), the cross effects
#'   (\code{A_offdiag}), the Cholesky factor for the drift matrix (\code{L_R}),
#'   the drift std. dev. parameters (\code{Q_sigma}), the continuous time
#'   intercepts (\code{b}), the ancestral states for the traits
#'   (\code{eta_anc}), the cutpoints for ordinal variables (\code{c}), the
#'   overdispersion parameters for negative binomial variables (\code{phi}),
#'   the shape parameters for gamma variables (\code{shape}), the sigma
#'   parameters for Gaussian Processes over locations (\code{sigma_dist}), the
#'   rho parameters for Gaussian Processes over locations (\code{rho_dist}), the
#'   standard deviation parameters for non-phylogenetic group-level varying
#'   effects (\code{sigma_group}), and the Cholesky factor for the
#'   non-phylogenetic group-level correlation matrix (\code{L_group}). These
#'   must be entered with valid prior strings, e.g.
#'   \code{list(A_offdiag = "normal(0, 2)")}.
#' @param scale Logical. If \code{TRUE} (default), variables following the
#'   \code{normal} and \code{gamma_log} response distributions are scaled before
#'   fitting the model. Continuous variables following the \code{normal}
#'   distribution are standardised (e.g., mean centered and divided by their
#'   standard deviation) and positive real variables following the
#'   \code{gamma_log} distribution are divided by the mean value without
#'   centering. This approach is recommended when using default priors to
#'   improve efficiency and ensure accurate inferences. If \code{FALSE},
#'   variables are left unscaled for model fitting. In this case, users should
#'   take care to set sensible priors on variables.
#' @param estimate_Q_offdiag Logical. If \code{TRUE} (default), the model
#'   estimates the off-diagonals for the \deqn{Q} drift matrix (i.e., correlated
#'   drift). If \code{FALSE}, the off-diagonals for the \deqn{Q} drift matrix
#'   are set to zero.
#' @param log_lik Logical. Set to \code{FALSE} by default. If \code{TRUE}, the
#'   model returns the pointwise log likelihood, which can be used to calculate
#'   WAIC and LOO.
#' @param prior_only Logical. If \code{FALSE} (default), the model is fitted to
#'   the data and returns a posterior distribution. If \code{TRUE}, the model
#'   samples from the prior only, ignoring the likelihood.
#'
#' @return A list containing the data for fitting the dynamic coevolutionary
#'   model in \pkg{Stan}.
#'
#' @author Scott Claessens \email{scott.claessens@@gmail.com}, Erik Ringen
#'   \email{erikjacob.ringen@@uzh.ch}
#'
#' @details For further details, see \code{help(coev_fit)}
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
#' @seealso \code{\link{coev_make_stancode}}, \code{\link{coev_fit}}
#'
#' @examples
#' # make stan data
#' coev_make_standata(
#'   data = authority$data,
#'   variables = list(
#'     political_authority = "ordered_logistic",
#'     religious_authority = "ordered_logistic"
#'   ),
#'   id = "language",
#'   tree = authority$phylogeny
#'   )
#'
#' @export
coev_make_standata <- function(data, variables, id, tree,
                               effects_mat = NULL, complete_cases = FALSE,
                               dist_mat = NULL, dist_cov = "exp_quad",
                               measurement_error = NULL,
                               prior = NULL, scale = TRUE,
                               estimate_Q_offdiag = TRUE,
                               log_lik = FALSE,
                               prior_only = FALSE) {
  # check arguments
  run_checks(data, variables, id, tree, effects_mat, complete_cases, dist_mat,
             dist_cov, measurement_error, prior, scale, estimate_Q_offdiag,
             log_lik, prior_only)
  # coerce data argument to data frame
  data <- as.data.frame(data)
  # warning if scale = FALSE
  if (!scale) {
    warning2(
      paste0(
        "When scale = FALSE, continuous variables are left unstandardised in ",
        "the Stan data list. Users should take care to set sensible priors ",
        "for model fitting, rather than use default priors."
      )
    )
  }
  # if complete_cases = TRUE, remove data rows with NAs
  if (complete_cases) {
    any_missing <- apply(data[,names(variables)], 1, function(x) any(is.na(x)))
    data <- data[!any_missing,]
  }
  # coerce tree object to multiPhylo
  tree <- phytools::as.multiPhylo(tree)
  # ensure that all trees have same tip labels
  tree <- ape::.compressTipLabel(tree)
  # prune tree to dataset
  tree <- ape::keep.tip.multiPhylo(tree, data[,id])
  # match data ordering to tree tip label ordering
  matched_data <- data.frame()
  for (tip in tree[[1]]$tip.label) {
    matched_data <- rbind(matched_data, data[data[,id] == tip,])
  }
  data <- matched_data
  # match distance matrix to tree tip label ordering
  if (!is.null(dist_mat)) {
    dist_mat <- dist_mat[tree[[1]]$tip.label, tree[[1]]$tip.label]
  }
  # create effects matrix if not specified
  if (is.null(effects_mat)) {
    effects_mat <-
      matrix(TRUE, nrow = length(variables), ncol = length(variables),
             dimnames = list(names(variables), names(variables)))
  }
  # match effects_mat to vector of variables
  effects_mat <- effects_mat[names(variables),names(variables)]
  # convert effects_mat to integer matrix (unary conversion +)
  effects_mat <- +effects_mat
  # stop for internal mismatches
  if (!identical(tree[[1]]$tip.label, unique(data[,id]))) {
    stop2("Data and phylogeny tips do not match.")
  } else if (!is.null(dist_mat) & (
    !identical(tree[[1]]$tip.label, rownames(dist_mat)) |
    !identical(tree[[1]]$tip.label, colnames(dist_mat))
    )) {
    stop2("Distance matrix and phylogeny tips do not match.")
  } else if (!identical(names(variables), rownames(effects_mat)) |
             !identical(names(variables), colnames(effects_mat))) {
    stop2("Effects matrix and variable names do not match.")
  }
  # get number of trees
  N_tree <- length(tree)
  # get number of segments in trees
  N_seg <- length(ape::node.depth(tree[[1]]))
  # initialise tree variables for stan
  stan_node_seq <- matrix(NA, N_tree, N_seg)
  stan_parent <- matrix(NA, N_tree, N_seg)
  stan_ts <- matrix(NA, N_tree, N_seg)
  stan_tip <- matrix(NA, N_tree, N_seg)
  # loop over phylogenetic trees
  for (t in 1:length(tree)) {
    # cut up tree into segments
    times <- ape::node.depth.edgelength(tree[[t]])
    # line up date of each node with the split points in the tree
    split_points <- sort(unique(times))
    node_time <- match(times, split_points)
    # create a sequence of nodes, respecting temporal order
    node_seq <- seq(from = 1, to = length(node_time))
    node_seq <- node_seq[order(node_time)]
    # find the "parent" node for each node and tip of the tree
    parent <- phangorn::Ancestors(tree[[t]], node_seq, type = "parent")
    # parent time indicates amount of time since the parent node
    # scaled by the total depth of the tree
    parent_time <- rep(NA, length(node_seq))
    parent_time[1] <- -99 # placeholder for ancestral state
    for (i in 2:length(parent_time)) {
      parent_time[i] <-
        (ape::node.depth.edgelength(tree[[t]])[node_seq[i]] -
           ape::node.depth.edgelength(tree[[t]])[parent[i]]) /
        max(ape::node.depth.edgelength(tree[[t]]))
    }
    # indicate whether a node in the seq is a tip
    tip <- ifelse(node_seq > length(tree[[t]]$tip.label), 0, 1)
    # save variables for stan
    stan_node_seq[t,] <- node_seq
    stan_parent[t,] <- parent
    stan_ts[t,] <- parent_time
    stan_tip[t,] <- tip
  }
  # get data matrix
  y <- list()
  for (j in 1:length(variables)) {
    if (scale & variables[[j]] == "normal") {
      # standardised continuous variables
      y[[names(variables)[j]]] <- as.numeric(scale(data[,names(variables)[j]]))
    } else if (scale & variables[[j]] == "gamma_log") {
      # positive reals scaled by sample mean
      y[[names(variables)[j]]] <-
        as.numeric(
          data[,names(variables)[j]] /
            max(data[,names(variables)[j]], na.rm = TRUE)
          )
    } else {
      # unscaled binary/ordered/count variables
      y[[names(variables)[j]]] <- as.numeric(data[,names(variables)[j]])
    }
  }
  y <- as.matrix(as.data.frame(y))
  # get missing matrix
  miss <- ifelse(is.na(y), 1, 0)
  # replace y with -9999 if missing
  y[miss == 1] <- -9999
  # get matrix with squared standard errors
  if (!is.null(measurement_error)) {
    se <- list()
    for (j in 1:length(variables)) {
      if (names(variables)[j] %in% names(measurement_error)) {
        # if standard errors declared for this variable
        # get standard error column
        se_column <- measurement_error[[names(variables)[j]]]
        # get vector of standard error values
        se_values <- data[[se_column]]
        if (scale) {
          # if scale = TRUE, need to also scale the standard errors correctly
          se_values <- se_values / sd(data[[names(variables)[j]]], na.rm = TRUE)
        }
        # add squared standard errors to matrix, set any NAs to zero
        se[[names(variables)[j]]] <- ifelse(is.na(se_values), 0, se_values^2)
      } else {
        # if no standard errors declared, zero for entire column in matrix
        se[[names(variables)[j]]] <- rep(0, times = nrow(data))
      }
    }
    se <- as.matrix(as.data.frame(se))
  }
  # normalise distance matrix so that maximum distance = 1
  if (!is.null(dist_mat)) dist_mat <- dist_mat / max(dist_mat)
  # match tip ids
  tip_id <- match(data[,id], tree[[1]]$tip.label)
  # data list for stan
  sd <- list(
    N_tips = length(tree[[1]]$tip.label), # number of tips
    N_tree = N_tree,                      # number of tips
    N_obs = nrow(y),                      # number of observations
    J = length(variables),                # number of variables
    N_seg = N_seg,                        # number of segments in the tree
    node_seq = stan_node_seq,             # sequence of nodes
    parent = stan_parent,                 # parent node for each node
    ts = stan_ts,                         # amount of time since parent node
    tip = stan_tip,                       # is tip?
    effects_mat = effects_mat,            # which effects should be estimated?
    num_effects = sum(effects_mat),       # number of effects being estimated
    y = y,                                # observed data
    miss = miss,                          # are data points missing?
    tip_id = tip_id                       # tip ids
  )
  # add distance matrix if specified
  if (!is.null(dist_mat)) sd[["dist_mat"]] <- dist_mat
  # add squared standard errors if measurement_error specified
  if (!is.null(measurement_error)) sd[["se"]] <- se
  # add prior_only
  sd[["prior_only"]] <- as.numeric(prior_only)
  # produce warnings for missing data
  if (sum(miss) > 0) {
    message(
      paste0(
        "Note: Missing values (NAs) detected. These values have been included ",
        "in the Stan data list and will be imputed by the model. Set ",
        "complete_cases = TRUE to exclude taxa with missing values."
      )
    )
  }
  # return stan data
  return(sd)
}
