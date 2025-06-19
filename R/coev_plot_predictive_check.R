#' Plot predictive check from a fitted \code{coevfit} object
#'
#' Plot posterior (or prior) predictive checks for variables from a fitted
#' \code{coevfit} object, with help from the \pkg{bayesplot} package.
#'
#' @param object An object of class \code{coevfit}
#' @param variables If NULL (default), the function returns a list of plots for
#'   all coevolving variables from the model. Otherwise, a character vector
#'   declaring the variables to be included in the list of plots.
#' @param ndraws An integer indicating the number of draws to return. The
#'   default and maximum number of draws is the size of the posterior sample.
#' @param tree_id An integer indicating the tree ID to use when making
#'   posterior predictions. Set to \code{NULL} by default, which will use draws
#'   from every tree, integrating phylogenetic uncertainty.
#'
#' @returns A list of \code{ggplot} objects
#'
#' @author Scott Claessens \email{scott.claessens@@gmail.com}, Erik Ringen
#'   \email{erikjacob.ringen@@uzh.ch}
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
#' # plot predictive checks for all variables
#' coev_plot_predictive_check(fit)
#' }
#'
#' @export
coev_plot_predictive_check <- function(object, variables = NULL,
                                       ndraws = NULL, tree_id = NULL) {
  # stop if object is not of class coevfit
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
      # stop if variables is not a character string
      stop2(
        paste0(
          "Argument 'variables' must be a character string or a ",
          "vector of character strings."
        )
      )
    } else if (!all(variables %in% names(object$variables))) {
      # stop if variables not included in the fitted model
      stop2(
        paste0(
          "Argument 'variables' contains variable names that are not ",
          "included in the fitted model."
        )
      )
    }
  }
  # stop if ndraws is not a single integer between 1 and the total num draws
  if (!is.null(ndraws)) {
    if (!is.numeric(ndraws)) {
      stop2("Argument 'ndraws' must be numeric.")
    } else if (!all(as.integer(ndraws) == ndraws) | length(ndraws) != 1) {
      stop2("Argument 'ndraws' must be a single integer.")
    } else if (ndraws < 1 | ndraws > nrow(object$fit$draws())) {
      stop2(
        "Argument 'ndraws' must be between 1 and the total number of draws."
        )
    }
  }
  # stop if tree_id is not a single integer between 1 and the total num trees
  if (!is.null(tree_id)) {
    if (!is.numeric(tree_id)) {
      stop2("Argument 'tree_id' must be numeric.")
    } else if (!all(as.integer(tree_id) == tree_id) | length(tree_id) != 1) {
      stop2("Argument 'tree_id' must be a single integer.")
    } else if (tree_id < 1 | tree_id > object$stan_data$N_tree) {
      stop2(
        "Argument 'tree_id' must be between 1 and the total number of trees."
      )
    }
  }
  # get posterior predictions
  post <- posterior::as_draws_rvars(object$fit$draws(variables = "yrep"))
  # get draws ids
  if (is.null(ndraws)) {
    draws_ids <- 1:nrow(object$fit$draws())
  } else {
    draws_ids <- sample(1:nrow(object$fit$draws()), size = ndraws)
  }
  # get tree ids
  if (is.null(tree_id)) {
    tree_ids <- 1:object$stan_data$N_tree
  } else{
    tree_ids <- tree_id
  }
  # function for choosing plot type
  choose_plot_type <- function(variable) {
    # response distribution
    resp_dist <- object$variables[[variable]]
    # choose plot type based on response distribution
    if (resp_dist %in% c("normal", "poisson_softplus",
                         "negative_binomial_softplus", "gamma_log")) {
      bayesplot::ppc_dens_overlay
    } else if (resp_dist %in% c("bernoulli_logit", "ordered_logistic")) {
      bayesplot::ppc_bars
    }
  }
  # list to populate
  out <- list()
  # loop over variables
  if (is.null(variables)) {
    variables <- names(object$variables)
  }
  for (variable in variables) {
    # variable index
    var_id <- which(names(object$variables) == variable)
    # get y and yrep
    y <- object$stan_data$y[,variable]
    miss <- object$stan_data$miss[,variable]
    y[miss == 1] <- NA
    # combine draws from multiple trees, if applicable
    yrep <- posterior::as_draws_matrix(
      post$yrep[tree_ids[1],,var_id]
      )[draws_ids,]
    if (length(tree_ids) > 1) {
      for (i in 2:length(tree_ids)) {
        yrep <-
          rbind(
            yrep,
            posterior::as_draws_matrix(
              post$yrep[tree_ids[i],,var_id]
              )[draws_ids,]
            )
        }
    }
    # remove missing data for this variable only
    if (any(is.na(y))) {
      yrep <- yrep[,!is.na(y)]
      y <- y[!is.na(y)]
    }
    # posterior predictive check
    pp <-
      bayesplot::pp_check(
        object = as.numeric(y),
        yrep = yrep,
        choose_plot_type(variable)
      ) +
      ggplot2::ggtitle(
        paste0("Variable = '", variable, "'")
        ) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank()
        )
    # add to plot list
    out[[variable]] <- pp
  }
  return(out)
}
