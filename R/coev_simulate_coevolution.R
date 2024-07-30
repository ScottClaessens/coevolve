#' Simulate the coevolution of multiple variables in discrete time steps
#'
#' @param n Number of data points in the resulting data frame.
#' @param variables A character vector of variable names (e.g.,
#'   \code{c("x","y")})
#' @param selection_matrix A numeric matrix determining the strength of
#'   selection between variables. The matrix must have a number of rows and
#'   columns equal to the number of variables and its row and column names must
#'   contain all specified variables. Each cell determines the strength of
#'   selection from the column variable to the row variable. For example, the
#'   cell on the "x" column and the "x" row indicates how much previous values
#'   of x influence future values of x (autocorrelation). By contrast, the cell
#'   on the "x" column and the "y" row indicates how much previous values of x
#'   influence future values of y (cross-lagged effect).
#' @param drift A named numeric vector specifying the strength of drift for
#'   different variables. Names must include all specified variables.
#' @param prob_split A numeric probability of a species split in any given
#'   timestep.
#'
#' @return List with dataset at final timestep (data), full simulation log
#'   (simulation), and pruned phylogenetic tree (tree).
#' @export
#'
#' @examples
#' # simulate coevolution of x and y
#' n <- 100
#' variables <- c("x","y")
#' # x -> y but not vice versa
#' selection_matrix <- matrix(
#'   c(
#'     0.95, 0.00,
#'     0.80, 0.95
#'   ),
#'   nrow = 2,
#'   byrow = TRUE,
#'   dimnames = list(variables, variables)
#' )
#' drift <- c("x" = 0.05, "y" = 0.05)
#' prob_split <- 0.05
#' # run simulation
#' sim <- coev_simulate_coevolution(n, variables, selection_matrix,
#'                                  drift, prob_split)
coev_simulate_coevolution <- function(n,
                                      variables,
                                      selection_matrix,
                                      drift,
                                      prob_split) {
  # check arguments
  # stop if n is not numeric
  if (!methods::is(n, "numeric")) {
    stop2("Argument 'n' is not numeric.")
  }
  # stop if variables is not a character vector or has length < 2
  # must also have non-reserved names
  if (!methods::is(variables, "character")) {
    stop2("Argument 'variables' is not a character vector.")
  } else if (length(variables) < 2) {
    stop2("Argument 'variables' must specify at least two variable names.")
  } else if (any(c("ts","species","parent","split") %in% variables)) {
    stop2(
      paste0(
        "Argument 'variables' uses variable names reserved internally for ",
        "simulation ('ts','species','parent','split')."
        )
      )
  }
  # stop if selection matrix is not a numeric matrix with correct dims and names
  if (!methods::is(selection_matrix, "matrix")) {
    stop2("Argument 'selection_matrix' is not a matrix.")
  } else if (!methods::is(as.vector(selection_matrix), "numeric")) {
    stop2("Argument 'selection_matrix' is not a numeric matrix.")
  } else if (!identical(as.numeric(dim(selection_matrix)),
                        rep(as.numeric(length(variables)), 2))) {
    stop2(
      paste0(
        "Argument 'selection_matrix' has number of rows or columns not equal ",
        "to the number of variables."
        )
      )
  } else if (!(all(variables %in% rownames(selection_matrix)) &
              all(variables %in% colnames(selection_matrix)))) {
    stop2(
      paste0(
        "Argument 'selection_matrix' has row or column names not equal to ",
        "variable names."
        )
      )
  }
  # stop if drift is not named numeric vector with correct dims and var names
  if (!methods::is(drift, "numeric")) {
    stop2("Argument 'drift' is not numeric.")
  } else if (length(drift) != length(variables)) {
    stop2(
      "Argument 'drift' has length different to specified number of variables."
      )
  } else if (!all(variables %in% names(drift))) {
    stop2("Argument 'drift' has names different to specified variable names.")
  }
  # stop if prob_split is not numeric of length 1
  if (!methods::is(prob_split, "numeric")) {
    stop2("Argument 'prob_split' is not numeric.")
  } else if (length(prob_split) != 1) {
    stop2("Argument 'prob_split' must be of length 1.")
  }
  # ensure selection_matrix and drift vector have names
  # in the same order as "variables" vector
  selection_matrix <- selection_matrix[variables,variables]
  drift <- drift[variables]
  if (!(identical(variables, rownames(selection_matrix)) &
        identical(variables, colnames(selection_matrix)))) {
    stop2("Selection matrix names are not equal to variable names.")
  } else if (!identical(variables, names(drift))) {
    stop2("Drift vector names are not equal to variable names.")
  }
  # begin simulation
  # first timestep
  sim <- data.frame(
    ts = 1,
    species = "t1",
    parent = NA,
    split = FALSE
    )
  for (i in variables) sim[i] <- 0
  # initial tree in text form
  tree <- "(t1:1);"
  # current timestep, species list, and number of species
  ts <- 1
  current_species <- unique(sim$species)
  current_number <- length(current_species)
  # run timesteps until there are at least n species
  while (current_number < as.integer(n)) {
    # increment timestep
    ts <- ts + 1
    # update species list and number of species
    current_species <- unique(sim$species)
    current_number <- length(current_species)
    # loop over current species
    for (sp in current_species) {
      # get previous values of variables
      prev <- c()
      for (i in variables) {
        prev[i] <- sim[sim$ts == (ts-1) & sim$species == sp, i]
      }
      # randomly generate split
      if (stats::runif(1) < prob_split) {
        # if there is a split, a new species branches off
        new_sp <- paste0("t", current_number + 1)
        # both species get new values of variables from shared parent
        new_values <- data.frame(
          ts = ts,
          species = c(sp, new_sp),
          parent = sp,
          split = TRUE
        )
        # loop over outcome variables
        for (i in variables) {
          mean <- 0
          # loop over predictor variables
          for (j in variables) mean <- mean + selection_matrix[i,j]*prev[j]
          new_values[i] <- stats::rnorm(1, mean, drift[i])
        }
        sim <- rbind(sim, new_values)
        # update current number of species
        current_number <- current_number + 1
        # update tree
        tree <- stringr::str_replace(
          string = tree,
          pattern = stringr::fixed(paste0(sp, ":")),
          replacement = paste0("(", sp, ":1,", new_sp, ":1):")
        )
      } else {
        # else if there is no split, species gets new values from itself
        new_values <- data.frame(
          ts = ts,
          species = sp,
          parent = sp,
          split = FALSE
        )
        # loop over outcome variables
        for (i in variables) {
          mean <- 0
          # loop over predictor variables
          for (j in variables) mean <- mean + selection_matrix[i,j]*prev[j]
          new_values[i] <- stats::rnorm(1, mean, drift[i])
        }
        sim <- rbind(sim, new_values)
        # update tree
        current_count <- stringr::str_extract(
          string = stringr::str_extract(
            string = tree,
            pattern = paste0(sp, ":\\d+")
            ),
          pattern = ":\\d+"
        )
        tree <- stringr::str_replace(
          string = tree,
          pattern = paste0(sp, ":\\d+"),
          replacement = paste0(sp, ":", readr::parse_number(current_count) + 1)
        )
      }
    }
  }
  # get final values
  d <- sim[sim$ts == max(sim$ts), c("species", variables)]
  # if there are more species than n, randomly prune dataset
  if (nrow(d) > n) d <- dplyr::slice_sample(d, n = as.integer(n))
  # remove rownames from d
  rownames(d) <- NULL
  # prune phylogenetic tree as well
  tree <- ape::keep.tip(ape::read.tree(text = tree), d$species)
  # return a list with data of the run and phylogenetic tree
  out <- list(
    data = d,
    simulation = sim,
    tree = tree
    )
  return(out)
}
