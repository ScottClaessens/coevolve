#' Internal helper function for checking data argument
#'
#' @srrstats {G1.4a} Non-exported function documented here
#' @srrstats {G5.2, G5.2a} Unique error messages for each input
#'
#' @description Checks the data argument for the functions coev_make_stancode(),
#'   coev_make_standata(), and coev_fit().
#'
#' @returns Error message if any of the checks fail
#'
#' @noRd
run_checks_data <- function(data) {
  # coerce data argument to data frame
  data <- try(as.data.frame(data), silent = TRUE)
  # stop if data not coercible to data frame
  #' @srrstats {G2.7} Accepts multiple tabular forms
  if (methods::is(data, "try-error") || !is.data.frame(data)) {
    stop2("Argument 'data' must be coercible to a data.frame.")
  }
  # stop if data does not contain observations
  #' @srrstats {G2.0} Assertion on number of data rows
  if (!(nrow(data) > 0L)) {
    stop2("Argument 'data' does not contain observations.")
  }
}

#' Internal helper function for checking variables argument
#'
#' @srrstats {G1.4a} Non-exported function documented here
#' @srrstats {G5.2, G5.2a} Unique error messages for each input
#'
#' @description Checks the variables argument for the functions
#'   coev_make_stancode(), coev_make_standata(), and coev_fit().
#'
#' @returns Error message if any of the checks fail
#'
#' @noRd
run_checks_variables <- function(data, variables) {
  # as data frame
  data <- as.data.frame(data)
  # stop if variables argument is not a named list
  #' @srrstats {G2.1} Assertion on type of input
  #' @srrstats {G2.8} Input to sub-function is of class 'list'
  if (!is.list(variables) || is.null(names(variables))) {
    stop2("Argument 'variables' is not a named list.")
  }
  # extract distributions and variable names from named list
  distributions <- as.character(variables)
  variables <- names(variables)
  # stop if variables are not valid column names in data
  if (!all(variables %in% colnames(data))) {
    stop2("Some variable names are not valid column names in the data.")
  }
  # stop if response distributions are not valid
  if (!all(distributions %in% c("bernoulli_logit", "ordered_logistic",
                                "poisson_softplus", "normal", "gamma_log",
                                "negative_binomial_softplus"))) {
    stop2(
      paste0(
        "Response distributions other than 'bernoulli_logit', ",
        "'ordered_logistic', 'poisson_softplus', ",
        "'negative_binomial_softplus', 'normal', and 'gamma_log' are ",
        "not yet supported."
      )
    )
  }
  # stop if not at least two variables
  #' @srrstats {G2.0} Assertion on length of input
  if (!(length(variables) >= 2)) {
    stop2("Must be at least two coevolving variables.")
  }
  # stop if data contains list columns
  #' @srrstats {G2.12} Checking for list columns
  if ("list" %in% sapply(data[, variables], class)) {
    stop2("Data must not contain list columns.")
  }
  # stop if data contains undefined values
  #' @srrstats {G2.16} Checking for undefined values
  if (any(data[, variables] == Inf, na.rm = TRUE) ||
        any(data[, variables] == -Inf, na.rm = TRUE)) {
    stop2("Data must not contain undefined values (i.e., Inf or -Inf).")
  }
  # check distributional constraints
  #' @srrstats {G2.11} The following checks ensure that data columns meet the
  #' requirements for different response distributions and are thus
  #' appropriately processed, will also ensure that columns are not lists
  for (i in seq_along(distributions)) {
    #' @srrstats {G2.15} Software does not assume non-missingness
    if (distributions[i] == "bernoulli_logit" &&
          (!is.integer(data[, variables[i]]) ||
             !all(data[, variables[i]] %in% c(0, 1, NA)))) {
      # stop if any bernoulli variables are not 0/1 integers
      stop2(
        paste0(
          "Variables following the 'bernoulli_logit' response distribution ",
          "must be integers with values of 0/1 in the data. Try using the ",
          "as.integer() function to convert variables to integers."
        )
      )
    } else if (distributions[i] == "ordered_logistic" &&
                 !is.ordered(data[, variables[i]])) {
      # stop if any ordinal variables are not ordered factors in data
      #' @srrstats {G2.5} Checking ordered factors
      stop2(
        paste0(
          "Variables following the 'ordered_logistic' response distribution ",
          "must be ordered factors in the data. Try using the as.ordered() ",
          "function to convert variables to ordered factors."
        )
      )
    } else if (distributions[i] == "poisson_softplus" &&
        (!is.integer(data[, variables[i]]) ||
           !all(data[, variables[i]] >= 0 | is.na(data[, variables[i]])))
    ) {
      # stop if any count variables are not integers greater than or equal to 0
      stop2(
        paste0(
          "Variables following the 'poisson_softplus' response distribution ",
          "must be integers greater than or equal to zero in the data. Try ",
          "using the as.integer() function to convert variables to integers."
        )
      )
    } else if (distributions[i] == "negative_binomial_softplus" &&
        (!is.integer(data[, variables[i]]) ||
           !all(data[, variables[i]] >= 0 | is.na(data[, variables[i]])))
    ) {
      # stop if any count variables are not integers greater than or equal to 0
      stop2(
        paste0(
          "Variables following the 'negative_binomial_softplus' response ",
          "distribution must be integers greater than or equal to zero in ",
          "the data. Try using the as.integer() function to convert variables ",
          "to integers."
        )
      )
    } else if (distributions[i] == "negative_binomial_softplus") {
      # stop if any negative binomial variables are not overdispersed
      #' @srrstats {G2.15} Software does not assume non-missingness (na.rm)
      if (stats::sd(data[, variables[i]], na.rm = TRUE)^2 <=
            mean(data[, variables[i]], na.rm = TRUE)) {
        stop2(
          paste0(
            "No overdispersion or potentially underdispersion for ",
            "variable '", variables[i], "' (sd^2 <= mean), do not use ",
            "the 'negative_binomial_softplus' response distribution ",
            "for this variable."
          )
        )
      }
    } else if (distributions[i] == "normal" &&
                 !is.numeric(data[, variables[i]])) {
      # stop if any normal variables are not numeric
      stop2(
        paste0(
          "Variables following the 'normal' response distribution must be ",
          "numeric in the data."
        )
      )
    } else if (distributions[i] == "gamma_log") {
      #' @srrstats {G2.15} Software does not assume non-missingness (na.rm)
      if (!is.numeric(data[, variables[i]]) ||
            !all(data[, variables[i]] > 0, na.rm = TRUE)) {
        # stop if any gamma variables are not numeric and positive
        stop2(
          paste0(
            "Variables following the 'gamma_log' response distribution must ",
            "be positive reals in the data."
          )
        )
      }
    }
  }
}

#' Internal helper function for checking id argument
#'
#' @srrstats {G1.4a} Non-exported function documented here
#' @srrstats {G5.2, G5.2a} Unique error messages for each input
#'
#' @description Checks the id argument for the functions coev_make_stancode(),
#'   coev_make_standata(), and coev_fit().
#'
#' @returns Error message if any of the checks fail
#'
#' @noRd
run_checks_id <- function(data, id) {
  # as data frame
  data <- as.data.frame(data)
  # stop if id is not a character of length one
  #' @srrstats {G2.0, G2.1, G2.2} Assertion on type and length of input
  #' @srrstats {G2.8} Input to sub-function is of class 'character'
  if (length(id) != 1 || !is.character(id)) {
    stop2("Argument 'id' must be a character of length one.")
  }
  # stop if id is not a valid column name
  #' @srrstats {G2.3, G2.3a} Permit only expected character input
  if (!(id %in% colnames(data))) {
    stop2("Argument 'id' is not a valid column name in the data.")
  }
  # stop if id in data contains missing values
  if (any(is.na(data[, id]))) {
    stop2("The id variable in the data must not contain NAs.")
  }
}

#' Internal helper function for checking tree argument
#'
#' @srrstats {G1.4a} Non-exported function documented here
#' @srrstats {G5.2, G5.2a} Unique error messages for each input
#'
#' @description Checks the tree argument for the functions coev_make_stancode(),
#'   coev_make_standata(), and coev_fit().
#'
#' @returns Error message if any of the checks fail
#'
#' @noRd
run_checks_tree <- function(data, id, tree) {
  # as data frame
  data <- as.data.frame(data)
  # stop if tree is not a phylo or multiPhylo object
  #' @srrstats {G2.1} Assertion on type of input
  #' @srrstats {G2.8} Input to sub-function is of class 'phylo' or 'multiPhylo'
  if (!(methods::is(tree, "phylo") || methods::is(tree, "multiPhylo"))) {
    stop2(
      paste0(
        "Argument 'tree' must be a phylogenetic tree object of class phylo ",
        "or multiPhylo."
      )
    )
  }
  tree <- phytools::as.multiPhylo(tree)
  for (t in seq_along(tree)) {
    # stop if tree does not have branch length information
    if (is.null(tree[[t]]$edge.length)) {
      stop2("All trees in 'tree' argument must include branch lengths.")
    }
    # stop if tree is not rooted
    if (!ape::is.rooted(tree[[t]])) {
      stop2("All trees in 'tree' argument must be rooted.")
    }
    # stop if tree contains any branch lengths not > 0
    if (!all(tree[[t]]$edge.length > 0)) {
      stop2(
        paste0(
          "All trees in 'tree' argument must have positive non-zero branch ",
          "lengths."
        )
      )
    }
    # stop if id in data does not match tree tip labels exactly
    if (!identical(sort(unique(data[, id])), sort(tree[[t]]$tip.label))) {
      stop2(
        "The id variable in the data does not match tree tip labels exactly."
      )
    }
  }
  # stop if trees have different numbers of internal nodes or branches
  if (length(unique(lapply(tree, function(x) x$Nnode))) != 1 ||
        length(unique(lapply(tree, function(x) length(x$edge.length)))) != 1) {
    stop2(
      paste0(
        "All trees in 'tree' argument must have the same number of ",
        "internal nodes and branches."
      )
    )
  }
}

#' Internal helper function for checking effects_mat argument
#'
#' @srrstats {G1.4a} Non-exported function documented here
#' @srrstats {G5.2, G5.2a} Unique error messages for each input
#'
#' @description Checks the effects_mat argument for the functions
#'   coev_make_stancode(), coev_make_standata(), and coev_fit().
#'
#' @returns Error message if any of the checks fail
#'
#' @noRd
run_checks_effects_mat <- function(variables, effects_mat) {
  # if user entered an effects matrix
  if (!is.null(effects_mat)) {
    # extract variable names from named list
    variables <- names(variables)
    # stop if effects_mat is not a matrix
    #' @srrstats {G2.1, G2.8} Assertion on matrix input
    if (!methods::is(effects_mat, "matrix")) {
      stop2("Argument 'effects_mat' must be a matrix.")
    }
    # stop if effects_mat is not logical
    #' @srrstats {G2.1} Assertion on type of input
    if (!is.logical(effects_mat)) {
      stop2("Argument 'effects_mat' must be a boolean matrix.")
    }
    # stop if effects_mat has missing row or column names
    if (is.null(rownames(effects_mat)) || is.null(colnames(effects_mat))) {
      stop2("Argument 'effects_mat' does not have valid row or column names.")
    }
    # stop if row or column names do not match variable names exactly
    if (!identical(sort(variables), sort(rownames(effects_mat))) ||
          !identical(sort(variables), sort(colnames(effects_mat)))) {
      stop2(
        paste0(
          "Row and column names for argument 'effects_mat' do not match ",
          "variable names exactly."
        )
      )
    }
    # stop if diagonal of effects_mat is ever FALSE
    for (i in variables) {
      if (!effects_mat[i, i]) {
        stop2(
          paste0(
            "Argument 'effects_mat' must specify TRUE for ",
            "all autoregressive effects."
          )
        )
      }
    }
  }
}

#' Internal helper function for checking dist_mat argument
#'
#' @srrstats {G1.4a} Non-exported function documented here
#' @srrstats {G5.2, G5.2a} Unique error messages for each input
#'
#' @description Checks the dist_mat argument for the functions
#'   coev_make_stancode(), coev_make_standata(), and coev_fit().
#'
#' @returns Error message if any of the checks fail
#'
#' @noRd
run_checks_dist_mat <- function(data, dist_mat, id) {
  # if user entered a distance matrix
  if (!is.null(dist_mat)) {
    # as data frame
    data <- as.data.frame(data)
    # stop if dist_mat is not a matrix
    #' @srrstats {G2.1, G2.8} Assertion on matrix input
    if (!methods::is(dist_mat, "matrix")) {
      stop2("Argument 'dist_mat' must be a matrix.")
    }
    # stop if dist_mat is not numeric
    if (!is.numeric(dist_mat)) {
      stop2("Argument 'dist_mat' must be a numeric matrix.")
    }
    # stop if dist_mat is not symmetric
    if (!isSymmetric(dist_mat)) {
      stop2("Argument 'dist_mat' must be a symmetric matrix.")
    }
    # stop if diagonal of dist_mat is not 0
    #' @srrstats {G2.4, G2.4b} Convert to continuous for checking
    if (!identical(as.numeric(diag(dist_mat)), rep(0, nrow(dist_mat)))) {
      stop2(
        "Argument 'dist_mat' must have zeroes on the diagonal of the matrix."
      )
    }
    # stop if dist_mat has missing row or column names
    if (is.null(rownames(dist_mat)) || is.null(colnames(dist_mat))) {
      stop2("Argument 'dist_mat' does not have valid row or column names.")
    }
    # stop if row and column names do not match tip labels exactly
    if (!identical(sort(unique(data[, id])), sort(rownames(dist_mat))) ||
          !identical(sort(unique(data[, id])), sort(colnames(dist_mat)))) {
      stop2(
        paste0(
          "Row and column names for argument 'dist_mat' do not match tree ",
          "tip labels exactly."
        )
      )
    }
  }
}

#' Internal helper function for checking dist_cov argument
#'
#' @srrstats {G1.4a} Non-exported function documented here
#' @srrstats {G5.2, G5.2a} Unique error messages for each input
#'
#' @description Checks the dist_cov argument for the functions
#'   coev_make_stancode(), coev_make_standata(), and coev_fit().
#'
#' @returns Error message if any of the checks fail
#'
#' @noRd
run_checks_dist_cov <- function(dist_cov) {
  # stop if dist_cov is not a character string
  #' @srrstats {G2.1, G2.8} Assertion on character input
  if (!methods::is(dist_cov, "character")) {
    stop2("Argument 'dist_cov' is not a character string.")
  }
  # stop if dist_cov is not of length 1
  #' @srrstats {G2.0, G2.2} Assertion on length of input
  if (length(dist_cov) != 1) {
    stop2("Argument 'dist_cov' is not of length 1.")
  }
  # stop if specified dist_cov is not supported
  #' @srrstats {G2.3, G2.3a} Permit only expected character input
  if (!(dist_cov %in% c("exp_quad", "exponential", "matern32"))) {
    stop2(
      paste0(
        "Argument 'dist_cov' currently only supports 'exp_quad', ",
        "'exponential', and 'matern32'."
      )
    )
  }
}

#' Internal helper function for checking measurement_error argument
#'
#' @srrstats {G1.4a} Non-exported function documented here
#' @srrstats {G5.2, G5.2a} Unique error messages for each input
#'
#' @description Checks the measurement_error argument for the functions
#'   coev_make_stancode(), coev_make_standata(), and coev_fit().
#'
#' @returns Error message if any of the checks fail
#'
#' @noRd
run_checks_measurement_error <- function(data, variables, measurement_error) {
  # if user declares measurement error
  if (!is.null(measurement_error)) {
    # as data frame
    data <- as.data.frame(data)
    # stop if measurement_error argument is not a named list
    #' @srrstats {G2.1, G2.8} Assertion on named list input
    if (!is.list(measurement_error) || is.null(names(measurement_error))) {
      stop2("Argument 'measurement_error' is not a named list.")
    }
    # extract standard error columns and variable names from named list
    #' @srrstats {G2.4, G2.4b} Convert to character for checking
    error_columns <- as.character(measurement_error)
    error_variables <- names(measurement_error)
    # extract distributions and variable names from named list
    #' @srrstats {G2.4, G2.4b} Convert to character for checking
    distributions <- as.character(variables)
    variables <- names(variables)
    # stop if any variables are not valid or normally distributed
    if (!all(error_variables %in% variables[distributions == "normal"])) {
      stop2(
        paste0(
          "Argument 'measurement_error' contains variables that were not ",
          "declared as normally-distributed variables in the model."
        )
      )
    }
    # stop if any columns are not actual columns in the dataset
    if (!all(error_columns %in% colnames(data))) {
      stop2(
        paste0(
          "Argument 'measurement_error' refers to measurement error columns ",
          "that are not valid column names in the data."
        )
      )
    }
    # loop over all error columns
    for (i in seq_along(error_columns)) {
      # stop if error column is non-numeric or contains any non-positive-reals
      #' @srrstats {G2.15} Software does not assume non-missingness (na.rm)
      if (!is.numeric(data[[error_columns[i]]]) ||
            !all(data[[error_columns[i]]] >= 0, na.rm = TRUE)) {
        stop2(
          paste0(
            "Standard errors in measurement error columns must be zero or ",
            "positive reals."
          )
        )
      }
      # check if error column contains NAs where there are
      # observed values for the focal variable
      check <- is.na(
        data[[error_columns[i]]][!is.na(data[[error_variables[i]]])]
      )
      if (any(check)) {
        stop2(
          paste0(
            "Standard errors in measurement error columns must not be NA ",
            "in rows where there is observed data for the focal variable."
          )
        )
      }
    }
  }
}

#' Internal helper function for checking prior argument
#'
#' @srrstats {G1.4a} Non-exported function documented here
#' @srrstats {G5.2, G5.2a} Unique error messages for each input
#'
#' @description Checks the prior argument for the functions
#'   coev_make_stancode(), coev_make_standata(), and coev_fit().
#'
#' @returns Error message if any of the checks fail
#'
#' @noRd
run_checks_prior <- function(prior) {
  # if user entered a list of priors
  if (!is.null(prior)) {
    # stop if prior not a list
    #' @srrstats {G2.1, G2.8} Assertion on list input
    if (!methods::is(prior, "list")) {
      stop2("Argument 'prior' is not a list.")
    }
    # stop if prior not a named list
    if (is.null(names(prior))) {
      stop2("Argument 'prior' is not a named list.")
    }
    # stop if prior names not allowed
    if (!all(names(prior) %in% c("b", "eta_anc", "A_offdiag", "A_diag",
                                 "L_R", "Q_sigma", "c", "phi", "shape",
                                 "sigma_dist", "rho_dist", "sigma_residual",
                                 "L_residual"))) {
      stop2(
        paste0(
          "Argument 'prior' list contains names that are not allowed. Please ",
          "use only the following names: 'b', 'eta_anc', 'A_offdiag', ",
          "'A_diag', 'L_R', 'Q_sigma', 'c', 'phi', 'shape', 'sigma_dist', ",
          "'rho_dist', 'sigma_residual', and 'L_residual'"
        )
      )
    }
    # stop if prior names contains duplicates
    if (length(unique(names(prior))) != length(names(prior))) {
      stop2("Argument 'prior' contains duplicate names.")
    }
  }
}

#' Internal helper function for checking all arguments
#'
#' @srrstats {G1.4a} Non-exported function documented here
#' @srrstats {G5.2, G5.2a} Unique error messages for each input
#'
#' @description Checks all arguments for the functions coev_make_stancode(),
#'   coev_make_standata(), and coev_fit().
#'
#' @returns Error message if any of the checks fail
#'
#' @noRd
run_checks <- function(data, variables, id, tree, effects_mat, complete_cases,
                       dist_mat, dist_cov, measurement_error, prior, scale,
                       estimate_correlated_drift, estimate_residual, log_lik,
                       prior_only) {
  # run checks from previous functions
  run_checks_data(data)
  run_checks_variables(data, variables)
  run_checks_id(data, id)
  run_checks_tree(data, id, tree)
  run_checks_effects_mat(variables, effects_mat)
  run_checks_dist_mat(data, dist_mat, id)
  run_checks_dist_cov(dist_cov)
  run_checks_measurement_error(data, variables, measurement_error)
  run_checks_prior(prior)
  # check that other arguments are logical of length one
  #' @srrstats {G2.0, G2.1, G2.2, G2.8} Assertion on length and type of inputs
  if (!is.logical(complete_cases) || length(complete_cases) != 1) {
    stop2("Argument 'complete_cases' must be a logical of length one.")
  }
  if (!is.logical(scale) || length(scale) != 1) {
    stop2("Argument 'scale' must be a logical of length one.")
  }
  if (!is.logical(estimate_correlated_drift) ||
        length(estimate_correlated_drift) != 1) {
    stop2(
      "Argument 'estimate_correlated_drift' must be a logical of length one."
    )
  }
  if (!is.logical(estimate_residual) || length(estimate_residual) != 1) {
    stop2("Argument 'estimate_residual' must be a logical of length one.")
  }
  if (!is.logical(log_lik) || length(log_lik) != 1) {
    stop2("Argument 'log_lik' must be a logical of length one.")
  }
  if (!is.logical(prior_only) || length(prior_only) != 1) {
    stop2("Argument 'prior_only' must be a logical of length one.")
  }
}

#' Internal helper function for producing errors
#'
#' @srrstats {G1.4a} Non-exported function documented here
#'
#' @description Stops execution of the current expression and executes an
#'   error action, without the call becoming part of the error message.
#'
#' @returns Error message
#'
#' @noRd
stop2 <- function(...) {
  stop(..., call. = FALSE)
}

#' Internal helper function for producing warnings
#'
#' @srrstats {G1.4a} Non-exported function documented here
#'
#' @description Generates a warning message that corresponds to its argument,
#'   without the call becoming part of the warning message.
#'
#' @returns Warning message
#'
#' @noRd
warning2 <- function(...) {
  warning(..., call. = FALSE)
}
