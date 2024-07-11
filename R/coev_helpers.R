# helper function for checking arguments
run_checks <- function(data, variables, id, tree, effects_mat,
                       dist_mat, prior, prior_only) {
  # coerce data argument to data frame
  data <- try(as.data.frame(data), silent = TRUE)
  # stop if data not coercible to data frame
  if (methods::is(data, "try-error") | !is.data.frame(data)) {
    stop2("Argument 'data' must be coercible to a data.frame.")
  }
  # stop if data does not contain observations
  if (!(nrow(data) > 0L)) {
    stop2("Argument 'data' does not contain observations.")
  }
  # stop if variables argument is not a named list
  if (!is.list(variables) | is.null(names(variables))) {
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
                                "poisson_softplus", "normal", "lognormal",
                                "negative_binomial_softplus"))) {
    stop2(
      paste0(
        "Response distributions other than 'bernoulli_logit', ",
        "'ordered_logistic', 'poisson_softplus', 'normal', 'lognormal', and ",
        "'negative_binomial_softplus' are not yet supported."
        )
      )
  }
  # stop if not at least two variables
  if (!(length(variables) >= 2)) {
    stop2("Must be at least two coevolving variables.")
  }
  # stop if any bernoulli variables are not 0/1 integers
  for (i in 1:length(distributions)) {
    if (distributions[i] == "bernoulli_logit" &
        (!is.integer(data[,variables[i]]) |
         !all(data[,variables[i]] %in% 0:1))) {
      stop2(
        paste0(
          "Variables following the 'bernoulli_logit' response distribution ",
          "must be integers with values of 0/1 in the data."
          )
        )
    }
  }
  # stop if any ordinal variables are not ordered factors in data
  for (i in 1:length(distributions)) {
    if (distributions[i] == "ordered_logistic" &
        !is.ordered(data[,variables[i]])) {
      stop2(
        paste0(
          "Variables following the 'ordered_logistic' response distribution ",
          "must be ordered factors in the data."
          )
        )
    }
  }
  # stop if any count variables are not integers greater than or equal to 0
  for (i in 1:length(distributions)) {
    if (distributions[i] == "poisson_softplus" &
        (!is.integer(data[,variables[i]]) | !all(data[,variables[i]] >= 0))) {
      stop2(
        paste0(
          "Variables following the 'poisson_softplus' response distribution ",
          "must be integers greater than or equal to zero in the data."
          )
        )
    }
    if (distributions[i] == "negative_binomial_softplus" &
        (!is.integer(data[,variables[i]]) | !all(data[,variables[i]] >= 0))) {
      stop2(
        paste0(
          "Variables following the 'negative_binomial_softplus' response ",
          "distribution must be integers greater than or equal to zero in ",
          "the data."
        )
      )
    }
  }
  # stop if any negative binomial variables are not overdispersed
  for (i in 1:length(distributions)) {
    if (distributions[i] == "negative_binomial_softplus") {
      # if variance <= mean
      if (sd(data[,variables[i]])^2 <= mean(data[,variables[i]])) {
        stop2(
          paste0(
            "No overdispersion or potentially underdispersion for ",
            "variable '", variables[i], "' (sd^2 <= mean), do not use ",
            "the 'negative_binomial_softplus' response distribution ",
            "for this variable."
          )
        )
      }
    }
  }
  # stop if any normal variables are not numeric
  for (i in 1:length(distributions)) {
    if (distributions[i] == "normal" & !is.numeric(data[,variables[i]])) {
      stop2(
        paste0(
          "Variables following the 'normal' response distribution must be ",
          "numeric in the data."
          )
        )
    }
  }
  # stop if any lognormal variables are not numeric or are
  # equal to or less than zero
  for (i in 1:length(distributions)) {
    if (distributions[i] == "lognormal" & (!is.numeric(data[,variables[i]]) |
                                           any(data[,variables[i]] <= 0))) {
      stop2(
        paste0(
          "Variables following the 'lognormal' response distribution must be ",
          "numeric in the data and must be greater than zero."
        )
      )
    }
  }
  # stop if id is not a character of length one
  if (length(id) != 1 | !is.character(id)) {
    stop2("Argument 'id' must be a character of length one.")
  }
  # stop if id is not a valid column name
  if (!(id %in% colnames(data))) {
    stop2("Argument 'id' is not a valid column name in the data.")
  }
  # stop if tree is not a phylo object
  if (!methods::is(tree, "phylo")) {
    stop2("Argument 'tree' must be an phylogenetic tree object of class phylo.")
  }
  # stop if id in data does not match tree tip labels exactly
  if (!identical(sort(data[,id]), sort(tree$tip.label))) {
    stop2("The id variable in the data does not match tree tip labels exactly.")
  }
  # stop if id in data contains missing values
  if (any(is.na(data[,id]))) {
    stop2("The id variable in the data must not contain NAs.")
  }
  # stop if coevolving variables contain missing data
  if (any(is.na(data[,variables]))) {
    stop2("Coevolving variables in the data must not contain NAs.")
  }
  # if user entered an effects matrix
  if (!is.null(effects_mat)) {
    # stop if effects_mat is not a matrix
    if (!methods::is(effects_mat, "matrix")) {
      stop2("Argument 'effects_mat' must be a matrix.")
    }
    # stop if effects_mat is not logical
    if (!is.logical(effects_mat)) {
      stop2("Argument 'effects_mat' must be a boolean matrix.")
    }
    # stop if effects_mat has missing row or column names
    if (is.null(rownames(effects_mat)) | is.null(colnames(effects_mat))) {
      stop2("Argument 'effects_mat' does not have valid row or column names.")
    }
    # stop if row or column names do not match variable names exactly
    if (!identical(sort(variables), sort(rownames(effects_mat))) |
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
      if (!effects_mat[i,i]) {
        stop2(
          paste0(
            "Argument 'effects_mat' must specify TRUE for ",
            "all autoregressive effects."
          )
        )
      }
    }
  }
  # if user entered a geographic distance matrix
  if (!is.null(dist_mat)) {
    # stop if dist_mat is not a matrix
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
    if (!identical(as.numeric(diag(dist_mat)), rep(0, nrow(dist_mat)))) {
      stop2(
        "Argument 'dist_mat' must have zeroes on the diagonal of the matrix."
        )
    }
    # stop if dist_mat has missing row or column names
    if (is.null(rownames(dist_mat)) | is.null(colnames(dist_mat))) {
      stop2("Argument 'dist_mat' does not have valid row or column names.")
    }
    # stop if row and column names do not match tip labels exactly
    if (!identical(sort(data[,id]), sort(rownames(dist_mat))) |
        !identical(sort(data[,id]), sort(colnames(dist_mat)))) {
      stop2(
        paste0(
          "Row and column names for argument 'dist_mat' do not match tree ",
          "tip labels exactly."
          )
        )
    }
  }
  # if user entered a list of priors
  if (!is.null(prior)) {
    # stop if prior not a list
    if (!methods::is(prior, "list")) {
      stop2("Argument 'prior' is not a list.")
    }
    # stop if prior not a named list
    if (is.null(names(prior))) {
      stop2("Argument 'prior' is not a named list.")
    }
    # stop if prior names not allowed
    if (!all(names(prior) %in% c("b", "eta_anc", "A_offdiag", "A_diag",
                                 "Q_diag", "c", "phi", "sigma_dist",
                                 "rho_dist"))) {
      stop2(
        paste0(
          "Argument 'prior' list contains names that are not allowed. Please ",
          "use only the following names: 'b', 'eta_anc', 'A_offdiag', ",
          "'A_diag', 'Q_diag', 'c', 'phi', 'sigma_dist', and 'rho_dist'"
          )
        )
    }
    # stop if prior names contains duplicates
    if (length(unique(names(prior))) != length(names(prior))) {
      stop2("Argument 'prior' contains duplicate names.")
    }
  }
  # stop if prior_only is not logical
  if (!is.logical(prior_only)) {
    stop2("Argument 'prior_only' is not logical.")
  }
}

# helper function for producing errors
stop2 <- function(...) {
  stop(..., call. = FALSE)
}

# helper function for producing warnings
warning2 <- function(...) {
  warning(..., call. = FALSE)
}
