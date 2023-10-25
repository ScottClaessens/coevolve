# helper function for checking arguments
run_checks <- function(data, variables, id, tree) {
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
  if (!all(distributions %in% c("bernoulli_logit", "ordered_logistic", "poisson_log"))) {
    stop2("Response distributions other than 'bernoulli_logit', 'ordered_logistic', and 'poisson_log' are not yet supported.")
  }
  # stop if not at least two variables
  if (!(length(variables) >= 2)) {
    stop2("Must be at least two coevolving variables.")
  }
  # stop if any binary variables are not 0/1 integers
  for (i in 1:length(distributions)) {
    if (distributions[i] == "bernoulli_logit" & (!is.integer(data[,variables[i]]) | !all(data[,variables[i]] %in% 0:1))) {
      stop2("Variables following the 'bernoulli_logit' response distribution must be integers with values of 0/1 in the data.")
    }
  }
  # stop if any ordinal variables are not ordered factors in data
  for (i in 1:length(distributions)) {
    if (distributions[i] == "ordered_logistic" & !is.ordered(data[,variables[i]])) {
      stop2("Variables following the 'ordered_logistic' response distribution must be ordered factors in the data.")
    }
  }
  # stop if any count variables are not integers greater than or equal to 0
  for (i in 1:length(distributions)) {
    if (distributions[i] == "poisson_log" & (!is.integer(data[,variables[i]]) | !all(data[,variables[i]] >= 0))) {
      stop2("Variables following the 'poisson_log' response distribution must be integers greater than or equal to zero in the data.")
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
    stop2("Argument 'id' must be an phylogenetic tree object of class phylo.")
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
}

# helper function for producing errors
stop2 <- function(...) {
  stop(..., call. = FALSE)
}

# helper function for producing warnings
warning2 <- function(...) {
  warning(..., call. = FALSE)
}
