# helper function for checking arguments
run_checks <- function(data, variables, id, tree) {
  # coerce data argument to data frame
  data <- try(as.data.frame(data), silent = TRUE)
  # stop if data not coercible to data frame
  if (methods::is(data, "try-error") | !is.data.frame(data)) {
    stop("Argument 'data' must be coercible to a data.frame.")
  }
  # stop if data does not contain observations
  if (!(nrow(data) > 0L)) {
    stop("Argument 'data' does not contain observations.")
  }
  # stop if variables argument is not a named list
  if (!is.list(variables) | is.null(names(variables))) {
    stop("Argument 'variables' is not a named list.")
  }
  # extract distributions and variable names from named list
  distributions <- as.character(variables)
  variables <- names(variables)
  # stop if variables are not valid column names in data
  if (!all(variables %in% colnames(data))) {
    stop("Some variable names are not valid column names in the data.")
  }
  # stop if response distributions are not valid
  if (!all(distributions %in% c("bernoulli_logit", "ordered_logistic"))) {
    stop("Response distributions other than 'bernoulli_logit' and 'ordered_logistic' are not yet supported.")
  }
  # stop if not at least two variables
  if (!(length(variables) >= 2)) {
    stop("Must be at least two coevolving variables.")
  }
  # stop if any binary variables are not 0/1 integers
  for (i in 1:length(distributions)) {
    if (distributions[i] == "bernoulli_logit" & (!is.integer(data[,variables[i]]) | !all(data[,variables[i]] %in% 0:1))) {
      stop("Variables following the 'bernoulli_logit' response distribution must be integers with values of 0/1 in the data.")
    }
  }
  # stop if any ordinal variables are not ordered factors in data
  for (i in 1:length(distributions)) {
    if (distributions[i] == "ordered_logistic" & !is.ordered(data[,variables[i]])) {
      stop("Variables following the 'ordered_logistic' response distribution must be ordered factors in the data.")
    }
  }
  # stop if id is not a character of length one
  if (length(id) != 1 | !is.character(id)) {
    stop("Argument 'id' must be a character of length one.")
  }
  # stop if id is not a valid column name
  if (!(id %in% colnames(data))) {
    stop("Argument 'id' is not a valid column name in the data.")
  }
  # stop if id in data does not match tree tip labels exactly
  if (!identical(sort(data[,id]), sort(tree$tip.label))) {
    stop("The id variable in the data does not match tree tip labels exactly.")
  }
  # stop if tree is not a phylo object
  if (!methods::is(tree, "phylo")) {
    stop("Argument 'id' must be an phylogenetic tree object of class phylo.")
  }
}
