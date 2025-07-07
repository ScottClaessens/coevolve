test_that("coev_simulate_coevolution() produces expected errors", {
  n <- 10
  variables <- c("x", "y")
  selection_matrix <- matrix(c(0.95, 0.00, 0.80, 0.95), nrow = 2,
                             dimnames = list(variables, variables))
  drift <- c("x" = 0.05, "y" = 0.05)
  prob_split <- 0.05
  # expect the following errors
  #' @srrstats {G5.2, G5.2b} Test all error messages
  expect_error(
    # n is not numeric
    coev_simulate_coevolution(n = "fail", variables, selection_matrix,
                              drift, prob_split),
    "Argument 'n' is not numeric.",
    fixed = TRUE
  )
  expect_error(
    # variables are not characters
    coev_simulate_coevolution(n, variables = c(1, 2), selection_matrix,
                              drift, prob_split),
    "Argument 'variables' is not a character vector.",
    fixed = TRUE
  )
  expect_error(
    # less than two variables
    coev_simulate_coevolution(n, variables = "fail", selection_matrix,
                              drift, prob_split),
    "Argument 'variables' must specify at least two variable names.",
    fixed = TRUE
  )
  expect_error(
    # variables use reserved names
    coev_simulate_coevolution(
      n, variables = c("ts", "species", "parent", "split"),
      selection_matrix, drift, prob_split
    ),
    paste0(
      "Argument 'variables' uses variable names reserved internally for ",
      "simulation ('ts','species','parent','split')."
    ),
    fixed = TRUE
  )
  expect_error(
    # selection matrix is not a matrix
    coev_simulate_coevolution(n, variables, selection_matrix = "fail",
                              drift, prob_split),
    "Argument 'selection_matrix' is not a matrix.",
    fixed = TRUE
  )
  expect_error(
    # selection matrix is not numeric
    coev_simulate_coevolution(n, variables, selection_matrix = matrix("fail"),
                              drift, prob_split),
    "Argument 'selection_matrix' is not a numeric matrix.",
    fixed = TRUE
  )
  expect_error(
    # selection matrix has incorrect dimensions
    coev_simulate_coevolution(n, variables, selection_matrix = matrix(1:4),
                              drift, prob_split),
    paste0(
      "Argument 'selection_matrix' has number of rows or columns not equal to ",
      "the number of variables."
    ),
    fixed = TRUE
  )
  expect_error(
    # selection matrix has incorrect row/col names
    coev_simulate_coevolution(n, variables,
                              selection_matrix = matrix(1:4, nrow = 2),
                              drift, prob_split),
    paste0(
      "Argument 'selection_matrix' has row or column names not equal to ",
      "variable names."
    ),
    fixed = TRUE
  )
  expect_error(
    # drift is not numeric
    coev_simulate_coevolution(n, variables, selection_matrix,
                              drift = "fail", prob_split),
    "Argument 'drift' is not numeric.",
    fixed = TRUE
  )
  expect_error(
    # drift has incorrect length
    coev_simulate_coevolution(n, variables, selection_matrix,
                              drift = 1, prob_split),
    "Argument 'drift' has length different to specified number of variables.",
    fixed = TRUE
  )
  expect_error(
    # drift has incorrect names
    coev_simulate_coevolution(n, variables, selection_matrix,
                              drift = c(1, 2), prob_split),
    "Argument 'drift' has names different to specified variable names.",
    fixed = TRUE
  )
  expect_error(
    # prob_split is not numeric
    coev_simulate_coevolution(n, variables, selection_matrix,
                              drift, prob_split = "fail"),
    "Argument 'prob_split' is not numeric.",
    fixed = TRUE
  )
  expect_error(
    # prob_split is not length 1
    coev_simulate_coevolution(n, variables, selection_matrix,
                              drift, prob_split = c(1, 2)),
    "Argument 'prob_split' must be of length 1.",
    fixed = TRUE
  )
  expect_error(
    # prob_split is not between 0 and 1
    coev_simulate_coevolution(n, variables, selection_matrix,
                              drift, prob_split = -0.05),
    "Argument 'prob_split' must be between 0 and 1.",
    fixed = TRUE
  )
  expect_error(
    # prob_split is not between 0 and 1
    coev_simulate_coevolution(n, variables, selection_matrix,
                              drift, prob_split = 1.05),
    "Argument 'prob_split' must be between 0 and 1.",
    fixed = TRUE
  )
  expect_error(
    # intercepts is not numeric
    coev_simulate_coevolution(n, variables, selection_matrix,
                              drift, prob_split, intercepts = "fail"),
    "Argument 'intercepts' is not numeric.",
    fixed = TRUE
  )
  expect_error(
    # intercepts is not correct length
    coev_simulate_coevolution(n, variables, selection_matrix,
                              drift, prob_split, intercepts = 0),
    paste0(
      "Argument 'intercepts' has length different to specified number of ",
      "variables."
    ),
    fixed = TRUE
  )
  expect_error(
    # intercepts does not have correct names
    coev_simulate_coevolution(n, variables, selection_matrix,
                              drift, prob_split,
                              intercepts = c("a" = 0, "b" = 0)),
    "Argument 'intercepts' has names different to specified variable names.",
    fixed = TRUE
  )
  expect_error(
    # ancestral_states is not numeric
    coev_simulate_coevolution(n, variables, selection_matrix,
                              drift, prob_split, ancestral_states = "fail"),
    "Argument 'ancestral_states' is not numeric.",
    fixed = TRUE
  )
  expect_error(
    # ancestral_states is not correct length
    coev_simulate_coevolution(n, variables, selection_matrix,
                              drift, prob_split, ancestral_states = 0),
    paste0(
      "Argument 'ancestral_states' has length different to specified number ",
      "of variables."
    ),
    fixed = TRUE
  )
  expect_error(
    # ancestral_states does not have correct names
    coev_simulate_coevolution(n, variables, selection_matrix,
                              drift, prob_split,
                              ancestral_states = c("a" = 0, "b" = 0)),
    paste0(
      "Argument 'ancestral_states' has names different to specified variable ",
      "names."
    ),
    fixed = TRUE
  )
})

test_that("coev_simulate_coevolution() returns named list", {
  n <- 10
  variables <- c("x", "y")
  selection_matrix <- matrix(c(0.95, 0.80, 0.00, 0.95), nrow = 2,
                             # test when name order does not match var order
                             dimnames = list(c("x", "y"), c("y", "x")))
  drift <- c("x" = 0.05, "y" = 0.05)
  prob_split <- 0.05
  sim <- coev_simulate_coevolution(n, variables, selection_matrix,
                                   drift, prob_split)
  # expect named list of length 3
  expect_no_error(sim)
  expect_type(sim, "list")
  expect_length(sim, 3)
  # first element is data frame
  expect_true(is.data.frame(sim$data))
  # which should have num columns = number of variables + 1
  expect_true(ncol(sim$data) == length(variables) + 1)
  # and should have no missing values
  #' @srrstats {G5.3} Ensure no missing values
  expect_equal(sum(is.na(sim$data)), 0)
  # second element is data frame
  expect_true(is.data.frame(sim$simulation))
  # third element is phylo
  expect_true(class(sim$tree) == "phylo")
})

test_that("Intercepts and anc states work in coev_simulate_coevolution()", {
  n <- 10
  variables <- c("x", "y")
  selection_matrix <- matrix(c(0.95, 0.00, 0.00, 0.95), nrow = 2,
                             dimnames = list(c("x", "y"), c("x", "y")))
  drift <- c("x" = 0.05, "y" = 0.05)
  prob_split <- 0.05
  # setting intercepts
  expect_no_error(
    coev_simulate_coevolution(
      n, variables, selection_matrix,
      drift, prob_split,
      intercepts = c("x" = 1, "y" = 4)
    )
  )
  # setting ancestral states
  expect_no_error(
    coev_simulate_coevolution(
      n, variables, selection_matrix,
      drift, prob_split,
      ancestral_states = c("x" = 1, "y" = 4)
    )
  )
  # setting both intercepts and ancestral states
  expect_no_error(
    coev_simulate_coevolution(
      n, variables, selection_matrix,
      drift, prob_split,
      intercepts = c("x" = 1, "y" = 4),
      ancestral_states = c("x" = 1, "y" = 4)
    )
  )
})
