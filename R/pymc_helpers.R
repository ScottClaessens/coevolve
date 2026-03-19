#' Ensure PyMC and PyTensor are available via reticulate
#'
#' @description Uses py_require to install (if needed) and then verifies that
#'   pymc and pytensor can be imported.
#'
#' @param mlx Logical. If TRUE, also require mlx for Apple Silicon GPU.
#'
#' @returns Logical. TRUE if pymc is available, FALSE otherwise.
#'
#' @noRd
check_pymc_available <- function(mlx = FALSE) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    return(FALSE)
  }
  tryCatch({
    pkgs <- c("pymc", "pytensor", "arviz")
    if (mlx) pkgs <- c(pkgs, "mlx")
    reticulate::py_require(pkgs)
    reticulate::import("pymc", convert = FALSE)
    reticulate::import("pytensor", convert = FALSE)
    TRUE
  }, error = function(e) {
    FALSE
  })
}

#' Stop if PyMC is not available via reticulate
#'
#' @description Stops with an informative error if pymc or pytensor cannot be
#'   imported via reticulate.
#'
#' @returns Error message if pymc is not available.
#'
#' @noRd
stop_if_pymc_not_available <- function() {
  if (!check_pymc_available()) {
    stop2(
      "PyMC/PyTensor is not available. Please install with: ",
      "pip install pymc pytensor arviz. ",
      "For Apple Silicon GPU support via MLX, also install mlx. ",
      "You may also need to configure reticulate to find your Python ",
      "installation."
    )
  }
}

#' Convert a Stan prior string to PyMC distribution code
#'
#' @description Maps a Stan prior specification to the equivalent PyMC
#'   distribution constructor call. Handles parameter constraints by selecting
#'   appropriate truncated or half distributions.
#'
#' @param prior_str Stan prior string (e.g., "std_normal()", "normal(0, 2)").
#' @param constraint One of "none", "lower_zero", "upper_zero".
#' @param shape_str Optional Python expression for shape argument
#'   (e.g., "J", "(N_tree, J)").
#'
#' @returns Character string of PyMC distribution code (without name argument).
#'
#' @noRd
stan_prior_to_pymc <- function(prior_str, constraint = "none",
                               shape_str = NULL) {
  parsed <- parse_stan_prior(prior_str)

  shape_arg <- if (!is.null(shape_str)) {
    sprintf(", shape=%s", shape_str)
  } else {
    ""
  }

  pymc_str <- switch(parsed$dist,
    "std_normal" = {
      if (constraint == "upper_zero") {
        sprintf("pm.TruncatedNormal(%%s, mu=0.0, sigma=1.0, upper=0.0%s)",
                shape_arg)
      } else if (constraint == "lower_zero") {
        sprintf("pm.HalfNormal(%%s, sigma=1.0%s)", shape_arg)
      } else {
        sprintf("pm.Normal(%%s, mu=0.0, sigma=1.0%s)", shape_arg)
      }
    },
    "normal" = {
      mu <- parsed$args[1]
      sigma <- parsed$args[2]
      if (constraint == "upper_zero") {
        sprintf("pm.TruncatedNormal(%%s, mu=%s, sigma=%s, upper=0.0%s)",
                mu, sigma, shape_arg)
      } else if (constraint == "lower_zero") {
        sprintf("pm.TruncatedNormal(%%s, mu=%s, sigma=%s, lower=0.0%s)",
                mu, sigma, shape_arg)
      } else {
        sprintf("pm.Normal(%%s, mu=%s, sigma=%s%s)", mu, sigma, shape_arg)
      }
    },
    "exponential" = {
      sprintf("pm.Exponential(%%s, lam=%s%s)", parsed$args[1], shape_arg)
    },
    "gamma" = {
      sprintf("pm.Gamma(%%s, alpha=%s, beta=%s%s)",
              parsed$args[1], parsed$args[2], shape_arg)
    },
    "lkj_corr_cholesky" = {
      sprintf("__LKJ__%s", parsed$args[1])
    },
    stop2("Unknown prior distribution for PyMC mapping: ", parsed$dist)
  )
  pymc_str
}

#' Convert Stan data list to PyMC-friendly format
#'
#' @description Takes the output of \code{coev_make_standata()} and converts
#'   it for use by the PyMC backend: 1-based indices become 0-based,
#'   observed means for count variables are precomputed, and ordinal variable
#'   metadata is attached.
#'
#' @param sd Named list from \code{coev_make_standata()}.
#' @param distributions Character vector of response distributions.
#'
#' @returns Modified named list suitable for the PyMC model.
#'
#' @noRd
standata_to_pymc <- function(sd, distributions) {
  sd$node_seq <- sd$node_seq - 1L
  sd$parent <- sd$parent - 1L
  sd$tip_id <- sd$tip_id - 1L
  sd$length_index <- sd$length_index - 1L
  sd$tip_to_seg <- sd$tip_to_seg - 1L

  count_idx <- which(distributions %in%
                       c("poisson_softplus", "negative_binomial_softplus"))
  obs_means <- rep(0.0, length(distributions))
  for (j in count_idx) {
    non_miss <- sd$miss[, j] == 0
    obs_means[j] <- mean(sd$y[non_miss, j])
  }
  sd$obs_means <- obs_means

  nb_idx <- which(distributions == "negative_binomial_softplus")
  inv_overdisp <- rep(0.0, length(distributions))
  for (j in nb_idx) {
    non_miss <- sd$miss[, j] == 0
    obs_j <- sd$y[non_miss, j]
    m <- mean(obs_j)
    v <- stats::var(obs_j)
    inv_overdisp[j] <- ifelse(v > m, m^2 / (v - m), m)
  }
  sd$inv_overdisp <- inv_overdisp
  sd
}

#' Convert R data list to Python dict for PyMC
#'
#' @param data_list Named list of model data.
#'
#' @returns Python dict (via reticulate).
#'
#' @noRd
convert_r_to_python_data_pymc <- function(data_list) {
  np <- reticulate::import("numpy", convert = FALSE)

  integer_vars <- c("node_seq", "parent", "tip", "effects_mat",
                     "tip_id", "N_tips", "N_tree", "N_obs", "J", "N_seg",
                     "num_effects", "prior_only", "miss", "length_index",
                     "tip_to_seg", "N_unique_lengths")

  python_data <- list()
  for (name in names(data_list)) {
    value <- data_list[[name]]
    is_int <- name %in% integer_vars || is.integer(value)

    if (is.matrix(value) || is.array(value)) {
      dtype <- if (is_int) "int32" else "float64"
      if (is_int) {
        value <- array(as.integer(value), dim = dim(value))
      } else {
        value <- array(as.double(value), dim = dim(value))
      }
      python_data[[name]] <- np$array(value, dtype = dtype)
    } else if (length(value) == 1) {
      python_data[[name]] <- if (is_int) as.integer(value) else as.double(value)
    } else {
      dtype <- if (is_int) "int32" else "float64"
      value <- if (is_int) as.integer(value) else as.double(value)
      python_data[[name]] <- np$array(value, dtype = dtype)
    }
  }
  reticulate::r_to_py(python_data)
}
