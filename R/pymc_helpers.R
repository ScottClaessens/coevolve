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

#' Parse a Stan-style prior string
#'
#' @param prior_str Stan prior string such as \code{"normal(0, 2)"}.
#'
#' @returns A list with \code{dist} and \code{args}.
#'
#' @noRd
parse_stan_prior <- function(prior_str) {
  if (!is.character(prior_str) || length(prior_str) != 1 || is.na(prior_str)) {
    stop2("Prior must be a single non-missing character string.")
  }

  prior_str <- stringr::str_trim(prior_str)
  matches <- stringr::str_match(prior_str, "^([A-Za-z0-9_]+)\\((.*)\\)$")

  if (is.na(matches[1, 1])) {
    stop2("Could not parse Stan prior string: ", prior_str)
  }

  dist <- matches[1, 2]
  arg_string <- stringr::str_trim(matches[1, 3])
  args <- if (!nzchar(arg_string)) {
    character(0)
  } else {
    stringr::str_split(arg_string, "\\s*,\\s*", simplify = FALSE)[[1]]
  }

  list(dist = dist, args = args)
}

#' Convert a Stan prior string to a Python-compatible prior spec dict
#'
#' @param prior_str Stan prior string (e.g., "std_normal()", "normal(0, 2)").
#' @param constraint One of "none", "lower_zero", "upper_zero".
#'
#' @returns Named list with \code{dist}, \code{args}, and \code{constraint}.
#'   Passed to data_dict and consumed by \code{CoevPymcModel._make_rv()}.
#'
#' @noRd
prior_spec_from_stan <- function(prior_str, constraint = "none") {
  p <- parse_stan_prior(prior_str)
  # Use as.list() so reticulate always converts args to a Python list, not a
  # scalar float (which would make `for a in args` fail in _make_rv).
  list(
    dist       = p$dist,
    args       = as.list(as.numeric(p$args)),
    constraint = constraint
  )
}

#' Embed PyMC config flags into a Stan data dict
#'
#' Merges the config list returned by \code{coev_make_pymc()} into the
#' PyMC-ready data dict returned by \code{standata_to_pymc()}, producing
#' a single dict that \code{CoevPymcModel().build()} can consume.
#'
#' @param sd_pymc Named list from \code{standata_to_pymc()}.
#' @param config Named list from \code{coev_make_pymc()}.
#'
#' @returns Extended named list suitable for \code{CoevPymcModel().build()}.
#'
#' @noRd
embed_pymc_config <- function(sd_pymc, config) {
  for (nm in names(config)) sd_pymc[[nm]] <- config[[nm]]
  sd_pymc
}

#' Load the coev_pymc_model Python module via reticulate
#'
#' Imports \code{inst/python/coev_pymc_model.py} from the installed package
#' and returns the module object.
#'
#' @param convert Logical. Passed to \code{reticulate::import_from_path}.
#'
#' @returns Python module object with \code{build_model} and
#'   \code{CoevPymcModel} attributes.
#'
#' @noRd
load_pymc_model_module <- function(convert = FALSE) {
  py_file <- system.file("python", "coev_pymc_model.py", package = "coevolve")
  if (!nzchar(py_file)) {
    stop2("inst/python/coev_pymc_model.py not found (broken package install?).")
  }
  reticulate::import_from_path(
    tools::file_path_sans_ext(basename(py_file)),
    path    = normalizePath(dirname(py_file), winslash = "/", mustWork = TRUE),
    convert = convert
  )
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

  lvl <- compute_tree_levels(
    node_seq_0    = sd$node_seq,
    parent_0      = sd$parent,
    tip_0         = sd$tip,
    length_index_0 = sd$length_index,
    N_tree        = sd$N_tree,
    N_seg         = sd$N_seg
  )
  sd$n_levels        <- lvl$n_levels
  sd$max_level_size  <- lvl$max_level_size
  sd$root_ids        <- lvl$root_ids
  sd$level_node_ids  <- lvl$level_node_ids
  sd$level_parent_ids <- lvl$level_parent_ids
  sd$level_length_idx <- lvl$level_length_idx
  sd$level_is_internal <- lvl$level_is_internal
  sd$level_drift_idx <- lvl$level_drift_idx
  sd$level_sizes     <- lvl$level_sizes

  sd
}

#' Pre-compute level-batched tree traversal data
#'
#' Groups tree nodes by topological depth so the PyMC graph can process
#' each level as a single batched operation instead of one op per node.
#'
#' @param node_seq_0 Integer matrix (N_tree, N_seg), 0-indexed node IDs in
#'   traversal order.
#' @param parent_0 Integer matrix (N_tree, N_seg), 0-indexed parent node IDs.
#' @param tip_0 Integer matrix (N_tree, N_seg), 1 = tip, 0 = internal.
#' @param length_index_0 Integer matrix (N_tree, N_seg), 0-indexed index into
#'   unique_lengths.
#' @param N_tree Integer, number of trees.
#' @param N_seg Integer, number of segments per tree.
#'
#' @returns Named list with padded level-batched arrays.
#'
#' @noRd
compute_tree_levels <- function(node_seq_0, parent_0, tip_0,
                                length_index_0, N_tree, N_seg) {
  all_depths <- matrix(0L, N_tree, N_seg)
  root_ids <- integer(N_tree)

  for (t in seq_len(N_tree)) {
    ns <- node_seq_0[t, ]
    pa <- parent_0[t, ]
    root_ids[t] <- ns[1]

    node_to_step <- integer(N_seg)
    node_to_step[ns + 1L] <- seq_len(N_seg) - 1L

    depth <- integer(N_seg)
    for (i in 2:N_seg) {
      parent_step_r <- node_to_step[pa[i] + 1L] + 1L
      depth[i] <- depth[parent_step_r] + 1L
    }
    all_depths[t, ] <- depth
  }

  n_levels_total <- max(all_depths) + 1L
  n_nonroot <- n_levels_total - 1L

  max_level_size <- 0L
  for (t in seq_len(N_tree)) {
    for (l in seq_len(n_nonroot)) {
      max_level_size <- max(max_level_size, sum(all_depths[t, ] == l))
    }
  }

  internal_slot <- matrix(-1L, N_tree, N_seg)
  for (t in seq_len(N_tree)) {
    cnt <- 0L
    for (i in 2:N_seg) {
      if (tip_0[t, i] == 0L) {
        internal_slot[t, i] <- cnt
        cnt <- cnt + 1L
      }
    }
  }

  lvl_node_ids    <- array(0L, dim = c(N_tree, n_nonroot, max_level_size))
  lvl_parent_ids  <- array(0L, dim = c(N_tree, n_nonroot, max_level_size))
  lvl_length_idx  <- array(0L, dim = c(N_tree, n_nonroot, max_level_size))
  lvl_is_internal <- array(0L, dim = c(N_tree, n_nonroot, max_level_size))
  lvl_drift_idx   <- array(0L, dim = c(N_tree, n_nonroot, max_level_size))
  lvl_sizes       <- matrix(0L, N_tree, n_nonroot)

  for (t in seq_len(N_tree)) {
    ns <- node_seq_0[t, ]
    pa <- parent_0[t, ]
    tp <- tip_0[t, ]
    li <- length_index_0[t, ]
    rid <- root_ids[t]

    for (l in seq_len(n_nonroot)) {
      steps <- which(all_depths[t, ] == l)
      sz <- length(steps)
      lvl_sizes[t, l] <- sz

      for (k in seq_len(sz)) {
        s <- steps[k]
        lvl_node_ids[t, l, k]    <- ns[s]
        lvl_parent_ids[t, l, k]  <- pa[s]
        lvl_length_idx[t, l, k]  <- li[s]
        lvl_is_internal[t, l, k] <- if (tp[s] == 0L) 1L else 0L
        lvl_drift_idx[t, l, k]   <- if (tp[s] == 0L) internal_slot[t, s] else 0L
      }
      if (sz < max_level_size) {
        for (k in (sz + 1L):max_level_size) {
          lvl_node_ids[t, l, k]   <- rid
          lvl_parent_ids[t, l, k] <- rid
        }
      }
    }
  }

  list(
    n_levels        = n_nonroot,
    max_level_size  = max_level_size,
    root_ids        = root_ids,
    level_node_ids  = lvl_node_ids,
    level_parent_ids = lvl_parent_ids,
    level_length_idx = lvl_length_idx,
    level_is_internal = lvl_is_internal,
    level_drift_idx = lvl_drift_idx,
    level_sizes     = lvl_sizes
  )
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
                     "tip_to_seg", "N_unique_lengths",
                     "n_levels", "max_level_size", "root_ids",
                     "level_node_ids", "level_parent_ids",
                     "level_length_idx", "level_is_internal",
                     "level_drift_idx", "level_sizes")

  python_data <- list()
  for (name in names(data_list)) {
    value <- data_list[[name]]
    is_int <- name %in% integer_vars || is.integer(value)

    # Character vectors, named lists (e.g. prior_specs, distributions): pass
    # through reticulate's default conversion without numeric coercion.
    if (is.character(value) || (is.list(value) && !is.null(names(value)))) {
      python_data[[name]] <- reticulate::r_to_py(value)
    } else if (is.matrix(value) || is.array(value)) {
      dtype <- if (is_int) "int32" else "float64"
      value <- if (is_int) {
        array(as.integer(value), dim = dim(value))
      } else {
        array(as.double(value), dim = dim(value))
      }
      python_data[[name]] <- np$array(value, dtype = dtype)
    } else if (length(value) == 1) {
      python_data[[name]] <- if (is_int) as.integer(value) else as.double(value)
    } else {
      # Handles length-0 and length > 1 numeric/integer vectors
      dtype <- if (is_int) "int32" else "float64"
      value <- if (is_int) as.integer(value) else as.double(value)
      python_data[[name]] <- np$array(value, dtype = dtype)
    }
  }
  reticulate::r_to_py(python_data)
}
