#' Get pinned Python dependency specs for the JAX backend
#'
#' @description Returns a character vector of pip-style version specs
#'   read from \code{inst/python/requirements.txt}.
#'
#' @returns Character vector, e.g. \code{c("jax==0.5.3", ...)}.
#'
#' @noRd
jax_python_deps <- function() {
  req_file <- system.file(
    "python", "requirements.txt", package = "coevolve"
  )
  if (!nzchar(req_file)) {
    return(c("jax", "numpyro", "nutpie"))
  }
  lines <- readLines(req_file, warn = FALSE)
  lines <- trimws(lines)
  lines[nzchar(lines) & !startsWith(lines, "#")]
}

#' Ensure JAX and numpyro are available via reticulate
#'
#' @description Checks that jax, numpyro, and nutpie can be imported
#'   at the pinned versions from \code{inst/python/requirements.txt}.
#'
#' @returns Logical. TRUE if all packages are available, FALSE otherwise.
#'
#' @noRd
check_jax_available <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    return(FALSE)
  }
  tryCatch({
    reticulate::py_require(jax_python_deps())
    reticulate::import("jax", convert = FALSE)
    reticulate::import("numpyro", convert = FALSE)
    reticulate::import("nutpie", convert = FALSE)
    TRUE
  }, error = function(e) {
    FALSE
  })
}

#' Stop if JAX backend is not available
#'
#' @returns Error message if JAX is not available.
#'
#' @noRd
stop_if_jax_not_available <- function() {
  if (!check_jax_available()) {
    stop2(
      "JAX backend is not available. Install pinned dependencies with:\n",
      "  pip install -r ",
      system.file("python", "requirements.txt", package = "coevolve"),
      "\nOr let reticulate manage it automatically by ensuring ",
      "RETICULATE_PYTHON is configured."
    )
  }
}

#' Load the coev_jax_model Python module via reticulate
#'
#' @param convert Logical. Passed to \code{reticulate::import_from_path}.
#'
#' @returns Python module object with \code{build_nutpie_model} function.
#'
#' @noRd
load_jax_model_module <- function(convert = FALSE) {
  py_file <- system.file("python", "coev_jax_model.py", package = "coevolve")
  if (!nzchar(py_file)) {
    stop2("inst/python/coev_jax_model.py not found (broken package install?).")
  }
  reticulate::import_from_path(
    tools::file_path_sans_ext(basename(py_file)),
    path    = normalizePath(dirname(py_file), winslash = "/", mustWork = TRUE),
    convert = convert
  )
}

#' Convert R data list to Python dict for JAX model
#'
#' @param data_list Named list of model data.
#'
#' @returns Python dict (via reticulate).
#'
#' @noRd
convert_r_to_python_data_jax <- function(data_list) {
  np <- reticulate::import("numpy", convert = FALSE)

  integer_vars <- c(
    "node_seq", "parent", "tip", "effects_mat",
    "tip_id", "N_tips", "N_tree", "N_obs", "J", "N_seg",
    "num_effects", "prior_only", "miss", "length_index",
    "tip_to_seg", "N_unique_lengths",
    "n_levels", "max_level_size", "root_ids",
    "level_node_ids", "level_parent_ids",
    "level_length_idx", "level_is_internal",
    "level_drift_idx", "level_sizes",
    "NBgp"
  )

  python_data <- list()
  for (name in names(data_list)) {
    value <- data_list[[name]]
    is_int <- name %in% integer_vars || is.integer(value)

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
      dtype <- if (is_int) "int32" else "float64"
      value <- if (is_int) as.integer(value) else as.double(value)
      python_data[[name]] <- np$array(value, dtype = dtype)
    }
  }
  reticulate::r_to_py(python_data)
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
#'
#' @noRd
prior_spec_from_stan <- function(prior_str, constraint = "none") {
  p <- parse_stan_prior(prior_str)
  list(
    dist       = p$dist,
    args       = as.list(as.numeric(p$args)),
    constraint = constraint
  )
}

#' Embed model config flags into a Stan data dict
#'
#' Merges the config list returned by \code{coev_make_model_config()} into the
#' JAX-ready data dict returned by \code{standata_to_jax()}, producing
#' a single dict that the model builder can consume.
#'
#' @param sd_jax Named list from \code{standata_to_jax()}.
#' @param config Named list from \code{coev_make_model_config()}.
#'
#' @returns Extended named list suitable for model building.
#'
#' @noRd
embed_model_config <- function(sd_jax, config) {
  for (nm in names(config)) sd_jax[[nm]] <- config[[nm]]
  sd_jax
}

#' Convert Stan data list to JAX-friendly format
#'
#' @description Takes the output of \code{coev_make_standata()} and converts
#'   it for use by the JAX backend: 1-based indices become 0-based,
#'   observed means for count variables are precomputed, and ordinal variable
#'   metadata is attached.
#'
#' @param sd Named list from \code{coev_make_standata()}.
#' @param distributions Character vector of response distributions.
#'
#' @returns Modified named list suitable for the JAX model.
#'
#' @noRd
standata_to_jax <- function(sd, distributions) {
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
#' Groups tree nodes by topological depth so the model graph can process
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
                                length_index_0, # nolint
                                N_tree, N_seg) { # nolint
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

  # Map each segment position to its z_drift index (= segment_pos - 1,
  # matching Stan's z_drift[t, i-1] indexing for i in 2:N_seg).
  seg_drift_slot <- matrix(0L, N_tree, N_seg)
  for (t in seq_len(N_tree)) {
    for (i in 2:N_seg) {
      seg_drift_slot[t, i] <- i - 2L
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
        lvl_drift_idx[t, l, k]   <- seg_drift_slot[t, s]
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
