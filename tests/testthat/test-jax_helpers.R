# Unit tests for internal helpers in R/jax_helpers.R that are not directly
# exercised by the user-facing API.

# Helper: build the 0-indexed inputs that compute_tree_levels expects from
# a phylo or multiPhylo, mirroring what standata_to_jax does.
make_tree_inputs <- function(tree) {
  tip_labels <- if (inherits(tree, "multiPhylo")) {
    tree[[1]]$tip.label
  } else {
    tree$tip.label
  }
  n <- length(tip_labels)
  d <- data.frame(
    id = tip_labels,
    x = rnorm(n),
    y = rnorm(n)
  )
  sd <- coev_make_standata(
    data = d,
    variables = list(x = "normal", y = "normal"),
    id = "id",
    tree = tree
  )
  list(
    node_seq_0     = sd$node_seq - 1L,
    parent_0       = sd$parent - 1L,
    tip_0          = sd$tip,           # already 0/1 indicator
    length_index_0 = sd$length_index - 1L,
    N_tree         = sd$N_tree,
    N_seg          = sd$N_seg
  )
}

call_levels <- function(inp) {
  coevolve:::compute_tree_levels(
    inp$node_seq_0, inp$parent_0, inp$tip_0,
    inp$length_index_0, inp$N_tree, inp$N_seg
  )
}

test_that("compute_tree_levels: output shapes are consistent", {
  withr::with_seed(1, tree <- ape::rcoal(5))
  inp <- make_tree_inputs(tree)
  lvl <- call_levels(inp)

  expect_named(
    lvl,
    c("n_levels", "max_level_size", "root_ids",
      "level_node_ids", "level_parent_ids", "level_length_idx",
      "level_is_internal", "level_drift_idx", "level_sizes")
  )

  level_dims <- c(inp$N_tree, lvl$n_levels, lvl$max_level_size)
  expect_equal(dim(lvl$level_node_ids), level_dims)
  expect_equal(dim(lvl$level_parent_ids), level_dims)
  expect_equal(dim(lvl$level_length_idx), level_dims)
  expect_equal(dim(lvl$level_is_internal), level_dims)
  expect_equal(dim(lvl$level_drift_idx), level_dims)
  expect_equal(dim(lvl$level_sizes), c(inp$N_tree, lvl$n_levels))
  expect_equal(length(lvl$root_ids), inp$N_tree)
})

test_that("compute_tree_levels: each non-root seg placed exactly once", {
  withr::with_seed(2, tree <- ape::rcoal(6))
  inp <- make_tree_inputs(tree)
  lvl <- call_levels(inp)

  # Sum of real (unpadded) sizes for each tree must equal N_seg - 1
  # (every segment except the root traversal entry).
  for (t in seq_len(inp$N_tree)) {
    expect_equal(sum(lvl$level_sizes[t, ]), inp$N_seg - 1L)
  }
})

test_that("compute_tree_levels: root_id matches first traversal entry", {
  withr::with_seed(3, tree <- ape::rcoal(5))
  inp <- make_tree_inputs(tree)
  lvl <- call_levels(inp)
  for (t in seq_len(inp$N_tree)) {
    expect_equal(lvl$root_ids[t], inp$node_seq_0[t, 1])
  }
})

test_that("compute_tree_levels: padding slots point at the tree's root_id", {
  withr::with_seed(4, tree <- ape::rcoal(6))
  inp <- make_tree_inputs(tree)
  lvl <- call_levels(inp)

  for (t in seq_len(inp$N_tree)) {
    rid <- lvl$root_ids[t]
    for (l in seq_len(lvl$n_levels)) {
      sz <- lvl$level_sizes[t, l]
      if (sz < lvl$max_level_size) {
        pad <- (sz + 1L):lvl$max_level_size
        expect_true(all(lvl$level_node_ids[t, l, pad] == rid))
        expect_true(all(lvl$level_parent_ids[t, l, pad] == rid))
      }
    }
  }
})

test_that("compute_tree_levels: multiphylo gets N_tree dimension right", {
  withr::with_seed(5, {
    t1 <- ape::rcoal(5)
    t2 <- ape::rcoal(5)
    t2$tip.label <- t1$tip.label
  })
  trees <- structure(list(t1, t2), class = "multiPhylo")
  inp <- make_tree_inputs(trees)
  lvl <- call_levels(inp)

  expect_equal(inp$N_tree, 2L)
  expect_equal(length(lvl$root_ids), 2L)
  expect_equal(dim(lvl$level_node_ids)[1], 2L)

  # The two trees can have different topologies, so each must independently
  # account for all its own non-root segments.
  for (t in seq_len(inp$N_tree)) {
    expect_equal(sum(lvl$level_sizes[t, ]), inp$N_seg - 1L)
  }
})

test_that("compute_tree_levels: drift_idx covers 0..(N_seg-2) exactly once", {
  # The Stan model expects every non-root segment to consume exactly one
  # z_drift slot, with positions 0..(N_seg-2) covered exactly once.
  withr::with_seed(6, tree <- ape::rcoal(7))
  inp <- make_tree_inputs(tree)
  lvl <- call_levels(inp)

  for (t in seq_len(inp$N_tree)) {
    used_slots <- integer(0)
    for (l in seq_len(lvl$n_levels)) {
      sz <- lvl$level_sizes[t, l]
      if (sz > 0L) {
        used_slots <- c(used_slots, lvl$level_drift_idx[t, l, seq_len(sz)])
      }
    }
    expect_setequal(used_slots, seq.int(0L, inp$N_seg - 2L))
  }
})
