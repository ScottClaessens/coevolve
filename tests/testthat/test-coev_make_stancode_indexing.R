# Test suite for verifying correct array indexing in generated Stan code
#
# This test suite catches bugs where segment indices are confused with node
# indices. Key insight: VCV_tips and L_VCV_tips are indexed by NODE NUMBER
# (1 to N_seg), not by segment traversal order.
#
# The tip_to_seg mapping converts tip node numbers to segment indices,
# so using tip_to_seg to index VCV_tips would access the WRONG elements.

test_that("VCV_tips and L_VCV_tips do not use tip_to_seg for access", {
  # Generate Stan code for a simple model

  sc <- coev_make_stancode(
    data = authority$data,
    variables = list(
      political_authority = "ordered_logistic",
      religious_authority = "ordered_logistic"
    ),
    id = "language",
    tree = authority$phylogeny
  )

  # VCV_tips should NOT be accessed via tip_to_seg

  # (tip_to_seg returns segment indices, but VCV_tips is indexed by node number)
  expect_false(
    grepl("VCV_tips\\[t,\\s*tip_to_seg", sc),
    info = "VCV_tips should not use tip_to_seg indexing"
  )
  expect_false(
    grepl("L_VCV_tips\\[t,\\s*tip_to_seg", sc),
    info = "L_VCV_tips should not use tip_to_seg indexing"
  )
})

test_that("VCV_tips storage uses node_seq during tree traversal", {
  sc <- coev_make_stancode(
    data = authority$data,
    variables = list(
      political_authority = "ordered_logistic",
      religious_authority = "ordered_logistic"
    ),
    id = "language",
    tree = authority$phylogeny
  )

  # When storing VCV_tips during tree traversal, should use node_seq
  expect_true(
    grepl("VCV_tips\\[t, node_seq\\[t,", sc),
    info = "VCV_tips storage should use node_seq indexing"
  )
  expect_true(
    grepl("L_VCV_tips\\[t, node_seq\\[t,", sc),
    info = "L_VCV_tips storage should use node_seq indexing"
  )
})

test_that("segment-indexed arrays use correct traversal indexing", {
  sc <- coev_make_stancode(
    data = authority$data,
    variables = list(
      political_authority = "ordered_logistic",
      religious_authority = "ordered_logistic"
    ),
    id = "language",
    tree = authority$phylogeny
  )

  # z_drift should be accessed by segment index (i-1 in the loop)
  expect_true(
    grepl("z_drift\\[t, i-1\\]", sc),
    info = "z_drift should use segment indexing"
  )

  # length_index should be accessed by segment index
  expect_true(
    grepl("length_index\\[t, i\\]", sc),
    info = "length_index should use segment indexing"
  )
})

test_that("tip_id maps observations to tips correctly in model block", {
  sc <- coev_make_stancode(
    data = authority$data,
    variables = list(
      political_authority = "ordered_logistic",
      religious_authority = "ordered_logistic"
    ),
    id = "language",
    tree = authority$phylogeny
  )

  # eta should be accessed via tip_id for observations
  expect_true(
    grepl("eta\\[t,\\s*tip_id\\[i\\]\\]", sc),
    info = "eta should use tip_id to access tip states for observations"
  )
})

test_that("measurement error model uses correct VCV_tips indexing", {
  skip_on_cran()

  # Create test data with measurement error
  d <- authority$data
  d$x <- rnorm(nrow(d))
  d$y <- rnorm(nrow(d))
  d$x_se <- rexp(nrow(d), 5)
  d$y_se <- rexp(nrow(d), 5)

  sc <- coev_make_stancode(
    data = d,
    variables = list(x = "normal", y = "normal"),
    id = "language",
    tree = authority$phylogeny,
    measurement_error = list(x = "x_se", y = "y_se")
  )

  # With measurement error, VCV_tips access should still not use tip_to_seg
  expect_false(
    grepl("VCV_tips\\[t,\\s*tip_to_seg", sc),
    info = "Measurement error model should not use tip_to_seg for VCV_tips"
  )

  # Should use tip_id for observation-level access in model block
  expect_true(
    grepl("VCV_tips\\[t,\\s*tip_id\\[i\\]\\]", sc),
    info = "Measurement error model should use tip_id for VCV_tips access"
  )
})

test_that("normal distribution model uses correct L_VCV_tips indexing", {
  skip_on_cran()

  # Create test data with normal variables (no measurement error)
  d <- authority$data
  d$x <- rnorm(nrow(d))
  d$y <- rnorm(nrow(d))

  sc <- coev_make_stancode(
    data = d,
    variables = list(x = "normal", y = "normal"),
    id = "language",
    tree = authority$phylogeny
  )

  # L_VCV_tips access should not use tip_to_seg

  expect_false(
    grepl("L_VCV_tips\\[t,\\s*tip_to_seg", sc),
    info = "Normal model should not use tip_to_seg for L_VCV_tips"
  )

  # Should use tip_id for observation-level access
  expect_true(
    grepl("L_VCV_tips\\[t,\\s*tip_id\\[i\\]\\]", sc),
    info = "Normal model should use tip_id for L_VCV_tips access"
  )
})

test_that("generated quantities uses correct indexing for yrep", {
  sc <- coev_make_stancode(
    data = authority$data,
    variables = list(
      political_authority = "ordered_logistic",
      religious_authority = "ordered_logistic"
    ),
    id = "language",
    tree = authority$phylogeny
  )

  # In generated quantities, tip_id should be used to link observations to tips
  expect_true(
    grepl("terminal_drift_rep\\[t,\\s*tip_id\\[i\\]\\]", sc),
    info = "Generated quantities should use tip_id for terminal_drift_rep"
  )
})
