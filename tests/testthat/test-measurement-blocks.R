test_that("plain measurement block lists are normalized with singleton unblocked variables", {
  blocks <- as_measurement_blocks(
    measurement_blocks = list(
      alpha = c("a1", "a2"),
      beta = c("b1", "b2")
    ),
    variables = c("a1", "a2", "b1", "b2", "x")
  )

  expect_equal(names(blocks), c("alpha", "beta", "x"))
  expect_equal(blocks$x, "x")
})

test_that("cscore-style higher-order composites are merged into final blocks", {
  composite_list <- list(
    lower = list(
      consistency_interest = c("gs2", "gs3"),
      perseverance_effort = c("gs1", "gs4"),
      extraversion = c("e1", "e2")
    ),
    higher = list(
      grit = c("consistency_interest", "perseverance_effort")
    ),
    order = c("consistency_interest", "perseverance_effort", "extraversion", "grit")
  )

  blocks <- as_measurement_blocks(
    measurement_blocks = composite_list,
    variables = c("gs1", "gs2", "gs3", "gs4", "e1", "e2", "x")
  )

  expect_equal(sort(blocks$grit), sort(c("gs1", "gs2", "gs3", "gs4")))
  expect_equal(blocks$extraversion, c("e1", "e2"))
  expect_equal(blocks$x, "x")
  expect_false("consistency_interest" %in% names(blocks))
  expect_false("perseverance_effort" %in% names(blocks))
})

test_that("measurement block validation rejects missing and duplicated items", {
  expect_error(
    as_measurement_blocks(
      measurement_blocks = list(alpha = c("a1", "missing")),
      variables = c("a1", "a2")
    ),
    "not present in `variables`"
  )

  expect_error(
    as_measurement_blocks(
      measurement_blocks = list(alpha = c("a1", "a2"), beta = c("a2", "b1")),
      variables = c("a1", "a2", "b1")
    ),
    "more than one final block"
  )

  expect_error(
    as_measurement_blocks(
      measurement_blocks = list(
        lower = list(alpha = c("a1", "a2")),
        higher = list(total = c("alpha", "missing_lower"))
      ),
      variables = c("a1", "a2")
    ),
    "references lower-order composites"
  )
})

test_that("block shrinkage matrix preserves within-block entries and shrinks between-block entries", {
  blocks <- list(alpha = c("a1", "a2"), beta = c("b1", "b2"))
  mat <- make_block_shrinkage_matrix(
    variables = c("a1", "a2", "b1", "b2"),
    measurement_blocks = blocks,
    between_block_shrink = 0.25
  )

  expect_equal(unname(diag(mat)), rep(1, 4))
  expect_equal(mat["a1", "a2"], 1)
  expect_equal(mat["b1", "b2"], 1)
  expect_equal(mat["a1", "b1"], 0.25)
  expect_equal(rownames(mat), c("a1", "a2", "b1", "b2"))
})

test_that("block shrinkage attenuates between-block correlations only", {
  cor_mat <- matrix(
    0.8,
    nrow = 4,
    ncol = 4,
    dimnames = list(c("a1", "a2", "b1", "b2"), c("a1", "a2", "b1", "b2"))
  )
  diag(cor_mat) <- 1
  shrink <- make_block_shrinkage_matrix(
    variables = c("a1", "a2", "b1", "b2"),
    measurement_blocks = list(alpha = c("a1", "a2"), beta = c("b1", "b2")),
    between_block_shrink = 0.5
  )

  out <- apply_block_shrinkage(cor_mat, shrink)

  expect_equal(unname(diag(out)), rep(1, 4))
  expect_gt(abs(out["a1", "a2"]), abs(out["a1", "b1"]))
})

test_that("within-block correlations can be restored from an observed reference", {
  variables <- c("a1", "a2", "b1", "b2")
  blocks <- list(alpha = c("a1", "a2"), beta = c("b1", "b2"))

  target <- matrix(
    0.1,
    nrow = 4,
    ncol = 4,
    dimnames = list(variables, variables)
  )
  diag(target) <- 1
  target["a1", "a2"] <- target["a2", "a1"] <- 0.2
  target["b1", "b2"] <- target["b2", "b1"] <- 0.3

  reference <- target
  reference["a1", "a2"] <- reference["a2", "a1"] <- 0.75
  reference["b1", "b2"] <- reference["b2", "b1"] <- 0.65

  out <- restore_within_block_correlations(
    cor_mat = target,
    reference_cor = reference,
    measurement_blocks = blocks
  )

  expect_equal(out["a1", "a2"], 0.75)
  expect_equal(out["b1", "b2"], 0.65)
  expect_equal(out["a1", "b1"], 0.1)
  expect_equal(unname(diag(out)), rep(1, 4))
})
