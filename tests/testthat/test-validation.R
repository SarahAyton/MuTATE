features <- c("age", "sex", "biomarker")
outcomes <- c("response", "tumor_size", "ae_count", "OS_definition_time_status")
outcome_defs <- c("Cat", "Cont", "Count", "Surv")

test_that("MTPart errors clearly on a missing feature column", {
  expect_error(
    MTPart(c(features, "nonexistent_col"), outcomes, outcome_defs, mutate_example),
    "features.*not found.*nonexistent_col"
  )
})

test_that("MTPart errors clearly on a missing outcome column", {
  expect_error(
    MTPart(features, c(outcomes, "nonexistent_outcome"), c(outcome_defs, "Cat"), mutate_example),
    "outcomes.*not found.*nonexistent_outcome"
  )
})

test_that("MTPart errors when outcome_defs length does not match outcomes", {
  expect_error(
    MTPart(features, outcomes, outcome_defs[1:2], mutate_example),
    "outcome_defs.*same length"
  )
})

test_that("MTPart errors on an invalid outcome_defs value", {
  expect_error(
    MTPart(features, outcomes, c("Cat", "Cont", "Count", "NotAType"), mutate_example),
    "outcome_defs.*NotAType"
  )
})

test_that("MTPart errors when a column appears in both features and outcomes", {
  expect_error(
    MTPart(c(features, "response"), outcomes, outcome_defs, mutate_example),
    "both.*features.*outcomes"
  )
})

test_that("MTPart errors on non-data-frame input", {
  expect_error(
    MTPart(features, outcomes, outcome_defs, as.matrix(mutate_example[, c(features, outcomes)])),
    "data frame"
  )
})

test_that("MTPart errors on empty data", {
  expect_error(
    MTPart(features, outcomes, outcome_defs, mutate_example[0, ]),
    "zero rows"
  )
})

test_that("MTPart errors on a mismatched wt length", {
  expect_error(
    MTPart(features, outcomes, outcome_defs, mutate_example, wt = c(1, 2)),
    "wt.*one weight per outcome"
  )
})

test_that("MTPart errors on invalid depth/nodesize/cp/continuous", {
  expect_error(MTPart(features, outcomes, outcome_defs, mutate_example, depth = 0), "depth")
  expect_error(MTPart(features, outcomes, outcome_defs, mutate_example, depth = 1.5), "depth")
  expect_error(MTPart(features, outcomes, outcome_defs, mutate_example, nodesize = -5), "nodesize")
  expect_error(MTPart(features, outcomes, outcome_defs, mutate_example, cp = "a"), "cp")
  expect_error(MTPart(features, outcomes, outcome_defs, mutate_example, continuous = "bogus"), "continuous")
})

test_that("MTTest errors on a malformed tree object", {
  expect_error(
    MTTest(list(1, 2, 3), features, outcomes, outcome_defs, mutate_example),
    "does not look like an MTPart"
  )
  expect_error(
    MTTest(list(data.frame(x = 1), list()), features, outcomes, outcome_defs, mutate_example),
    "partitions table"
  )
})

test_that("MTTest reuses the same core input validation as MTPart", {
  tree <- MTPart(features, outcomes, outcome_defs, mutate_example, depth = 2, nodesize = 30)
  expect_error(
    MTTest(tree, c(features, "nonexistent_col"), outcomes, outcome_defs, mutate_example),
    "features.*not found"
  )
})

test_that("CV_Tune errors clearly on a bad Y or kfolds", {
  expect_error(
    CV_Tune(features, outcomes, outcome_defs, mutate_example, Y = "nonexistent_col"),
    "Y.*must be a single column name"
  )
  expect_error(
    CV_Tune(features, outcomes, outcome_defs, mutate_example, Y = "response", kfolds = 1),
    "kfolds"
  )
})
