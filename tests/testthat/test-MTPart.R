features <- c("age", "sex", "biomarker")
outcomes <- c("response", "tumor_size", "ae_count", "OS_definition_time_status")
outcome_defs <- c("Cat", "Cont", "Count", "Surv")

test_that("MTPart fits on synthetic multi-target data without error", {
  tree <- MTPart(features, outcomes, outcome_defs, mutate_example,
                  depth = 2, nodesize = 30)
  expect_type(tree, "list")
  expect_named(tree, c("partitions", "tree_nodes"))
  expect_s3_class(tree$partitions, "data.frame")
  expect_true(is.list(tree$tree_nodes))
})

test_that("MTPart root node summarizes the full dataset", {
  tree <- MTPart(features, outcomes, outcome_defs, mutate_example,
                  depth = 2, nodesize = 30)
  root <- tree$tree_nodes[[1]]
  expect_equal(as.numeric(root$N), nrow(mutate_example))
  expect_equal(root$SplitVar, "Root")
})

test_that("MTPart does not leak wt/reuse/cp/nodesize into the calling environment", {
  # Regression test: MTPart previously used `<<-` to write these into the
  # global environment on every call, corrupting state across calls/users.
  rm(list = intersect(c("wt", "reuse", "cp", "nodesize"), ls(envir = .GlobalEnv)),
     envir = .GlobalEnv)
  invisible(MTPart(features, outcomes, outcome_defs, mutate_example,
                    depth = 2, nodesize = 30, cp = 0.02, wt = NULL, reuse = FALSE))
  expect_false(any(c("wt", "reuse", "cp", "nodesize") %in% ls(envir = .GlobalEnv)))
})

test_that("MTPart runs with a single outcome of each supported type", {
  for (i in seq_along(outcomes)) {
    tree <- MTPart(features, outcomes[i], outcome_defs[i], mutate_example,
                    depth = 2, nodesize = 30)
    expect_true(is.list(tree$tree_nodes[[1]]$Targets),
                info = paste("outcome type:", outcome_defs[i]))
  }
})

test_that("MTPart respects the nodesize floor (no split below nodesize)", {
  # nodesize larger than the whole dataset should prevent any split:
  # the root is marked terminal ("1 *") and no child nodes are created.
  tree <- MTPart(features, outcomes, outcome_defs, mutate_example,
                  depth = 3, nodesize = 1000)
  expect_equal(nrow(tree$partitions), 1)
  expect_match(tree$tree_nodes[[1]]$NodeID, "\\*")
  expect_null(tree$tree_nodes[[2]])
})

test_that("MTPart runs with reuse = TRUE (allow reusing the parent split variable)", {
  tree <- MTPart(features, outcomes, outcome_defs, mutate_example,
                  depth = 2, nodesize = 30, reuse = TRUE)
  expect_s3_class(tree$partitions, "data.frame")
})

test_that("MTPart runs with explicit outcome weights", {
  tree <- MTPart(features, outcomes, outcome_defs, mutate_example,
                  depth = 2, nodesize = 30, wt = c(1, 1, 1, 2))
  expect_s3_class(tree$partitions, "data.frame")
})

test_that("MTPart runs with parallelpart = TRUE (lookahead split search)", {
  tree <- MTPart(features, outcomes, outcome_defs, mutate_example,
                  depth = 2, nodesize = 30,
                  parallelpart = TRUE, parallelsplit = 2, paralleldepth = 1)
  expect_s3_class(tree$partitions, "data.frame")
  expect_true(nrow(tree$partitions) >= 1)
})

test_that("MTPart runs with parallelpart = TRUE and parallelsplit = \"all\"", {
  tree <- MTPart(features, outcomes, outcome_defs, mutate_example,
                  depth = 2, nodesize = 30,
                  parallelpart = TRUE, parallelsplit = "all", paralleldepth = 1)
  expect_s3_class(tree$partitions, "data.frame")
})

test_that("MTPart runs with parallelpart = TRUE at a deeper paralleldepth", {
  tree <- MTPart(features, outcomes, outcome_defs, mutate_example,
                  depth = 3, nodesize = 20,
                  parallelpart = TRUE, parallelsplit = 2, paralleldepth = 2)
  expect_s3_class(tree$partitions, "data.frame")
  expect_true(nrow(tree$partitions) >= 1)
})

test_that("MTPart runs with depth = 1 (root split only)", {
  tree <- MTPart(features, outcomes, outcome_defs, mutate_example,
                  depth = 1, nodesize = 30)
  expect_true(length(tree$tree_nodes) <= 2^1)
})
