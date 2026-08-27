features <- c("age", "sex", "biomarker")
outcomes <- c("response", "tumor_size", "ae_count", "OS_definition_time_status")
outcome_defs <- c("Cat", "Cont", "Count", "Surv")

test_that("MTPrune with cp = 0 keeps the full tree", {
  tree <- MTPart(features, outcomes, outcome_defs, mutate_example,
                  depth = 3, nodesize = 20)
  n_full <- sum(!vapply(tree$tree_nodes, is.null, logical(1)))
  pruned <- MTPrune(tree, cp = 0)
  n_pruned <- sum(!vapply(pruned$tree_nodes, is.null, logical(1)))
  expect_equal(n_pruned, n_full)
})

test_that("MTPrune with cp = 1 prunes back to nothing but the root", {
  tree <- MTPart(features, outcomes, outcome_defs, mutate_example,
                  depth = 3, nodesize = 20)
  pruned <- MTPrune(tree, cp = 1)
  expect_equal(nrow(pruned$partitions), 0)
})
