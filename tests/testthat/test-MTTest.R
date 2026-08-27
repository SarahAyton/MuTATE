features <- c("age", "sex", "biomarker")
outcomes <- c("response", "tumor_size", "ae_count", "OS_definition_time_status")
outcome_defs <- c("Cat", "Cont", "Count", "Surv")

test_that("MTTest scores new data against a fitted tree without error", {
  tree <- MTPart(features, outcomes, outcome_defs, mutate_example,
                  depth = 2, nodesize = 30)
  test_result <- MTTest(tree, features, outcomes, outcome_defs, mutate_example)
  expect_type(test_result, "list")
  expect_length(test_result, 2)
  expect_equal(as.numeric(test_result[[2]][[1]]$N), nrow(mutate_example))
})
