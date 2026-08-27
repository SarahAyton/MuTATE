features <- c("age", "sex", "biomarker")
outcomes <- c("response", "tumor_size", "ae_count", "OS_definition_time_status")
outcome_defs <- c("Cat", "Cont", "Count", "Surv")

test_that("CVsplitEval reports the expected metric fields per outcome type", {
  targets <- list(Definitions = outcome_defs, Z = mutate_example[, outcomes])
  res <- CVsplitEval("age < 60", mutate_example[, features], targets, mutate_example)

  expect_named(res, outcomes)
  expect_setequal(names(res$response), c("Accuracy", "AccuracySD", "Kappa", "KappaSD"))
  expect_setequal(names(res$tumor_size), c("RMSE", "RMSESD", "Rsquared", "RsquaredSD", "MAE", "MAESD"))
  expect_setequal(names(res$ae_count), c("RMSE", "RMSESD", "Rsquared", "RsquaredSD", "MAE", "MAESD"))
})

test_that("CVsplitEval accepts a categorical split condition", {
  targets <- list(Definitions = outcome_defs, Z = mutate_example[, outcomes])
  res <- CVsplitEval("sex in F", mutate_example[, features], targets, mutate_example)
  expect_named(res, outcomes)
})
