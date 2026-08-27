outcomes <- c("response", "tumor_size", "ae_count", "OS_definition_time_status")
outcome_defs <- c("Cat", "Cont", "Count", "Surv")

test_that("MTSummary reports the expected fields per outcome type", {
  expected_fields <- list(
    Cat   = "expectedloss",
    Cont  = "MSE",
    Count = "deviance",
    Surv  = "deviance"
  )
  for (i in seq_along(outcomes)) {
    targets <- list(Definitions = outcome_defs[i],
                     Z = mutate_example[, outcomes[i], drop = FALSE])
    summ <- MTSummary(targets, mutate_example)
    expect_true(expected_fields[[outcome_defs[i]]] %in% names(summ[[1]]),
                info = paste("outcome type:", outcome_defs[i]))
    expect_equal(summ[[1]]$absent, 0)
  }
})

test_that("MTPartSummary produces a summary table with the documented columns", {
  features <- c("age", "sex", "biomarker")
  tree <- MTPart(features, outcomes, outcome_defs, mutate_example,
                  depth = 2, nodesize = 30)
  summ <- MTPartSummary(tree)
  expect_named(summ, c("node_data", "summary_table"))
  expect_named(summ$summary_table,
               c("nsplit", "leaves", "CP", "AvgRelError", "TotRelError", "Xerror", "Xstd", "Eval"))
  expect_equal(nrow(summ$summary_table), length(unique(summ$node_data$parent)))
})
