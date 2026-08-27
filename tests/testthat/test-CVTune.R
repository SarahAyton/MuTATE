test_that("CV_Tune runs a minimal single-point grid without error", {
  skip_on_cran()
  features <- c("age", "sex", "biomarker")
  outcomes <- c("response", "tumor_size", "ae_count", "OS_definition_time_status")
  outcome_defs <- c("Cat", "Cont", "Count", "Surv")

  cv_results <- CV_Tune(features, outcomes, outcome_defs, mutate_example,
                         kfolds = 2, Y = "response",
                         drange = 2, noderange = 30, splitmin_div = 2,
                         method = "avgIG", alpharange = 0.05, igrange = 0.95,
                         psplitrange = 1, pdepthrange = 1, cp_val = 0)

  expect_length(cv_results, 2)
  expect_s3_class(cv_results[[1]], "data.frame")
  # one row per evaluation metric, all from the single grid point requested
  expect_equal(length(unique(cv_results[[1]]$depth)), 1)
  expect_equal(length(unique(cv_results[[1]]$nodesize)), 1)
})
