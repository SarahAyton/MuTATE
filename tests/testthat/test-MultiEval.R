features <- c("age", "sex", "biomarker")
outcomes <- c("response", "tumor_size", "ae_count", "OS_definition_time_status")
outcome_defs <- c("Cat", "Cont", "Count", "Surv")

test_that("MultiEval returns an Eval matrix and a ranked ScaledInfGain table", {
  Ztype <- list(Definitions = outcome_defs, Z = mutate_example[, outcomes])
  res <- suppressWarnings(MultiEval(X = mutate_example[, features], Ztype = Ztype, data = mutate_example))

  expect_length(res, 2)
  expect_setequal(colnames(res[[1]]),
                   c("Zk", "Ztype", "Xm", "Xtype", "Split", "Eval", "SplitEntropy",
                     "Inf_Gain", "Inf_GainProp", "P_Value", "Rank", "Scale"))
  expect_true(all(c("Xm", "Split", "Rank") %in% colnames(res[[2]])))
  # Rank should be a complete permutation of 1..nrow with no ties collapsed to NA
  expect_false(any(is.na(res[[2]]$Rank)))
})

test_that("Cat_Eval/Cont_Eval/Count_Eval/Surv_Eval each return a well-formed evaluation matrix", {
  X <- mutate_example[, features]
  Z <- mutate_example[, outcomes]

  cat_res   <- Cat_Eval(X, 1, Z, 1, mutate_example, continuous = "quantile")
  cont_res  <- Cont_Eval(X, 1, Z, 2, mutate_example, continuous = "quantile")
  count_res <- Count_Eval(X, 1, Z, 3, mutate_example, continuous = "quantile")
  surv_res  <- suppressWarnings(Surv_Eval(X, 1, Z, 4, mutate_example, continuous = "quantile"))

  for (res in list(cat_res, cont_res, count_res, surv_res)) {
    expect_true(is.matrix(res))
    expect_true(all(c("Xm", "Split", "Inf_Gain", "Rank") %in% colnames(res)))
    expect_true(nrow(res) > 0)
  }
})
