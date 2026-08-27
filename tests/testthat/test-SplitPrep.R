test_that("SplitPrep drops single-level factors and constant continuous features", {
  df <- mutate_example
  df$single_level <- factor("A", levels = c("A", "B"))
  df$constant_num <- 5
  Xdf <- df[, c("age", "sex", "biomarker", "single_level", "constant_num")]

  prep <- SplitPrep(Xdf = Xdf, df = df, data_splt = df, parentsplit = "skip")
  expect_setequal(colnames(prep[[1]]), c("age", "sex", "biomarker"))
})

test_that("SplitPrep drops the parent split variable", {
  df <- mutate_example
  Xdf <- df[, c("age", "sex", "biomarker")]

  prep <- SplitPrep(Xdf = Xdf, df = df, data_splt = df, parentsplit = "age")
  expect_setequal(colnames(prep[[1]]), c("sex", "biomarker"))
})

test_that("SplitPrep drops categorical features with a level below 5% prevalence", {
  df <- mutate_example
  df$rare_level <- factor(ifelse(seq_len(nrow(df)) <= 3, "rare", "common"))
  Xdf <- df[, c("age", "sex", "biomarker", "rare_level")]

  prep <- SplitPrep(Xdf = Xdf, df = df, data_splt = df, parentsplit = "skip")
  expect_false("rare_level" %in% colnames(prep[[1]]))
  expect_setequal(colnames(prep[[1]]), c("age", "sex", "biomarker"))
})

test_that("SplitPrep keeps all features when none are degenerate", {
  df <- mutate_example
  Xdf <- df[, c("age", "sex", "biomarker")]
  prep <- SplitPrep(Xdf = Xdf, df = df, data_splt = df, parentsplit = "skip")
  expect_setequal(colnames(prep[[1]]), c("age", "sex", "biomarker"))
})
