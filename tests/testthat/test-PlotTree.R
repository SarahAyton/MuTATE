test_that("PlotTree runs without error on a fitted tree", {
  features <- c("age", "sex", "biomarker")
  outcomes <- c("response", "tumor_size", "ae_count", "OS_definition_time_status")
  outcome_defs <- c("Cat", "Cont", "Count", "Surv")

  tree <- MTPart(features, outcomes, outcome_defs, mutate_example,
                  depth = 2, nodesize = 30)

  plot_file <- tempfile(fileext = ".png")
  grDevices::png(plot_file)
  on.exit({
    grDevices::dev.off()
    unlink(plot_file)
  })
  expect_no_error(PlotTree(tree))
})
