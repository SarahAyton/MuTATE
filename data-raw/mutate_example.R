## Generates the small synthetic dataset shipped with the package for use in
## function examples and tests. Not a realistic clinical dataset - purely
## illustrative of the four supported outcome types (Cat, Cont, Count, Surv).

set.seed(42)
n <- 200
mutate_example <- data.frame(
  age        = round(rnorm(n, 60, 10)),
  sex        = factor(sample(c("F", "M"), n, replace = TRUE)),
  biomarker  = rnorm(n, 5, 2),
  response   = factor(sample(c("responder", "non_responder"), n, replace = TRUE)),
  tumor_size = rnorm(n, 30, 8),
  ae_count   = rpois(n, 2),
  time       = rexp(n, 0.05),
  status     = rbinom(n, 1, 0.6)
)
# MTPart's Surv outcome convention: the outcome column name itself must
# encode the time/event column names as its 3rd/4th "_"-separated tokens
# (see MTSummary.R), and must also exist as an actual data column.
mutate_example$OS_definition_time_status <- mutate_example$status

usethis::use_data(mutate_example, overwrite = TRUE)
