#' Validate the core features/outcomes/outcome_defs/data/wt inputs shared by
#' MTPart, MTTest, and CV_Tune.
#'
#' Not exported - internal helper only. Raises an informative error via
#' \code{stop()} instead of letting malformed input fail deep inside the
#' recursive partitioning loop with a cryptic subscript/NA error.
#'
#' @noRd
validate_core_inputs <- function(features, outcomes, outcome_defs, data, wt = NULL,
                                  data_arg = "data") {
  if (!is.data.frame(data)) {
    stop("`", data_arg, "` must be a data frame, not ", class(data)[1], ".", call. = FALSE)
  }
  if (nrow(data) == 0) {
    stop("`", data_arg, "` has zero rows.", call. = FALSE)
  }
  if (!is.character(features) || length(features) == 0) {
    stop("`features` must be a non-empty character vector of column names.", call. = FALSE)
  }
  if (!is.character(outcomes) || length(outcomes) == 0) {
    stop("`outcomes` must be a non-empty character vector of column names.", call. = FALSE)
  }
  missing_features <- setdiff(features, colnames(data))
  if (length(missing_features) > 0) {
    stop("`features` not found in `", data_arg, "`: ",
         paste(missing_features, collapse = ", "), call. = FALSE)
  }
  missing_outcomes <- setdiff(outcomes, colnames(data))
  if (length(missing_outcomes) > 0) {
    stop("`outcomes` not found in `", data_arg, "`: ",
         paste(missing_outcomes, collapse = ", "), call. = FALSE)
  }
  overlap <- intersect(features, outcomes)
  if (length(overlap) > 0) {
    stop("The following columns appear in both `features` and `outcomes`: ",
         paste(overlap, collapse = ", "), call. = FALSE)
  }
  if (length(outcome_defs) != length(outcomes)) {
    stop("`outcome_defs` must have the same length as `outcomes` (",
         length(outcomes), "), got length ", length(outcome_defs), ".", call. = FALSE)
  }
  valid_defs <- c("Cat", "Cont", "Count", "Surv")
  bad_defs <- setdiff(unique(outcome_defs), valid_defs)
  if (length(bad_defs) > 0) {
    stop("`outcome_defs` must only contain ", paste(dQuote(valid_defs, q = FALSE), collapse = ", "),
         ". Invalid value(s): ", paste(bad_defs, collapse = ", "), call. = FALSE)
  }
  if (!is.null(wt)) {
    if (!is.numeric(wt) || length(wt) != length(outcomes)) {
      stop("`wt` must be a numeric vector with one weight per outcome (length ",
           length(outcomes), "), got length ", length(wt), ".", call. = FALSE)
    }
  }
  invisible(TRUE)
}

#' Validate MTPart-specific tuning parameters.
#' @noRd
validate_mtpart_params <- function(depth, nodesize, cp, continuous) {
  if (!is.numeric(depth) || length(depth) != 1 || depth < 1 || depth != round(depth)) {
    stop("`depth` must be a single positive integer, got ", depth, ".", call. = FALSE)
  }
  if (!is.numeric(nodesize) || length(nodesize) != 1 || nodesize < 1) {
    stop("`nodesize` must be a single positive number, got ", nodesize, ".", call. = FALSE)
  }
  if (!is.numeric(cp) || length(cp) != 1) {
    stop("`cp` must be a single numeric value, got ", cp, ".", call. = FALSE)
  }
  if (!continuous %in% c("all", "quantile")) {
    stop('`continuous` must be "all" or "quantile", got "', continuous, '".', call. = FALSE)
  }
  invisible(TRUE)
}

#' Validate that a `tree` argument looks like an MTPart() result.
#' @noRd
validate_tree_object <- function(tree) {
  if (!is.list(tree) || length(tree) != 2) {
    stop("`tree` does not look like an MTPart() result: expected a list of length 2 ",
         "(partitions, tree_nodes), got a ", class(tree)[1], " of length ", length(tree), ".",
         call. = FALSE)
  }
  if (!is.data.frame(tree[[1]]) || !all(c("parent", "child") %in% colnames(tree[[1]]))) {
    stop("`tree`'s first element does not look like an MTPart() partitions table ",
         "(expected columns `parent` and `child`).", call. = FALSE)
  }
  if (!is.list(tree[[2]])) {
    stop("`tree`'s second element does not look like an MTPart() tree_nodes list.", call. = FALSE)
  }
  invisible(TRUE)
}
