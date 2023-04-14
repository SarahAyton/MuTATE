#' CV_Tune: Cross-Validated Hyperparameter Tuning for MTPart Models
#'
#' This function performs hyperparameter tuning for MTPart models using k-fold cross-validation.
#'
#' @param features A character vector of feature names.
#' @param outcomes A character vector of outcome names.
#' @param outcome_defs A list of outcome definitions.
#' @param data A data frame containing the input data.
#' @param continuous A string indicating the method for handling continuous outcomes. Default is "quantile".
#' @param quantseq A numeric vector specifying the quantile sequence for quantile-based methods. Default is seq(0,1,0.25).
#' @param wt A numeric vector of weights for the samples. Default is NULL.
#' @param reuse A logical value indicating whether to reuse previous MTPart model fits. Default is FALSE.
#' @param kfolds An integer specifying the number of folds for k-fold cross-validation. Default is 10.
#' @param Y A character string specifying the outcome variable for partitioning.
#' @param seedno An integer specifying the seed for random number generation. Default is 1234.
#' @param drange A numeric vector specifying the range of values for the 'depth' parameter. Default is seq(2, 6, by=2).
#' @param noderange A numeric vector specifying the range of values for the 'nodesize' parameter. Default is seq(10, 30, by=10).
#' @param splitmin_div A numeric vector specifying the range of values for the 'splitmin_div' parameter. Default is seq(2, 4, by=1).
#' @param method A character vector specifying the method for selecting the best split. Default is c("avgIG","maxIG","mostIG","avgPVal","minPVal","mostPVal","splitError").
#' @param alpharange A numeric vector specifying the range of values for the 'alpha' parameter. Default is seq(0.05, 0.15, by=0.05).
#' @param igrange A numeric vector specifying the range of values for the 'IGcutoff' parameter. Default is seq(0.85, 0.95, by=0.05).
#' @param psplitrange A numeric vector specifying the range of values for the 'psplit' parameter. Default is seq(1, 2, by=1).
#' @param pdepthrange A numeric vector specifying the range of values for the 'pdepth' parameter. Default is seq(1, 2, by=1).
#' @param cp_val A numeric vector specifying the range of values for the 'CP' parameter. Default is c(seq(-1, 0, by=0.5)).
#'
#' @usage CV_Tune(features, outcomes, outcome_defs, data,
#'         continuous, quantseq, wt, reuse, kfolds, Y, seedno,
#'         drange, noderange, splitmin_div, method, alpharange,
#'         igrange, psplitrange, pdepthrange, cp_val)
#'
#' @return A data frame containing the results of hyperparameter tuning for the MTPart model.
#' @export
#' @import caret
#' @exclude caret precision
#' @importFrom caret::cluster exclude=TRUE
#' @importFrom caret::recall exclude=TRUE

CV_Tune <- function(features, outcomes, outcome_defs, data,
                    continuous = "quantile", quantseq = seq(0,1,0.25), wt=NULL, reuse=FALSE,
                    kfolds=10, Y, seedno = 1234,  #Y = outcome variable for partitioning
                    drange    = c(seq(2, 6, by=2)),
                    noderange = c(seq(10, 30, by=10)),
                    splitmin_div  = c(seq(2, 4, by=1)),
                    method = c("avgIG","maxIG","mostIG","avgPVal","minPVal","mostPVal","splitError"),
                    alpharange  = c(seq(0.05, 0.15, by=0.05)),
                    igrange  = c(seq(0.85, 0.95, by=0.05)),
                    psplitrange = c(seq(1, 2, by=1)),
                    pdepthrange = c(seq(1, 2, by=1)),
                    cp_val = c(seq(-1, 0, by=0.5)) ) {
  set.seed(seedno)
  # Create a grid of parameter combinations to search over
  paramGrid <- expand.grid(depth=drange, nodesize=noderange, splitmin_div=splitmin_div, method=method,
                           alpha=alpharange, IGcutoff=igrange, psplit=psplitrange, pdepth=pdepthrange,
                           CP=cp_val)
  paramGrid$ppart <- ifelse(paramGrid$method=="splitError", TRUE, FALSE)
  # Initialize table to store results
  results <- results2 <- data.frame()

  # Define function for fitting and evaluating a single model
  fit_and_eval <- function(test_index, params, data) {
    # Split data into training and testing sets
    train_data <- data[-test_index, ]
    test_data <- data[test_index, ]
    # Fit model on training data
    starttim <- Sys.time()
    model <- MTPart(features, outcomes, outcome_defs, train_data, continuous=continuous,
                    quantseq=quantseq, wt=wt, evalmethod=params$method, alpha=params$alpha,
                    IGcutoff=params$IGcutoff, depth=params$depth, nodesize=params$nodesize,
                    cp=params$CP, reuse=reuse, parallelpart=params$ppart, parallelsplit=params$psplit,
                    paralleldepth=params$pdepth, splitmin=floor(params$nodesize/params$splitmin_div) )
    endtim <- Sys.time()
    # Evaluate model on training data
    summarytable_train <- MTPartSummary(model)
    eval_table <- summarytable_train$summary_table
    eval_table$N_train <- nrow(train_data)
    eval_table$N_test <- nrow(test_data)
    eval_table$RunTime <- as.numeric(difftime(endtim, starttim), units="secs")
    # Test model on testing data
    starttim <- Sys.time()
    model_test <- MTTest(model, features, outcomes, outcome_defs, test_data)
    endtim <- Sys.time()
    eval_table$RunTime_test <- as.numeric(difftime(endtim, starttim), units="secs")
    # Evaluate model on testing data
    evaluations_train <- eval_table[(nrow(eval_table)),]
    colnames(evaluations_train) <- c("nsplit", "leaves", "CP_train", "AvgRelError_train", "TotRelError_train",
                                     "Xerror_train", "Xstd_train", "Eval_train", "N_train", "N_test", "RunTime_train", "RunTime_test")
    summarytable_test <- MTPartSummary(model_test)
    evaluations_test <- summarytable_test$summary_table[(nrow(summarytable_test$summary_table)),]
    colnames(evaluations_test) <- c("nsplit", "leaves", "CP_test", "AvgRelError_test", "TotRelError_test",
                                    "Xerror_test", "Xstd_test", "Eval_test")
    evaluations <- merge(evaluations_train, evaluations_test)
    # Return accuracy and parameter values
    return(evaluations)
  }

  # Create a list of indices for k-fold cross-validation
  folds <- caret::createFolds(data[,Y], k=kfolds, list = TRUE, returnTrain = FALSE)
  # Split data into k-folds
  fold_dfs <- lapply(seq_len(kfolds), function(k) {
    # Get fold indices
    fold_indices <- folds[[k]]
    # Split data into training and testing sets
    train_indices <- unname(unlist(folds[!(folds %in% folds[[k]])]))
    test_indices  <- folds[[k]]
    # Return a list of training and testing data frames
    list(train = train_indices, test = test_indices)
  })

  # Perform grid search with k-fold cross-validation
  for (i in 1:nrow(paramGrid)) {
    # Fit and evaluate model for current parameter combination across all folds
    results_i <- sapply(seq_along(folds), function(f) {
      fit_and_eval(fold_dfs[[f]]$test, paramGrid[i,], data)
    })

    # Convert all values in results_i to numeric
    results_i_2 <- apply(results_i, 2, function(x) as.numeric(as.character(x)))
    varnames <- rownames(results_i)

    # Take the mean of all rows in results_i_2
    result_average <- apply(results_i_2, MARGIN = 1, FUN = function(x) mean(x, na.rm=TRUE))

    # Add parameter values and accuracy to results
    result_row <- data.frame(result_average,
                             depth=paramGrid[i,]$depth, nodesize=paramGrid[i,]$nodesize,
                             splitmin_div=paramGrid[i,]$splitmin_div, method=paramGrid[i,]$method,
                             alpha=paramGrid[i,]$alpha, IGcutoff=paramGrid[i,]$IGcutoff,
                             psplit=paramGrid[i,]$psplit, pdepth=paramGrid[i,]$pdepth, CP=paramGrid[i,]$CP)
    result_row$EvalName <- varnames
    results <- rbind(results, result_row)

    result_row2 <- data.frame(results_i,
                              depth=paramGrid[i,]$depth, nodesize=paramGrid[i,]$nodesize,
                              splitmin_div=paramGrid[i,]$splitmin_div, method=paramGrid[i,]$method,
                              alpha=paramGrid[i,]$alpha, IGcutoff=paramGrid[i,]$IGcutoff,
                              psplit=paramGrid[i,]$psplit, pdepth=paramGrid[i,]$pdepth, CP=paramGrid[i,]$CP)
    result_row2$EvalName <- varnames
    results2 <- rbind(results2, result_row2)

  }

  # Return results
  return(list(results, results2))
}
