#' MTTest: Test Multi-Target Tree
#'
#' A function to perform tree-based evaluation on test data using a trained tree model.
#'
#' @param tree A list containing the trained tree model.
#' @param features A character vector specifying the names of the features used in the tree model.
#' @param outcomes A character vector specifying the names of the outcomes used in the tree model.
#' @param outcome_defs A list containing definitions of the outcomes used in the tree model.
#' @param data_test A data frame containing the test data.
#' @param wt An optional numeric vector of weights for the test data.
#'
#' @return A list of evaluation results for the test data.

MTTest <- function(tree, features, outcomes, outcome_defs, data_test, wt=NULL) {
  Z_test <- list("Definitions"=outcome_defs, "Z"=data_test[,c(outcomes)])
  X_test <- data_test[,c(features)]
  # Import trained tree
  df <- tree[[1]]
  mtpart_splits <- (tree[[2]])
  # Generate test tree
  mtpart_splits_test <- mtpart_splits
  i <- 1
  split_df <- splitX <- splitZ <- vector(mode = "list", length = length(mtpart_splits_test))
  split_df[[i]] <- temp  <- data_test
  splitZ[[i]]   <- tempZ <- Z_test
  splitX[[i]]   <- tempX <- X_test
  for (i in 1:length(mtpart_splits_test)) {
    if (i==1 & nrow(df) > 1 & !is.null(mtpart_splits_test[[i]])) {
      root_eval <- suppressWarnings(MTSummary(targets=Z_test, data=data_test))
      mtpart_splits_test[[i]]$N       <- paste(nrow(data_test))
      mtpart_splits_test[[i]]$Pnode   <- paste(nrow(data_test)/nrow(data_test))
      mtpart_splits_test[[i]]$Targets <- root_eval
      temp      <- split_df[[i]]
      splitvar  <- mtpart_splits[[2*i]]$SplitVar
      var       <- (unlist(strsplit(splitvar, " ")))[1]
      varnum    <- which(colnames(temp)==var)
      notation  <- (unlist(strsplit(splitvar, " ")))[2]
      thresh    <- (unlist(strsplit(splitvar, " ")))[3]
      ifelse(notation == "in", threshold <- str_remove(thresh, " in "),
             ifelse(notation == "<", threshold <- str_remove(thresh, " < "), NA))
      # Subset data based on selected split
      splitL_Z <- splitR_Z <- splitZ[[i]]
      ifelse(notation == "in", splitL <- (temp[which(temp[,c(varnum)] == threshold), c(1:ncol(temp))]),
             ifelse(notation == "<", splitL  <- (temp[which(temp[,c(varnum)] < as.numeric(threshold)), c(1:ncol(temp))]), NA))
      splitL_Z$Z <- subset(splitL_Z$Z, row.names(splitL_Z$Z) %in% row.names(splitL))
      ifelse(notation == "in", splitR <- (temp[which(temp[,c(varnum)] != threshold), c(1:ncol(temp))]),
             ifelse(notation == "<", splitR <- (temp[which(temp[,c(varnum)] >= as.numeric(threshold)), c(1:ncol(temp))]), NA))
      splitR_Z$Z <- subset(splitR_Z$Z, row.names(splitR_Z$Z) %in% row.names(splitR))
      ifelse(notation == "in", splitvarR <- paste(var, threshold, sep = " not in "),
             ifelse(notation == "<", splitvarR <- paste(var, threshold, sep = " >= "), NA))
      if ((!is.null(nrow(splitL_Z[[2]])) & nrow(splitL_Z[[2]]) > 0) &
          (!is.null(nrow(splitR_Z[[2]])) & nrow(splitR_Z[[2]]) > 0)) {
        # Evaluate child nodes across targets using MTSummary
        splitL_eval <- suppressWarnings(MTSummary(splitL_Z, data=splitL))
        splitR_eval <- suppressWarnings(MTSummary(splitR_Z, data=splitR))
        # Evaluate split using CV error & sd (caret) --> Record for parent node
        misclass_targets_L <- misclass_targets_R <- as.list(1:length(splitZ[[i]]$Z))
        for (j in 1:length(splitL_eval)) {
          if(splitZ[[i]]$Definitions[j]=="Cat") {
            splitL_eval[[j]]$relerr <- splitL_eval[[j]]$expectedloss/mtpart_splits_test[[1]]$Targets[[j]]$expectedloss
          }
          else if(splitZ[[i]]$Definitions[j] == "Surv") {
            splitL_eval[[j]]$relerr <- splitL_eval[[j]]$deviance/mtpart_splits_test[[1]]$Targets[[j]]$deviance
          }
          else if(splitZ[[i]]$Definitions[j] == "Count") {
            splitL_eval[[j]]$relerr <- splitL_eval[[j]]$deviance/mtpart_splits_test[[1]]$Targets[[j]]$deviance
          }
          else {
            splitL_eval[[j]]$relerr <- splitL_eval[[j]]$MSE/mtpart_splits_test[[1]]$Targets[[j]]$MSE
          }
          misclass_targets_L[[j]] <- splitL_eval[[j]]$relerr
        }
        ifelse(!is.null(wt), (L_relerr <- weighted.mean(as.numeric(misclass_targets_L), w=(wt), na.rm=TRUE)),
               (L_relerr <- mean(as.numeric(misclass_targets_L), na.rm = TRUE)))
        for (j in 1:length(splitR_eval)) {
          if(splitZ[[i]]$Definitions[j]=="Cat") {
            splitR_eval[[j]]$relerr <- splitR_eval[[j]]$expectedloss/mtpart_splits_test[[1]]$Targets[[j]]$expectedloss
          }
          else if(splitZ[[i]]$Definitions[j] == "Surv") {
            splitR_eval[[j]]$relerr <- splitR_eval[[j]]$deviance/mtpart_splits_test[[1]]$Targets[[j]]$deviance
          }
          else if(splitZ[[i]]$Definitions[j] == "Count") {
            splitR_eval[[j]]$relerr <- splitR_eval[[j]]$deviance/mtpart_splits_test[[1]]$Targets[[j]]$deviance
          }
          else {
            splitR_eval[[j]]$relerr <- splitR_eval[[j]]$MSE/mtpart_splits_test[[1]]$Targets[[j]]$MSE
          }
          misclass_targets_R[[j]] <- splitR_eval[[j]]$relerr
        }
        ifelse(!is.null(wt), (R_relerr <- weighted.mean(as.numeric(misclass_targets_R), w=(wt), na.rm=TRUE)),
               (R_relerr <- mean(as.numeric(misclass_targets_R), na.rm = TRUE)))
        mtpart_splits_test[[i]]$CP      <- mtpart_splits_test[[i]]$relerr - ((L_relerr + R_relerr)/2)
        mtpart_splits_test[[i]]$CVeval <- CVsplitEval(mtpart_splits[[2*i]]$SplitVar, X_test, Z_test, data_test)
        CVeval_Xerr_targets <- CVeval_Xstd_targets <- as.list(1:length(Z_test$Z))
        for (j in 1:length(Z_test$Z)) {
          if(Z_test$Definitions[j]=="Cat") {
            mtpart_splits_test[[i]]$CVeval[[j]]$Xerror <- (1 - mtpart_splits_test[[i]]$CVeval[[j]]$Accuracy)
            mtpart_splits_test[[i]]$CVeval[[j]]$Xstd   <- mtpart_splits_test[[i]]$CVeval[[j]]$AccuracySD
          }
          else if(Z_test$Definitions[j] == "Surv") {
            mtpart_splits_test[[i]]$CVeval[[j]]$Xerror <- (1 - mtpart_splits_test[[i]]$CVeval[[j]]$Rsquared)
            mtpart_splits_test[[i]]$CVeval[[j]]$Xstd   <- mtpart_splits_test[[i]]$CVeval[[j]]$RsquaredSD
          }
          else if(Z_test$Definitions[j] == "Count") {
            mtpart_splits_test[[i]]$CVeval[[j]]$Xerror <- (1 - mtpart_splits_test[[i]]$CVeval[[j]]$Rsquared)
            mtpart_splits_test[[i]]$CVeval[[j]]$Xstd   <- mtpart_splits_test[[i]]$CVeval[[j]]$RsquaredSD
          }
          else {
            mtpart_splits_test[[i]]$CVeval[[j]]$Xerror <- (1 - mtpart_splits_test[[i]]$CVeval[[j]]$Rsquared)
            mtpart_splits_test[[i]]$CVeval[[j]]$Xstd   <- mtpart_splits_test[[i]]$CVeval[[j]]$RsquaredSD
          }
          CVeval_Xerr_targets[[j]] <- mtpart_splits_test[[i]]$CVeval[[j]]$Xerror
          CVeval_Xstd_targets[[j]] <- mtpart_splits_test[[i]]$CVeval[[j]]$Xstd
        }
        ifelse(!is.null(wt), (mtpart_splits_test[[i]]$SXerror <- (weighted.mean(as.numeric(CVeval_Xerr_targets), w=(wt), na.rm=TRUE))),
               (mtpart_splits_test[[i]]$SXerror <- (mean(as.numeric(CVeval_Xerr_targets)))))
        ifelse(!is.null(wt), (mtpart_splits_test[[i]]$SXstd <- ( weighted.mean(as.numeric(CVeval_Xstd_targets), w=(wt), na.rm=TRUE))),
               (mtpart_splits_test[[i]]$SXstd <- ( mean(as.numeric(CVeval_Xstd_targets)))))
        mtpart_splits_test[[i]]$Eval <- NA
        # Record child dataframes
        split_df[[2*i+0]] <- splitL
        split_df[[2*i+1]] <- splitR
        splitZ[[2*i+0]] <- splitL_Z
        splitZ[[2*i+1]] <- splitR_Z
        splitX[[2*i+0]] <- subset(splitX[[i]], row.names(splitX[[i]]) %in% row.names(splitL))
        splitX[[2*i+1]] <- subset(splitX[[i]], row.names(splitX[[i]]) %in% row.names(splitR))
        # Record left-child node
        mtpart_splits_test[[2*i+0]]$N       <- paste(nrow(split_df[[2*i+0]]))
        mtpart_splits_test[[2*i+0]]$Pnode   <- paste(nrow(split_df[[2*i+0]])/nrow(split_df[[1]]))
        mtpart_splits_test[[2*i+0]]$Targets <- splitL_eval
        mtpart_splits_test[[2*i+0]]$relerr  <- L_relerr
        # Record right-child node
        mtpart_splits_test[[2*i+1]]$N       <- paste(nrow(split_df[[2*i+1]]))
        mtpart_splits_test[[2*i+1]]$Pnode   <- paste(nrow(split_df[[2*i+1]])/nrow(split_df[[1]]))
        mtpart_splits_test[[2*i+1]]$Targets <- splitR_eval
        mtpart_splits_test[[2*i+1]]$relerr  <- R_relerr
      }
      else {
        # Record child dataframes
        split_df[[2*i+0]] <- splitL
        split_df[[2*i+1]] <- splitR
        splitZ[[2*i+0]] <- splitL_Z
        splitZ[[2*i+1]] <- splitR_Z
        splitX[[2*i+0]] <- subset(splitX[[i]], row.names(splitX[[i]]) %in% row.names(splitL))
        splitX[[2*i+1]] <- subset(splitX[[i]], row.names(splitX[[i]]) %in% row.names(splitR))
        misclass_targets_L <- misclass_targets_R <- as.list(1:length(splitZ[[i]]$Z))

      }
    }
    else if (nrow(df) == 1 & !is.null(mtpart_splits_test[[i]])) {
      root_eval <- suppressWarnings(MTSummary(targets=Z_test, data=data_test))
      mtpart_splits_test[[i]]$N       <- paste(nrow(data_test))
      mtpart_splits_test[[i]]$Pnode   <- paste(nrow(data_test)/nrow(data_test))
      mtpart_splits_test[[i]]$Targets <- root_eval
      temp   <- split_df[[i]]
    }
    else if (2*i>length(mtpart_splits_test) & !is.null(mtpart_splits_test[[i]])) {
      mtpart_splits_test[[i]]$Eval <- NA
    }
    else if (!is.null(mtpart_splits_test[[i]]) && !is.null(split_df[[i]]) && !is.null(mtpart_splits[[2*i]]) &&
             length(levels(droplevels(as.factor(split_df[[i]][,c(which(colnames(split_df[[i]])==(unlist(strsplit(mtpart_splits[[2*i]]$SplitVar, " ")))[1]))] ))))>=2 &&
             2*i<=length(mtpart_splits_test)) {
      temp      <- split_df[[i]]
      splitvar  <- mtpart_splits[[2*i]]$SplitVar
      var       <- (unlist(strsplit(splitvar, " ")))[1]
      varnum    <- which(colnames(temp)==var)
      notation  <- (unlist(strsplit(splitvar, " ")))[2]
      thresh    <- (unlist(strsplit(splitvar, " ")))[3]
      ifelse(notation == "in", threshold <- str_remove(thresh, " in "),
             ifelse(notation == "<", threshold <- str_remove(thresh, " < "), NA))
      # Subset data based on selected split
      splitL_Z <- splitR_Z <- splitZ[[i]]
      ifelse(notation == "in", splitL <- (temp[which(temp[,c(varnum)] == threshold), c(1:ncol(temp))]),
             ifelse(notation == "<", splitL  <- (temp[which(temp[,c(varnum)] < as.numeric(threshold)), c(1:ncol(temp))]), NA))
      splitL_Z$Z <- subset(splitL_Z$Z, row.names(splitL_Z$Z) %in% row.names(splitL))
      ifelse(notation == "in", splitR <- (temp[which(temp[,c(varnum)] != threshold), c(1:ncol(temp))]),
             ifelse(notation == "<", splitR <- (temp[which(temp[,c(varnum)] >= as.numeric(threshold)), c(1:ncol(temp))]), NA))
      splitR_Z$Z <- subset(splitR_Z$Z, row.names(splitR_Z$Z) %in% row.names(splitR))
      ifelse(notation == "in", splitvarR <- paste(var, threshold, sep = " not in "),
             ifelse(notation == "<", splitvarR <- paste(var, threshold, sep = " >= "), NA))
      if (nrow(splitL_Z$Z)==0 | nrow(splitR_Z$Z)==0 | length(levels(as.factor(temp[,c(varnum)])))==0) {mtpart_splits_test[[i]]$Eval <- mtpart_splits_test[[i-1]]$Eval}
      else {
        # Evaluate child nodes across targets using MTSummary
        splitL_eval <- try(suppressWarnings(MTSummary(splitL_Z, data=splitL)))
        splitR_eval <- try(suppressWarnings(MTSummary(splitR_Z, data=splitR)))
        # Evaluate split using CV error & sd (caret) --> Record for parent node
        misclass_targets_L <- misclass_targets_R <- as.list(1:length(splitZ[[i]]$Z))
        for (j in 1:length(splitL_eval)) {
          if(splitZ[[i]]$Definitions[j]=="Cat") {
            splitL_eval[[j]]$relerr <- splitL_eval[[j]]$expectedloss/mtpart_splits_test[[1]]$Targets[[j]]$expectedloss
          }
          else if(splitZ[[i]]$Definitions[j] == "Surv") {
            splitL_eval[[j]]$relerr <- splitL_eval[[j]]$deviance/mtpart_splits_test[[1]]$Targets[[j]]$deviance
          }
          else if(splitZ[[i]]$Definitions[j] == "Count") {
            splitL_eval[[j]]$relerr <- splitL_eval[[j]]$deviance/mtpart_splits_test[[1]]$Targets[[j]]$deviance
          }
          else {
            splitL_eval[[j]]$relerr <- splitL_eval[[j]]$MSE/mtpart_splits_test[[1]]$Targets[[j]]$MSE
          }
          misclass_targets_L[[j]] <- splitL_eval[[j]]$relerr
        }
        ifelse(!is.null(wt), (L_relerr <- weighted.mean(as.numeric(misclass_targets_L), w=(wt), na.rm=TRUE)),
               (L_relerr <- mean(as.numeric(misclass_targets_L), na.rm = TRUE)))
        for (j in 1:length(splitR_eval)) {
          if(splitZ[[i]]$Definitions[j]=="Cat") {
            splitR_eval[[j]]$relerr <- splitR_eval[[j]]$expectedloss/mtpart_splits_test[[1]]$Targets[[j]]$expectedloss
          }
          else if(splitZ[[i]]$Definitions[j] == "Surv") {
            splitR_eval[[j]]$relerr <- splitR_eval[[j]]$deviance/mtpart_splits_test[[1]]$Targets[[j]]$deviance
          }
          else if(splitZ[[i]]$Definitions[j] == "Count") {
            splitR_eval[[j]]$relerr <- splitR_eval[[j]]$deviance/mtpart_splits_test[[1]]$Targets[[j]]$deviance
          }
          else {
            splitR_eval[[j]]$relerr <- splitR_eval[[j]]$MSE/mtpart_splits_test[[1]]$Targets[[j]]$MSE
          }
          misclass_targets_R[[j]] <- splitR_eval[[j]]$relerr
        }
        ifelse(!is.null(wt), (R_relerr <- weighted.mean(as.numeric(misclass_targets_R), w=(wt), na.rm=TRUE)),
               (R_relerr <- mean(as.numeric(misclass_targets_R), na.rm = TRUE)))
        mtpart_splits_test[[i]]$CP      <- mtpart_splits_test[[i]]$relerr - ((L_relerr + R_relerr)/2)
        splitCVeval <- CVsplitEval(splitvar, splitX[[i]], splitZ[[i]], split_df[[i]])
        mtpart_splits_test[[i]]$CVeval <- splitCVeval
        CVeval_Xerr_targets <- CVeval_Xstd_targets <- as.list(1:length(splitZ[[i]]$Z))
        for (j in 1:length(splitZ[[i]]$Z)) {
          if(splitZ[[i]]$Definitions[j]=="Cat") {
            mtpart_splits_test[[i]]$CVeval[[j]]$Xerror <- (1 - mtpart_splits_test[[i]]$CVeval[[j]]$Accuracy)
            mtpart_splits_test[[i]]$CVeval[[j]]$Xstd   <- mtpart_splits_test[[i]]$CVeval[[j]]$AccuracySD
          }
          else if(splitZ[[i]]$Definitions[j] == "Surv") {
            mtpart_splits_test[[i]]$CVeval[[j]]$Xerror <- (1 - mtpart_splits_test[[i]]$CVeval[[j]]$Rsquared)
            mtpart_splits_test[[i]]$CVeval[[j]]$Xstd   <- mtpart_splits_test[[i]]$CVeval[[j]]$RsquaredSD
          }
          else if(splitZ[[i]]$Definitions[j] == "Count") {
            mtpart_splits_test[[i]]$CVeval[[j]]$Xerror <- (1 - mtpart_splits_test[[i]]$CVeval[[j]]$Rsquared)
            mtpart_splits_test[[i]]$CVeval[[j]]$Xstd   <- mtpart_splits_test[[i]]$CVeval[[j]]$RsquaredSD
          }
          else {
            mtpart_splits_test[[i]]$CVeval[[j]]$Xerror <- (1 - mtpart_splits_test[[i]]$CVeval[[j]]$Rsquared)
            mtpart_splits_test[[i]]$CVeval[[j]]$Xstd   <- mtpart_splits_test[[i]]$CVeval[[j]]$RsquaredSD
          }
          CVeval_Xerr_targets[[j]] <- mtpart_splits_test[[i]]$CVeval[[j]]$Xerror
          CVeval_Xstd_targets[[j]] <- mtpart_splits_test[[i]]$CVeval[[j]]$Xstd
        }
        ifelse(!is.null(wt), (mtpart_splits_test[[i]]$SXerror <- (weighted.mean(as.numeric(CVeval_Xerr_targets), w=(wt), na.rm=TRUE))),
               (mtpart_splits_test[[i]]$SXerror <- (mean(as.numeric(CVeval_Xerr_targets)))))
        ifelse(!is.null(wt), (mtpart_splits_test[[i]]$SXstd <- ( weighted.mean(as.numeric(CVeval_Xstd_targets), w=(wt), na.rm=TRUE))),
               (mtpart_splits_test[[i]]$SXstd <- ( mean(as.numeric(CVeval_Xstd_targets)))))
        mtpart_splits_test[[i]]$Eval <- NA
        # Record child dataframes
        split_df[[2*i+0]] <- splitL
        split_df[[2*i+1]] <- splitR
        splitZ[[2*i+0]] <- splitL_Z
        splitZ[[2*i+1]] <- splitR_Z
        splitX[[2*i+0]] <- subset(splitX[[i]], row.names(splitX[[i]]) %in% row.names(splitL))
        splitX[[2*i+1]] <- subset(splitX[[i]], row.names(splitX[[i]]) %in% row.names(splitR))
        # Record left-child node
        mtpart_splits_test[[2*i+0]]$N       <- paste(nrow(split_df[[2*i+0]]))
        mtpart_splits_test[[2*i+0]]$Pnode   <- paste(nrow(split_df[[2*i+0]])/nrow(split_df[[1]]))
        mtpart_splits_test[[2*i+0]]$Targets <- splitL_eval
        mtpart_splits_test[[2*i+0]]$relerr  <- L_relerr
        # Record right-child node
        mtpart_splits_test[[2*i+1]]$N       <- paste(nrow(split_df[[2*i+1]]))
        mtpart_splits_test[[2*i+1]]$Pnode   <- paste(nrow(split_df[[2*i+1]])/nrow(split_df[[1]]))
        mtpart_splits_test[[2*i+1]]$Targets <- splitR_eval
        mtpart_splits_test[[2*i+1]]$relerr  <- R_relerr
      }
    }
    else if (!is.null(mtpart_splits_test[[i]]) & !is.null(split_df[[i]]) & 2*i<=length(mtpart_splits_test) & 2*i<length(mtpart_splits_test) &&
             is.null(mtpart_splits[[2*i]])) {
      mtpart_splits_test[[i]]$Eval <- NA
    }
    else { i <- i + 1 }
  }
  # Return test tree
  return(list(df, mtpart_splits_test))
}
