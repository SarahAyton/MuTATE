#' @name MTPart
#' @title Construct multi-target decision trees
#'
#' @description Evaluates multiple outcome variables of different types for a set of features.
#'
#' @param features Vector containing feature names to be evaluated.
#' @param outcomes Data frame object of all outcomes being evaluated. Used only in \code{Cat_Eval,
#' Cont_Eval, Count_Eval, Surv_Eval}.
#' @param outcome_defs Vector containing a list of outcome types ("Cat", "Cont", "Count", "Surv").
#' @param data The full dataset to be evaluated.
#' @param continuous How to handle partition evaluation of continuous features, in which
#' all potential split candidates may be examined ("all") or split candidates are derived
#' from quantiles ("quantile").
#' @param quantseq A user specified sequence ranging from 0 to 1 from which quantiles may
#' be specified (ex: seq(0,1,0.05)). If \code{continuous} and \code{quantseq} are not
#' specified the function will default to the "all" setting for continuous variables with
#' fewer than 100 unique values. If there are more than 100 unique values, the function
#' will use the "quantile" option with the default setting of \code{seq(0,1,0.05)}.
#' @param wt A numeric vector of the outcome variable weights to be applied in multi-target evaluations.
#' Default weights are set to \code{NULL}, indicating equal target variable importance.
#' @param evalmethod The method used to evaluate standardized proportion of information
#' gain across targets and target types.
#' @param alpha A numeric scalar specifying the significance level to be used in \code{evalmethod} if a p-value based
#' method (such as ) is selected. The default value for \code{alpha} is set to \code{0.05}.
#' @param IGcutoff A numeric scalar specifying the IG cutoff threshold to be used for \code{evalmethod} if methods
#' xxx or yyyy are selected. The default value for \code{IGcutoff} is \code{0.95}.
#' @param depth An integer scalar specifying the maximum depth of the tree (default is \code{6})
#' @param nodesize An integer scalar specifying the minimum number of samples per node (default is \code{20})
#' @param cp A numeric scalar specifying the complexity parameter for pruning (default is \code{0.02})
#' @param reuse A logical scalar specifying whether to reuse split variables at higher levels of the tree (default is \code{FALSE})
#' @param parallelpart A logical scalar specifying whether to use parallel computing for the partitioning step (default is \code{TRUE})
#' @param parallelsplit An integer scalar specifying the number of cores to use for parallel computing (default is \code{2})
#' @param paralleldepth An integer scalar specifying the minimum depth at which to start using parallel computing (default is \code{2})
#' @param splitmin An integer scalar specifying the minimum sample size \code{n} required for a node to be partitioned.
#' This is passed internally from \code{MTPart} and is used to assess feasibility of
#' feature binarization for evaluation (i.e., a binarized feature must result in
#' partitions with sufficient observations in both child nodes). The default value is of
#' \code{splitmin} is set to half of the \code{nodesize}.


MTPart <- function(features, outcomes, outcome_defs, data, continuous = "quantile", quantseq = seq(0,1,0.25),
                   wt=NULL, evalmethod = "avgIG", alpha = 0.05, IGcutoff = 0.95, depth=4,
                   nodesize=20, cp=-1.0, reuse=FALSE, parallelpart=FALSE, parallelsplit=NA,
                   paralleldepth=NA, splitmin=floor(nodesize/2)) { #and all inputs for MultiEval
  Z <- list("Definitions"=outcome_defs, "Z"=data[,c(outcomes)])
  X <- data[,c(features)]
  wt <<- wt
  reuse <<- reuse
  cp <<- cp
  nodesize <<- as.numeric(nodesize)
  mtpart_splits <- split_df <- splitX <- splitZ <- vector(mode = "list", length = 2^(depth))
  i <- 1
  split_df[[i]] <- temp  <- data
  splitZ[[i]]   <- tempZ <- Z
  splitX[[i]]   <- tempX <- X
  root_eval     <- suppressWarnings(MTSummary(targets=splitZ[[i]], data=split_df[[i]]))
  df <- data.frame(c(NA),c(NA))
  names(df) <- c("parent","child")
  mtpart_splits[[i]]$NodeID  <- paste(i)
  mtpart_splits[[i]]$SplitVar<- paste("Root")
  mtpart_splits[[i]]$Var     <- paste("Root")
  mtpart_splits[[i]]$Thresh  <- paste("Root")
  mtpart_splits[[i]]$N       <- paste(nrow(split_df[[i]]))
  mtpart_splits[[i]]$Pnode   <- paste(nrow(split_df[[i]])/nrow(split_df[[1]]))
  mtpart_splits[[i]]$Targets <- root_eval
  misclass_targets <- Xerr_targets <- Xstd_targets <- as.list(1:length(splitZ[[i]]$Z))
  for (j in 1:length(misclass_targets)) {
    if(splitZ[[i]]$Definitions[j]=="Cat") {
      mtpart_splits[[i]]$Targets[[j]]$relerr <- mtpart_splits[[i]]$Targets[[j]]$expectedloss/mtpart_splits[[1]]$Targets[[j]]$expectedloss
    }
    else if(splitZ[[i]]$Definitions[j] == "Surv") {
      mtpart_splits[[i]]$Targets[[j]]$relerr <- mtpart_splits[[i]]$Targets[[j]]$deviance/mtpart_splits[[1]]$Targets[[j]]$deviance
    }
    else if(splitZ[[i]]$Definitions[j] == "Count") {
      mtpart_splits[[i]]$Targets[[j]]$relerr <- mtpart_splits[[i]]$Targets[[j]]$deviance/mtpart_splits[[1]]$Targets[[j]]$deviance
    }
    else {
      mtpart_splits[[i]]$Targets[[j]]$relerr <- mtpart_splits[[i]]$Targets[[j]]$MSE/mtpart_splits[[1]]$Targets[[j]]$MSE
    }
    misclass_targets[[j]] <- mtpart_splits[[i]]$Targets[[j]]$relerr
  }
  ifelse(!is.null(wt), (mtpart_splits[[i]]$relerr <- weighted.mean(as.numeric(misclass_targets), w=(wt), na.rm=TRUE)),
         (mtpart_splits[[i]]$relerr <- mean(as.numeric(misclass_targets))))
  if (parallelpart==TRUE & !is.na(parallelsplit) & !is.na(paralleldepth)) {
    for (i in 1:2^(depth)) {
      # If the node i exists:
      if (!is.null(mtpart_splits[[i]]) ) {
        # And if the node i is sufficiently large (defined by nodesize argument):
        if (!is.null(data.frame(splitX[[i]])) & nrow(data.frame(splitX[[i]])) >= nodesize & ncol(data.frame(splitX[[i]])) > 0 ) {
          # Construct parallel trees of specified depth & select split var based on subtree relative error
          levlist <- vector(mode = "list", length = 0)
          for (o in 1:length(splitX[[i]])) {
            levlist <- append(levlist, levels(splitX[[i]][[o]]))
          }
          ifelse(parallelsplit=="all", parallelsplit <- length(levlist), parallelsplit <- parallelsplit)
          paralleltrees <- pareval <- vector(mode = "list", length = (parallelsplit))
          for (p in 1:length(paralleltrees)) {
            #Adjust to stop lookahead beyond tree depth?
            subtree <- subtree_df <- subtree_X <- subtree_Z <- vector(mode = "list", length = (2^(paralleldepth+1)-1))
            subtree[[1]] <- mtpart_splits[[i]]
            subtree_df[[1]] <- split_df[[i]]
            subtree_X[[1]]  <- splitX[[i]]
            subtree_Z[[1]]  <- splitZ[[i]]
            for (l in 1:length(subtree)) {
              # If the node i exists:
              if (!is.null(subtree[[l]]) ) {
                # And if the node i is sufficiently large (defined by nodesize argument):
                if (!is.null(data.frame(subtree_X[[l]])) & nrow(data.frame(subtree_X[[l]])) >= nodesize & ncol(data.frame(subtree_X[[l]])) > 0 ) {
                  ifelse(reuse==TRUE, parentsplit <- "skip", parentsplit <- subtree[[l]]$Var)
                  prep <- suppressWarnings(SplitPrep(Xdf=subtree_X[[l]], df=subtree_df[[l]], data_splt=subtree_df[[1]], parentsplit))
                  if (length(!is.na(data.frame(prep[[1]]))) >= 1) {
                    eval_temp <- NA
                    try(eval_temp <- suppressWarnings(MultiEval(X=data.frame(prep[[1]]), subtree_Z[[l]], data=prep[[2]], continuous = continuous, quantseq = quantseq, wt = wt, evalmethod = evalmethod, alpha = alpha, IGcutoff = IGcutoff, splitmin=splitmin)) , silent=TRUE)
                    # If the node i can be evaluated with the split variables and targets:
                    if (!is.na(eval_temp[2]) & (2*l+0 <= (2^(paralleldepth+1)-1)) & (2*l+1 <= (2^(paralleldepth+1)-1)) ) {
                      splitdf   <- eval_temp[[2]] #Extract splitvar with top rank
                      splitdf   <- (splitdf[!duplicated(splitdf[,c(1,3:ncol(splitdf))]), c(1:ncol(splitdf))])
                      temp      <- subtree_df[[l]]
                      splitdf   <- splitdf[order(splitdf$Rank),]
                      splitvar  <- as.character(head(splitdf[p,]$Split, n=1))
                      var       <- as.character(head(splitdf[p,]$Xm, n=1))
                      varnum    <- which(colnames(temp)==var)
                      notation  <- (unlist(strsplit(splitvar, " ")))
                      notation  <- notation[2]
                      thresh    <- str_remove(splitvar, var)
                      ifelse(notation == "in", threshold <- str_remove(thresh, " in "),
                             ifelse(notation == "<", threshold <- str_remove(thresh, " < "), NA))
                      # Subset data based on selected split
                      splitL_Z <- splitR_Z <- subtree_Z[[l]]
                      ifelse(notation == "in", splitL <- (temp[which(temp[,c(varnum)] == threshold), c(1:ncol(temp))]),
                             ifelse(notation == "<", splitL  <- (temp[which(temp[,c(varnum)] < as.numeric(threshold)), c(1:ncol(temp))]), NA))
                      splitL_Z$Z <- subset(splitL_Z$Z, row.names(splitL_Z$Z) %in% row.names(splitL))
                      ifelse(notation == "in", splitR <- (temp[which(temp[,c(varnum)] != threshold), c(1:ncol(temp))]),
                             ifelse(notation == "<", splitR <- (temp[which(temp[,c(varnum)] >= as.numeric(threshold)), c(1:ncol(temp))]), NA))
                      splitR_Z$Z <- subset(splitR_Z$Z, row.names(splitR_Z$Z) %in% row.names(splitR))
                      ifelse(notation == "in", splitvarR <- paste(var, threshold, sep = " not in "),
                             ifelse(notation == "<", splitvarR <- paste(var, threshold, sep = " >= "), NA))
                      # Evaluate child nodes across targets using MTSummary
                      splitL_eval <- suppressWarnings(MTSummary(splitL_Z, data=splitL))
                      splitR_eval <- suppressWarnings(MTSummary(splitR_Z, data=splitR))
                      # Evaluate split using CV error & sd (caret) --> Record for parent node
                      subtree[[l]]$PartVar <- var
                      misclass_targets_L <- misclass_targets_R <- as.list(1:length(subtree_Z[[2*l+0]]$Z))
                      for (j in 1:length(splitL_eval)) {
                        if(subtree_Z[[l]]$Definitions[j]=="Cat") {
                          splitL_eval[[j]]$relerr <- splitL_eval[[j]]$expectedloss/subtree[[1]]$Targets[[j]]$expectedloss
                        }
                        else if(subtree_Z[[l]]$Definitions[j] == "Surv") {
                          splitL_eval[[j]]$relerr <- splitL_eval[[j]]$deviance/subtree[[1]]$Targets[[j]]$deviance
                        }
                        else if(subtree_Z[[l]]$Definitions[j] == "Count") {
                          splitL_eval[[j]]$relerr <- splitL_eval[[j]]$deviance/subtree[[1]]$Targets[[j]]$deviance
                        }
                        else {
                          splitL_eval[[j]]$relerr <- splitL_eval[[j]]$MSE/subtree[[1]]$Targets[[j]]$MSE
                        }
                        misclass_targets_L[[j]] <- splitL_eval[[j]]$relerr
                      }
                      ifelse(!is.null(wt), (L_relerr <- weighted.mean(as.numeric(misclass_targets_L), w=(wt), na.rm=TRUE)),
                             (L_relerr <- mean(as.numeric(misclass_targets_L), na.rm = TRUE)))

                      for (j in 1:length(splitR_eval)) {
                        if(subtree_Z[[l]]$Definitions[j]=="Cat") {
                          splitR_eval[[j]]$relerr <- splitR_eval[[j]]$expectedloss/subtree[[1]]$Targets[[j]]$expectedloss
                        }
                        else if(subtree_Z[[l]]$Definitions[j] == "Surv") {
                          splitR_eval[[j]]$relerr <- splitR_eval[[j]]$deviance/subtree[[1]]$Targets[[j]]$deviance
                        }
                        else if(subtree_Z[[l]]$Definitions[j] == "Count") {
                          splitR_eval[[j]]$relerr <- splitR_eval[[j]]$deviance/subtree[[1]]$Targets[[j]]$deviance
                        }
                        else {
                          splitR_eval[[j]]$relerr <- splitR_eval[[j]]$MSE/subtree[[1]]$Targets[[j]]$MSE
                        }
                        misclass_targets_R[[j]] <- splitR_eval[[j]]$relerr
                      }
                      ifelse(!is.null(wt), (R_relerr <- weighted.mean(as.numeric(misclass_targets_R), w=(wt), na.rm=TRUE)),
                             (R_relerr <- mean(as.numeric(misclass_targets_R), na.rm = TRUE)))

                      subtree[[l]]$CP      <- (subtree[[l]]$relerr - ((L_relerr + R_relerr)/2))

                      subtree[[l]]$Eval <- splitdf
                      if (!is.na(subtree[[l]]$CP)) {
                        # Record child data frames, features, and targets
                        subtree_df[[2*l+0]] <- splitL
                        subtree_df[[2*l+1]] <- splitR
                        subtree_Z[[2*l+0]] <- splitL_Z
                        subtree_Z[[2*l+1]] <- splitR_Z
                        subtree_X[[2*l+0]] <- subset(subtree_X[[l]], row.names(subtree_X[[l]]) %in% row.names(splitL))
                        subtree_X[[2*l+1]] <- subset(subtree_X[[l]], row.names(subtree_X[[l]]) %in% row.names(splitR))
                        # Record left-child node
                        subtree[[2*l+0]]$NodeID  <- paste(2*l+0)
                        subtree[[2*l+0]]$SplitVar<- paste(splitvar)
                        subtree[[2*l+0]]$Var     <- paste(var)
                        subtree[[2*l+0]]$Thresh  <- paste(threshold)
                        subtree[[2*l+0]]$N       <- paste(nrow(subtree_df[[2*l+0]]))
                        subtree[[2*l+0]]$Pnode   <- paste(nrow(subtree_df[[2*l+0]])/nrow(subtree_df[[1]]))
                        subtree[[2*l+0]]$Targets <- splitL_eval
                        Xerr_targets <- Xstd_targets <- as.list(1:length(subtree_Z[[2*l+0]]$Z))
                        for (j in 1:length(misclass_targets)) {
                          Xerr_targets[[j]]     <- subtree[[2*l+0]]$Targets[[j]]$Xerror
                          Xstd_targets[[j]]     <- subtree[[2*l+0]]$Targets[[j]]$Xstd
                        }
                        subtree[[2*l+0]]$relerr <- L_relerr
                        # Record right-child node
                        subtree[[2*l+1]]$NodeID  <- paste(2*l+1)
                        subtree[[2*l+1]]$SplitVar<- paste(splitvarR)
                        subtree[[2*l+1]]$Var     <- paste(var)
                        subtree[[2*l+1]]$Thresh  <- paste(threshold)
                        subtree[[2*l+1]]$N       <- paste(nrow(subtree_df[[2*l+1]]))
                        subtree[[2*l+1]]$Pnode   <- paste(nrow(subtree_df[[2*l+1]])/nrow(subtree_df[[1]]))
                        subtree[[2*l+1]]$Targets <- splitR_eval
                        Xerr_targets <- Xstd_targets <- as.list(1:length(subtree_Z[[2*l+1]]$Z))
                        for (j in 1:length(misclass_targets)) {
                          Xerr_targets[[j]] <- subtree[[2*l+1]]$Targets[[j]]$Xerror
                          Xstd_targets[[j]] <- subtree[[2*l+1]]$Targets[[j]]$Xstd
                        }
                        subtree[[2*l+1]]$relerr <- R_relerr
                      }
                      else {
                        subtree[[l]]$NodeID  <- paste(subtree[[l]]$NodeID, "*")
                        subtree[[l]]$CP      <- NA
                      }
                    }
                    # If the node i cannot be evaluated, label node i as terminal (*) and move to next node:
                    else {
                      subtree[[l]]$NodeID  <- paste(subtree[[l]]$NodeID, "*")
                      subtree[[l]]$CP      <- NA
                    }
                  }
                  # Move on to next node: update i
                  l <- l + 1
                }
                # If node i is too small to partition (defined by nodesize), label node i as terminal (*) and move to next node:
                else {
                  subtree[[l]]$NodeID <- paste(subtree[[l]]$NodeID, "*")
                  subtree[[l]]$CP      <- NA
                  l <- l + 1
                }
              }
              # If node i doesn't exist, move on to next node:
              else {
                l <- l + 1
              }
            }
            paralleltrees[[p]] <- subtree
          }
          for (r in 1:length(pareval)) {
            relerrs <- vector(mode = "list", length = (2^(paralleldepth+1)-1))
            for (l in 1:length(relerrs)) {
              if (!is.null(paralleltrees[[r]][[l]])) {
                nodeid    <- (unlist(strsplit(paralleltrees[[r]][[l]]$NodeID, " ")))
                if (length(nodeid)>=2) {
                  nodeN     <- paralleltrees[[r]][[l]]$N
                  rootN     <- paralleltrees[[r]][[1]]$N
                  relerrs[[l]] <- paralleltrees[[r]][[l]]$relerr * (as.numeric(nodeN)/as.numeric(rootN))
                }
                else {
                  relerrs[[l]] <- NA
                }
              }
              else {
                relerrs[[l]] <- NA
              }
            }
            pareval[[r]] <- mean(as.numeric(unlist(relerrs)), na.rm = TRUE)
          }
          beststid  <- grep(min(as.numeric(unlist(pareval))), pareval)[1]
          bestst    <- paralleltrees[[beststid]]
          ifelse(!is.null(bestst[[2]]), splitvar  <- as.character(bestst[[2]]$SplitVar), splitvar <- NA)
          # If the node i can be evaluated with the split variables and targets:
          if (!is.na(splitvar) & (2*i+0 <= 2^(depth)) & (2*i+1 <= 2^(depth)) ) {
            temp      <- split_df[[i]]
            notation  <- (unlist(strsplit(splitvar, " ")))
            var       <- notation[1]
            varnum    <- which(colnames(temp)==var)
            notation  <- notation[2]
            thresh    <- str_remove(splitvar, var)
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
            # Evaluate child nodes across targets using MTSummary
            splitL_eval <- suppressWarnings(MTSummary(splitL_Z, data=splitL))
            splitR_eval <- suppressWarnings(MTSummary(splitR_Z, data=splitR))
            # Evaluate split using CV error & sd (caret) --> Record for parent node
            mtpart_splits[[i]]$PartVar <- var
            misclass_targets_L <- misclass_targets_R <- as.list(1:length(splitZ[[2*i+0]]$Z))

            for (j in 1:length(splitL_eval)) {
              if(splitZ[[i]]$Definitions[j]=="Cat") {
                splitL_eval[[j]]$relerr <- splitL_eval[[j]]$expectedloss/mtpart_splits[[1]]$Targets[[j]]$expectedloss
              }
              else if(splitZ[[i]]$Definitions[j] == "Surv") {
                splitL_eval[[j]]$relerr <- splitL_eval[[j]]$deviance/mtpart_splits[[1]]$Targets[[j]]$deviance
              }
              else if(splitZ[[i]]$Definitions[j] == "Count") {
                splitL_eval[[j]]$relerr <- splitL_eval[[j]]$deviance/mtpart_splits[[1]]$Targets[[j]]$deviance
              }
              else {
                splitL_eval[[j]]$relerr <- splitL_eval[[j]]$MSE/mtpart_splits[[1]]$Targets[[j]]$MSE
              }
              misclass_targets_L[[j]] <- splitL_eval[[j]]$relerr
            }
            ifelse(!is.null(wt), (L_relerr <- weighted.mean(as.numeric(misclass_targets_L), w=(wt), na.rm=TRUE)),
                   (L_relerr <- mean(as.numeric(misclass_targets_L), na.rm = TRUE)))

            for (j in 1:length(splitR_eval)) {
              if(splitZ[[i]]$Definitions[j]=="Cat") {
                splitR_eval[[j]]$relerr <- splitR_eval[[j]]$expectedloss/mtpart_splits[[1]]$Targets[[j]]$expectedloss
              }
              else if(splitZ[[i]]$Definitions[j] == "Surv") {
                splitR_eval[[j]]$relerr <- splitR_eval[[j]]$deviance/mtpart_splits[[1]]$Targets[[j]]$deviance
              }
              else if(splitZ[[i]]$Definitions[j] == "Count") {
                splitR_eval[[j]]$relerr <- splitR_eval[[j]]$deviance/mtpart_splits[[1]]$Targets[[j]]$deviance
              }
              else {
                splitR_eval[[j]]$relerr <- splitR_eval[[j]]$MSE/mtpart_splits[[1]]$Targets[[j]]$MSE
              }
              misclass_targets_R[[j]] <- splitR_eval[[j]]$relerr
            }
            ifelse(!is.null(wt), (R_relerr <- weighted.mean(as.numeric(misclass_targets_R), w=(wt), na.rm=TRUE)),
                   (R_relerr <- mean(as.numeric(misclass_targets_R), na.rm = TRUE)))

            propL <- (nrow(splitL_Z$Z)/nrow(splitZ[[i]]$Z))
            propR <- (nrow(splitR_Z$Z)/nrow(splitZ[[i]]$Z))

            mtpart_splits[[i]]$CP      <- mtpart_splits[[i]]$relerr - (((L_relerr*propL) + (R_relerr*propR))/2)
            # Evaluate split using CV error & sd (caret) --> Record for parent node
            splitCVeval <- CVsplitEval(splitvar, splitX[[i]], splitZ[[i]], split_df[[i]])
            mtpart_splits[[i]]$CVeval <- splitCVeval
            CVeval_Xerr_targets <- CVeval_Xstd_targets <- as.list(1:length(splitZ[[i]]$Z))
            for (j in 1:length(splitZ[[i]]$Z)) {
              if(splitZ[[i]]$Definitions[j]=="Cat") {
                mtpart_splits[[i]]$CVeval[[j]]$Xerror <- (1 - mtpart_splits[[i]]$CVeval[[j]]$Accuracy)
                mtpart_splits[[i]]$CVeval[[j]]$Xstd   <- mtpart_splits[[i]]$CVeval[[j]]$AccuracySD
              }
              else if(splitZ[[i]]$Definitions[j] == "Surv") {
                mtpart_splits[[i]]$CVeval[[j]]$Xerror <- (1 - mtpart_splits[[i]]$CVeval[[j]]$Rsquared)
                mtpart_splits[[i]]$CVeval[[j]]$Xstd   <- mtpart_splits[[i]]$CVeval[[j]]$RsquaredSD
              }
              else if(splitZ[[i]]$Definitions[j] == "Count") {
                mtpart_splits[[i]]$CVeval[[j]]$Xerror <- (1 - mtpart_splits[[i]]$CVeval[[j]]$Rsquared)
                mtpart_splits[[i]]$CVeval[[j]]$Xstd   <- mtpart_splits[[i]]$CVeval[[j]]$RsquaredSD
              }
              else {
                mtpart_splits[[i]]$CVeval[[j]]$Xerror <- (1 - mtpart_splits[[i]]$CVeval[[j]]$Rsquared)
                mtpart_splits[[i]]$CVeval[[j]]$Xstd   <- mtpart_splits[[i]]$CVeval[[j]]$RsquaredSD
              }
              CVeval_Xerr_targets[[j]] <- mtpart_splits[[i]]$CVeval[[j]]$Xerror
              CVeval_Xstd_targets[[j]] <- mtpart_splits[[i]]$CVeval[[j]]$Xstd
            }
            ifelse(!is.null(wt), (mtpart_splits[[i]]$SXerror <- (weighted.mean(as.numeric(CVeval_Xerr_targets), w=(wt), na.rm=TRUE))),
                   (mtpart_splits[[i]]$SXerror <- (mean(as.numeric(CVeval_Xerr_targets), na.rm=TRUE))))
            ifelse(!is.null(wt), (mtpart_splits[[i]]$SXstd <- ( weighted.mean(as.numeric(CVeval_Xstd_targets), w=(wt), na.rm=TRUE))),
                   (mtpart_splits[[i]]$SXstd <- ( mean(as.numeric(CVeval_Xstd_targets), na.rm=TRUE))))
            mtpart_splits[[i]]$Eval <- splitdf
            if (!is.null(mtpart_splits[[i]]$CP) && mtpart_splits[[i]]$CP >= cp) {
              # Record child data frames, features, and targets
              split_df[[2*i+0]] <- splitL
              split_df[[2*i+1]] <- splitR
              splitZ[[2*i+0]] <- splitL_Z
              splitZ[[2*i+1]] <- splitR_Z
              splitX[[2*i+0]] <- subset(splitX[[i]], row.names(splitX[[i]]) %in% row.names(splitL))
              splitX[[2*i+1]] <- subset(splitX[[i]], row.names(splitX[[i]]) %in% row.names(splitR))
              # Record left-child node
              mtpart_splits[[2*i+0]]$NodeID  <- paste(2*i+0)
              mtpart_splits[[2*i+0]]$SplitVar<- paste(splitvar)
              mtpart_splits[[2*i+0]]$Var     <- paste(var)
              mtpart_splits[[2*i+0]]$Thresh  <- paste(threshold)
              mtpart_splits[[2*i+0]]$N       <- paste(nrow(split_df[[2*i+0]]))
              mtpart_splits[[2*i+0]]$Pnode   <- paste(nrow(split_df[[2*i+0]])/nrow(split_df[[1]]))
              mtpart_splits[[2*i+0]]$Targets <- splitL_eval
              Xerr_targets <- Xstd_targets <- as.list(1:length(splitZ[[2*i+0]]$Z))
              for (j in 1:length(misclass_targets)) {
                Xerr_targets[[j]]     <- mtpart_splits[[2*i+0]]$Targets[[j]]$Xerror
                Xstd_targets[[j]]     <- mtpart_splits[[2*i+0]]$Targets[[j]]$Xstd
              }
              mtpart_splits[[2*i+0]]$relerr <- L_relerr
              # Record right-child node
              mtpart_splits[[2*i+1]]$NodeID  <- paste(2*i+1)
              mtpart_splits[[2*i+1]]$SplitVar<- paste(splitvarR)
              mtpart_splits[[2*i+1]]$Var     <- paste(var)
              mtpart_splits[[2*i+1]]$Thresh  <- paste(threshold)
              mtpart_splits[[2*i+1]]$N       <- paste(nrow(split_df[[2*i+1]]))
              mtpart_splits[[2*i+1]]$Pnode   <- paste(nrow(split_df[[2*i+1]])/nrow(split_df[[1]]))
              mtpart_splits[[2*i+1]]$Targets <- splitR_eval
              Xerr_targets <- Xstd_targets <- as.list(1:length(splitZ[[2*i+1]]$Z))
              for (j in 1:length(misclass_targets)) {
                Xerr_targets[[j]] <- mtpart_splits[[2*i+1]]$Targets[[j]]$Xerror
                Xstd_targets[[j]] <- mtpart_splits[[2*i+1]]$Targets[[j]]$Xstd
              }
              mtpart_splits[[2*i+1]]$relerr <- R_relerr
              # Record parent node id in child node id slots
              df[i,]$child         <- mtpart_splits[[i]]$NodeID
              df[(2*i+0),]$parent  <- mtpart_splits[[i]]$NodeID
              df[(2*i+1),]$parent  <- mtpart_splits[[i]]$NodeID
            }
            else {
              mtpart_splits[[i]]$NodeID  <- paste(mtpart_splits[[i]]$NodeID, "*")
              mtpart_splits[[i]]$SXerror <- NA
              mtpart_splits[[i]]$SXstd   <- NA
              mtpart_splits[[i]]$CP      <- NA
              df[i,]$child         <- mtpart_splits[[i]]$NodeID
            }
          }
          # If the node i cannot be evaluated, label node i as terminal (*) and move to next node:
          else {
            mtpart_splits[[i]]$NodeID  <- paste(mtpart_splits[[i]]$NodeID, "*")
            mtpart_splits[[i]]$SXerror <- NA
            mtpart_splits[[i]]$SXstd   <- NA
            mtpart_splits[[i]]$CP      <- NA

            df[i,]$child         <- mtpart_splits[[i]]$NodeID
          }
          # Move on to next node: update i
          i <- i + 1
        }
        # If node i is too small to partition (defined by nodesize), label node i as terminal (*) and move to next node:
        else {
          mtpart_splits[[i]]$NodeID <- paste(mtpart_splits[[i]]$NodeID, "*")
          mtpart_splits[[i]]$SXerror <- NA
          mtpart_splits[[i]]$SXstd   <- NA
          mtpart_splits[[i]]$CP      <- NA

          df[i,]$child         <- mtpart_splits[[i]]$NodeID
          i <- i + 1
        }
      }
      # If node i doesn't exist, move on to next node:
      else {
        i <- i + 1
      }
    }
  }
  else {
    for (i in 1:2^(depth)) {
      # If the node i exists:
      if (!is.null(mtpart_splits[[i]]) ) {
        # And if the node i is sufficiently large (defined by nodesize argument):
        if (!is.null(data.frame(splitX[[i]])) & nrow(data.frame(splitX[[i]])) >= nodesize & ncol(data.frame(splitX[[i]])) > 0 ) {
          ifelse(reuse==TRUE, parentsplit <- "skip", parentsplit <- mtpart_splits[[i]]$Var)
          prep <- suppressWarnings(SplitPrep(Xdf=splitX[[i]], df=split_df[[i]], data_splt=split_df[[1]], parentsplit))
          eval_temp <- NA
          suppressMessages(try(eval_temp <- suppressWarnings(MultiEval(X=data.frame(prep[[1]]), splitZ[[i]], data=prep[[2]], continuous = continuous, quantseq = quantseq, wt = wt, evalmethod = evalmethod, alpha = alpha, IGcutoff = IGcutoff, splitmin = splitmin)), silent=TRUE ))
          # If the node i can be evaluated with the split variables and targets:
          if (!is.na(eval_temp[2]) & (2*i+0 <= 2^(depth)) & (2*i+1 <= 2^(depth)) ) {
            splitdf   <- eval_temp[[2]] #Extract splitvar with top rank
            temp      <- split_df[[i]]
            splitvar  <- as.character(head(splitdf$Split[splitdf$Rank == min(splitdf$Rank)], n=1))
            var       <- as.character(head(splitdf$Xm[splitdf$Rank == min(splitdf$Rank)], n=1))
            varnum    <- which(colnames(temp)==var)
            notation  <- (unlist(strsplit(splitvar, " ")))
            notation  <- notation[2]
            thresh    <- str_remove(splitvar, var)
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
            # Evaluate child nodes across targets using MTSummary
            splitL_eval <- suppressWarnings(MTSummary(splitL_Z, data=splitL))
            splitR_eval <- suppressWarnings(MTSummary(splitR_Z, data=splitR))
            # Evaluate split using CV error & sd (caret) --> Record for parent node
            mtpart_splits[[i]]$PartVar <- var
            misclass_targets_L <- misclass_targets_R <- as.list(1:length(splitZ[[2*i+0]]$Z))

            for (j in 1:length(splitL_eval)) {
              if(splitZ[[i]]$Definitions[j]=="Cat") {
                splitL_eval[[j]]$relerr <- splitL_eval[[j]]$expectedloss/mtpart_splits[[1]]$Targets[[j]]$expectedloss
              }
              else if(splitZ[[i]]$Definitions[j] == "Surv") {
                splitL_eval[[j]]$relerr <- splitL_eval[[j]]$deviance/mtpart_splits[[1]]$Targets[[j]]$deviance
              }
              else if(splitZ[[i]]$Definitions[j] == "Count") {
                splitL_eval[[j]]$relerr <- splitL_eval[[j]]$deviance/mtpart_splits[[1]]$Targets[[j]]$deviance
              }
              else {
                splitL_eval[[j]]$relerr <- splitL_eval[[j]]$MSE/mtpart_splits[[1]]$Targets[[j]]$MSE
              }
              misclass_targets_L[[j]] <- splitL_eval[[j]]$relerr
            }
            ifelse(!is.null(wt), (L_relerr <- weighted.mean(as.numeric(misclass_targets_L), w=(wt), na.rm=TRUE)),
                   (L_relerr <- mean(as.numeric(misclass_targets_L), na.rm = TRUE)))

            for (j in 1:length(splitR_eval)) {
              if(splitZ[[i]]$Definitions[j]=="Cat") {
                splitR_eval[[j]]$relerr <- splitR_eval[[j]]$expectedloss/mtpart_splits[[1]]$Targets[[j]]$expectedloss
              }
              else if(splitZ[[i]]$Definitions[j] == "Surv") {
                splitR_eval[[j]]$relerr <- splitR_eval[[j]]$deviance/mtpart_splits[[1]]$Targets[[j]]$deviance
              }
              else if(splitZ[[i]]$Definitions[j] == "Count") {
                splitR_eval[[j]]$relerr <- splitR_eval[[j]]$deviance/mtpart_splits[[1]]$Targets[[j]]$deviance
              }
              else {
                splitR_eval[[j]]$relerr <- splitR_eval[[j]]$MSE/mtpart_splits[[1]]$Targets[[j]]$MSE
              }
              misclass_targets_R[[j]] <- splitR_eval[[j]]$relerr
            }
            ifelse(!is.null(wt), (R_relerr <- weighted.mean(as.numeric(misclass_targets_R), w=(wt), na.rm=TRUE)),
                   (R_relerr <- mean(as.numeric(misclass_targets_R), na.rm = TRUE)))

            propL <- (nrow(splitL_Z$Z)/nrow(splitZ[[i]]$Z))
            propR <- (nrow(splitR_Z$Z)/nrow(splitZ[[i]]$Z))

            mtpart_splits[[i]]$CP      <- mtpart_splits[[i]]$relerr - (((L_relerr*propL) + (R_relerr*propR))/2)
            # Evaluate split using CV error & sd (caret) --> Record for parent node
            splitCVeval <- CVsplitEval(splitvar, splitX[[i]], splitZ[[i]], split_df[[i]])
            mtpart_splits[[i]]$CVeval <- splitCVeval
            CVeval_Xerr_targets <- CVeval_Xstd_targets <- as.list(1:length(splitZ[[i]]$Z))
            for (j in 1:length(splitZ[[i]]$Z)) {
              if(splitZ[[i]]$Definitions[j]=="Cat") {
                mtpart_splits[[i]]$CVeval[[j]]$Xerror <- (1 - mtpart_splits[[i]]$CVeval[[j]]$Accuracy)
                mtpart_splits[[i]]$CVeval[[j]]$Xstd   <- mtpart_splits[[i]]$CVeval[[j]]$AccuracySD
              }
              else if(splitZ[[i]]$Definitions[j] == "Surv") {
                mtpart_splits[[i]]$CVeval[[j]]$Xerror <- (1 - mtpart_splits[[i]]$CVeval[[j]]$Rsquared)
                mtpart_splits[[i]]$CVeval[[j]]$Xstd   <- mtpart_splits[[i]]$CVeval[[j]]$RsquaredSD
              }
              else if(splitZ[[i]]$Definitions[j] == "Count") {
                mtpart_splits[[i]]$CVeval[[j]]$Xerror <- (1 - mtpart_splits[[i]]$CVeval[[j]]$Rsquared)
                mtpart_splits[[i]]$CVeval[[j]]$Xstd   <- mtpart_splits[[i]]$CVeval[[j]]$RsquaredSD
              }
              else {
                mtpart_splits[[i]]$CVeval[[j]]$Xerror <- (1 - mtpart_splits[[i]]$CVeval[[j]]$Rsquared)
                mtpart_splits[[i]]$CVeval[[j]]$Xstd   <- mtpart_splits[[i]]$CVeval[[j]]$RsquaredSD
              }
              CVeval_Xerr_targets[[j]] <- mtpart_splits[[i]]$CVeval[[j]]$Xerror
              CVeval_Xstd_targets[[j]] <- mtpart_splits[[i]]$CVeval[[j]]$Xstd
            }
            ifelse(!is.null(wt), (mtpart_splits[[i]]$SXerror <- (weighted.mean(as.numeric(CVeval_Xerr_targets), w=(wt), na.rm=TRUE))),
                   (mtpart_splits[[i]]$SXerror <- (mean(as.numeric(CVeval_Xerr_targets)))))
            ifelse(!is.null(wt), (mtpart_splits[[i]]$SXstd <- ( weighted.mean(as.numeric(CVeval_Xstd_targets), w=(wt), na.rm=TRUE))),
                   (mtpart_splits[[i]]$SXstd <- ( mean(as.numeric(CVeval_Xstd_targets)))))
            mtpart_splits[[i]]$Eval <- splitdf
            if (!is.null(mtpart_splits[[i]]$CP) && !is.na(mtpart_splits[[i]]$CP) && mtpart_splits[[i]]$CP >= cp) {
              # Record child data frames, features, and targets
              split_df[[2*i+0]] <- splitL
              split_df[[2*i+1]] <- splitR
              splitZ[[2*i+0]] <- splitL_Z
              splitZ[[2*i+1]] <- splitR_Z
              splitX[[2*i+0]] <- subset(splitX[[i]], row.names(splitX[[i]]) %in% row.names(splitL))
              splitX[[2*i+1]] <- subset(splitX[[i]], row.names(splitX[[i]]) %in% row.names(splitR))
              # Record left-child node
              mtpart_splits[[2*i+0]]$NodeID  <- paste(2*i+0)
              mtpart_splits[[2*i+0]]$SplitVar<- paste(splitvar)
              mtpart_splits[[2*i+0]]$Var     <- paste(var)
              mtpart_splits[[2*i+0]]$Thresh  <- paste(threshold)
              mtpart_splits[[2*i+0]]$N       <- paste(nrow(split_df[[2*i+0]]))
              mtpart_splits[[2*i+0]]$Pnode   <- paste(nrow(split_df[[2*i+0]])/nrow(split_df[[1]]))
              mtpart_splits[[2*i+0]]$Targets <- splitL_eval
              Xerr_targets <- Xstd_targets <- as.list(1:length(splitZ[[2*i+0]]$Z))
              for (j in 1:length(misclass_targets)) {
                Xerr_targets[[j]]     <- mtpart_splits[[2*i+0]]$Targets[[j]]$Xerror
                Xstd_targets[[j]]     <- mtpart_splits[[2*i+0]]$Targets[[j]]$Xstd
              }
              mtpart_splits[[2*i+0]]$relerr <- L_relerr
              # Record right-child node
              mtpart_splits[[2*i+1]]$NodeID  <- paste(2*i+1)
              mtpart_splits[[2*i+1]]$SplitVar<- paste(splitvarR)
              mtpart_splits[[2*i+1]]$Var     <- paste(var)
              mtpart_splits[[2*i+1]]$Thresh  <- paste(threshold)
              mtpart_splits[[2*i+1]]$N       <- paste(nrow(split_df[[2*i+1]]))
              mtpart_splits[[2*i+1]]$Pnode   <- paste(nrow(split_df[[2*i+1]])/nrow(split_df[[1]]))
              mtpart_splits[[2*i+1]]$Targets <- splitR_eval
              Xerr_targets <- Xstd_targets <- as.list(1:length(splitZ[[2*i+1]]$Z))
              for (j in 1:length(misclass_targets)) {
                Xerr_targets[[j]] <- mtpart_splits[[2*i+1]]$Targets[[j]]$Xerror
                Xstd_targets[[j]] <- mtpart_splits[[2*i+1]]$Targets[[j]]$Xstd
              }
              mtpart_splits[[2*i+1]]$relerr <- R_relerr
              # Record parent node id in child node id slots
              df[i,]$child         <- mtpart_splits[[i]]$NodeID
              df[(2*i+0),]$parent  <- mtpart_splits[[i]]$NodeID
              df[(2*i+1),]$parent  <- mtpart_splits[[i]]$NodeID
            }
            else {
              mtpart_splits[[i]]$NodeID  <- paste(mtpart_splits[[i]]$NodeID, "*")
              mtpart_splits[[i]]$SXerror <- NA
              mtpart_splits[[i]]$SXstd   <- NA
              mtpart_splits[[i]]$CP      <- NA

              df[i,]$child         <- mtpart_splits[[i]]$NodeID
            }
          }
          # If the node i cannot be evaluated, label node i as terminal (*) and move to next node:
          else {
            mtpart_splits[[i]]$NodeID  <- paste(mtpart_splits[[i]]$NodeID, "*")
            mtpart_splits[[i]]$SXerror <- NA
            mtpart_splits[[i]]$SXstd   <- NA
            mtpart_splits[[i]]$CP      <- NA

            df[i,]$child         <- mtpart_splits[[i]]$NodeID
          }
          # Move on to next node: update i
          i <- i + 1
        }
        # If node i is too small to partition (defined by nodesize), label node i as terminal (*) and move to next node:
        else {
          mtpart_splits[[i]]$NodeID <- paste(mtpart_splits[[i]]$NodeID, "*")
          mtpart_splits[[i]]$SXerror <- NA
          mtpart_splits[[i]]$SXstd   <- NA
          mtpart_splits[[i]]$CP      <- NA

          df[i,]$child         <- mtpart_splits[[i]]$NodeID
          i <- i + 1
        }
      }
      # If node i doesn't exist, move on to next node:
      else {
        i <- i + 1
      }
    }
  }
  return(list(partitions = df, tree_nodes = mtpart_splits)) #consider outputting eval directly instead of raw datasets
}

