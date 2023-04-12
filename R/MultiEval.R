#' @name MultiEval
#' @title Evaluate features across multiple outcomes and outcome types
#'
#' @import matrixStats
#' @import survival
#' @import reshape2
#'
#' @description Evaluates multiple outcome variables of different types for a set of features.
#' \code{MultiEval} evaluates multiple outcomes of varying types across all features.
#' \code{Cat_Eval} performs evaluation of a categorical outcome across all features.
#' \code{Cont_Eval} performs evaluation of a continuous outcome across all features.
#' \code{Count_Eval} performs evaluation of an event rate or count outcome across all features.
#' \code{Surv_Eval} performs evaluation of a time-to-event or survival outcome across all features
#'
#' @param X Data frame object of all features to be evaluated.
#' @param m Specific feature to evaluate. Used only in \code{Cat_Eval, Cont_Eval,
#' Count_Eval, Surv_Eval}.
#' @param Ztype Vector containing a list of outcome types ("Cat", "Cont", "Count", "Surv")
#' and a data frame object of all outcomes being evaluated.
#' @param Z Data frame object of all outcomes being evaluated. Used only in \code{Cat_Eval,
#' Cont_Eval, Count_Eval, Surv_Eval}.
#' @param k Specific outcome to evaluate. Used only in \code{Cat_Eval, Cont_Eval, Count_Eval,
#' Surv_Eval}.
#' @param data The full dataset to be evaluated.
#' @param continuous How to handle partition evaluation of continuous features, in which
#' all potential split candidates may be examined ("all") or split candidates are derived
#' from quantiles ("quantile").
#' @param quantseq A user specified sequence ranging from 0 to 1 from which quantiles may
#' be specified (ex: seq(0,1,0.05)). If \code{continuous} and \code{quantseq} are not
#' specified the function will default to the "all" setting for continuous variables with
#' fewer than 100 unique values. If there are more than 100 unique values, the function
#' will use the "quantile" option with the default setting of seq(0,1,0.05).
#' @param wt The outcome variable weights to be applied in multi-target evaluations.
#' Default weights are set to \code{NULL}, indicating equal target variable importance.
#' @param evalmethod The method used to evaluate standardized proportion of information
#' gain across targets and target types.
#' @param alpha The significance level to be used in \code{evalmethod} if a p-value based
#' method (such as ) is selected. The default value for \code{alpha} is set to 0.05.
#' @param IGcutoff The IG cutoff threshold to be used for \code{evalmethod} if methods
#' xxx or yyyy are selected. The default value for \code{IGcutoff} is 0.95.
#' @param splitmin The minimum sample size \code{n} required for a node to be partitioned.
#' This is passed internally from \code{MTPart} and is used to assess feasibility of
#' feature binarization for evaluation (i.e., a binarized feature must result in
#' partitions with sufficient observations in both child nodes). The default value is of
#' \code{splitmin} is set to 10.
#'
#' @details The following function requires inputs including all features to be evaluated
#' (\code{X}), which specific feature to evaluate (\code{m}), all outcomes being evaluated
#' (\code{Z}), the specific outcome to evaluate (\code{k}), and the data (\code{data}).
#' In addition, there are options for the user to specify how to handle partition evaluation
#' of continuous features, in which all potential split candidates may be examined ("all")
#' or split candidates are derived from quantiles ("quantile"); the option \code{quantseq}
#' takes a user specified sequence ranging from 0 to 1 from which quantiles may be specified
#' (ex: seq(0,1,0.05)). If \code{continuous} and \code{quantseq} are not manually specified
#' the function will default to the "all" setting for continuous variables with fewer than
#' 100 unique values. If there are more than 100 unique values, the function will use the
#' "quantile" option with the default setting of seq(0,1,0.05).
#'
#' The function outputs a matrix including information about the outcome (Zk), outcome
#' variable type (Ztype), the feature (Xm), feature type (Xtype), the split threshold
#' being evaluated (Split), evaluation method (Eval), and the entropy of a partition
#' using that split threshold (SplitEntropy), information gain relative to the root node
#' (Inf_Gain) and the proportion of information gain (Inf_GainProp), the significance of
#' association between the split threshold and outcome (P_Value), relative ranking of split
#' thresholds based on how much their respective partitions reduced entropy (Rank), and
#' the standard normalized proportion of information gain associated with each threshold.
#'
#' All evaluation functions operate internally within the multi-target evaluation
#' function in which all outcomes and all features are assessed at a node. Therefore,
#' m and k options are specified by a loop from the user supplied inputs: X, Z, and data.
#' Again, continuous and quantseq are optional arguments initialized by the user.
#'
#' @return The \code{MultiEval} function returns a vector including of the evaluation
#' metrics for each possible split of each feature in the dataset.
#'
#' @usage MultiEval(X, Ztype, data, continuous = "quantile",
#'           quantseq = seq(0,1,0.05), wt=NULL,
#'           evalmethod = "avgIG", alpha = 0.05,
#'           IGcutoff = 0.95, splitmin)
#' @export
#' @aliases MultiEval
#' @rdname MultiEval
MultiEval <- function(X, Ztype, data, continuous = "quantile", quantseq = seq(0,1,0.05), wt=NULL, evalmethod = "avgIG", alpha = 0.05, IGcutoff = 0.95, splitmin=10) {
  Z <- Ztype$Z
  Eval <- matrix(ncol = 12, nrow = 0)
  colnames(Eval) <- c("Zk","Ztype","Xm","Xtype","Split","Eval","SplitEntropy","Inf_Gain","Inf_GainProp","P_Value","Rank","Scale" )
  for (k in 1:length(Z)) {
    if (Ztype$Definitions[k]=="Cont") {
      for (m in 1:length(X)) {
        suppressMessages(try(evaltemp <- Cont_Eval(X, m, Z, k, data, continuous, quantseq, splitmin), silent = TRUE))
        if (!is.na(evaltemp[9])) {
          Eval <- rbind(Eval, evaltemp)
        }
        else {
          Eval <- Eval
        }
      }
    }
    else if (Ztype$Definitions[k]=="Cat") {
      for (m in 1:length(X)) {
        evaltemp <- Cat_Eval(X, m, Z, k, data, continuous, quantseq, splitmin)
        if (!is.na(evaltemp[9])) {
          Eval <- rbind(Eval, evaltemp)
        }
        else {
          Eval <- Eval
        }
      }
    }
    else if (Ztype$Definitions[k]=="Surv") {
      for (m in 1:length(X)) {
        evaltemp <- Surv_Eval(X, m, Z, k, data, continuous, quantseq, splitmin)
        if (!is.na(evaltemp[9])) {
          Eval <- rbind(Eval, evaltemp)
        }
        else {
          Eval <- Eval
        }
      }
    }
    else {
      for (m in 1:length(X)) {
        evaltemp <- Count_Eval(X, m, Z, k, data, continuous, quantseq, splitmin)
        if (!is.na(evaltemp[9])) {
          Eval <- rbind(Eval, evaltemp)
        }
        else {
          Eval <- Eval
        }
      }
    }
  }
  # Define feature/outcome ranking will be performed based on information gain
  # Apply standard normal scaling N(0,1) to proportion information gain
  Eval[,12] <- NA
  if (evalmethod != "partErr") {
    for (k in 1:length(Z)) {
      kScale <- Eval[which(Eval[,1]==names(Z[k])), ]
      if (ncol(data.frame(kScale)) <= 1) {
        kScale <- t(data.frame(kScale))
        forScale <- as.numeric(as.character(kScale[,9]))
        Eval[which(Eval[,1]==names(Z[k])), ][12] <- forScale
      }
      else {
        kScale <- kScale
        forScale <- as.numeric(as.character(kScale[,9]))
        scalelist <- scale(forScale)
        Eval[which(Eval[,1]==names(Z[k])), ][,12] <- scalelist[,1]
      }
    }
  }
  else {
    for (k in 1:length(Z)) {
      kScale <- Eval[which(Eval[,1]==names(Z[k])), ]
      if (ncol(data.frame(kScale)) <= 1) {
        kScale <- t(data.frame(kScale))
        forScale <- as.numeric(as.character(kScale[,7]))
        Eval[which(Eval[,1]==names(Z[k])), ][12] <- forScale
      }
      else {
        kScale <- kScale
        forScale <- as.numeric(as.character(kScale[,7]))
        scalelist <- scale(forScale)
        Eval[which(Eval[,1]==names(Z[k])), ][,12] <- scalelist[,1]
      }
    }
  }
  # Define total ranking (Rank)
  forRank <- as.numeric(as.character(Eval[,12]))
  Eval[,11] <- as.numeric(as.character(Eval[,11]))
  ifelse(evalmethod != "partErr", (Eval[,11] <- (rank(-forRank))), (Eval[,11] <- (rank(forRank))))
  df_Eval       <- data.frame(Eval)
  df_Eval$Scale <- as.numeric(as.character(df_Eval$Scale))
  suppressMessages(ScaledInfGain         <- reshape2::dcast(df_Eval, Xm + Split ~ Zk, value.var="Scale"))
  #Proposed methods to define best splits
  ScaledInfGain$AvgIG   <- matrixStats::rowWeightedMeans(as.matrix(ScaledInfGain[,c(3:ncol(ScaledInfGain))]), weights=wt, na.rm=TRUE)
  ScaledInfGain$MaxIG   <- NA
  tempIGVals <- ScaledInfGain[3:(which(names(ScaledInfGain)=="AvgIG")-1)]
  for (j in 1:ncol(ScaledInfGain[3:(which(names(ScaledInfGain)=="AvgIG")-1)])) {
    IGcutoff_val <- quantile(ScaledInfGain[,(2+j)], probs = IGcutoff, na.rm = TRUE)[[1]]
    tempIGVals[,j] <- ifelse(ScaledInfGain[,(2+j)] >= IGcutoff_val, 1*(ScaledInfGain[,(2+j)]), 0)
  }
  ScaledInfGain$MostIG  <- rowSums(tempIGVals)
  ScaledInfGain$AvgPVal  <- (1-(pnorm(ScaledInfGain$AvgIG)))
  ScaledInfGain$MinPVal <- NA
  for (j in 1:nrow(ScaledInfGain)) {
    ScaledInfGain[j,]$MaxIG   <- max(ScaledInfGain[j, 3:(which(names(ScaledInfGain)=="AvgIG")-1)])
    ScaledInfGain[j,]$MinPVal <- min(1-(pnorm(as.numeric(as.character(ScaledInfGain[j, 3:(which(names(ScaledInfGain)=="AvgIG")-1)])))))
  }
  ScaledInfGain$MostPVal    <- NA
  ScaledInfGain$AvgMostPVal <- NA
  tempPVals <- tempPValsSig <- ScaledInfGain[3:(which(names(ScaledInfGain)=="AvgIG")-1)]
  for (l in 1:ncol(ScaledInfGain[3:(which(names(ScaledInfGain)=="AvgIG")-1)])) {
    if (!is.na(sum(ScaledInfGain[,(2+l)]))) {
      vals <- ScaledInfGain[,(2+l)]
      tempPVals[,l] <- 1-pnorm(vals)
      for (r in 1:nrow(tempPVals)) {
        if (tempPVals[r,l] <= alpha) {
          tempPValsSig[r,l] <- tempPVals[r,l]
        }
        else {
          tempPValsSig[r,l] <- NA
        }
      }
    }
    else {
      next
    }
  }
  ScaledInfGain$MostPVal    <- rowSums(tempPVals <= alpha)
  ScaledInfGain$AvgMostPVal <- rowMeans(tempPValsSig, na.rm = TRUE)
  ScaledInfGain$InvMinPVal  <- (1 - ScaledInfGain$AvgMostPVal)
  ScaledInfGain$WtMinPval   <- (ScaledInfGain$InvMinPVal)*(ScaledInfGain$MostPVal)
  if (evalmethod == "maxIG") {
    ScaledInfGain$Rank    <- rank(-ScaledInfGain$MaxIG) #Maximum information gain on any one target
  }
  else if (evalmethod == "mostIG") {
    ScaledInfGain$Rank    <- rank(-ScaledInfGain$MostIG) #Most number of meaningful (>= IGcutoff) targets
  }
  else if (evalmethod == "avgPVal") {
    if (all(is.na(ScaledInfGain$AvgMostPVal))) {
      ScaledInfGain$Rank    <- rank(ScaledInfGain$AvgPVal) #Average p-value of all significant (<= alpha) targets
    }
    # If no targets are significant --> avg. sig = NA --> default to avg. p-val
    else {
      ScaledInfGain$Rank    <- rank(ScaledInfGain$AvgMostPVal) #Average p-value of all significant (<= alpha) targets
    }
  }
  else if (evalmethod == "minPVal") {
    if (all(is.na(ScaledInfGain$AvgMostPVal))) {
      ScaledInfGain$Rank    <- rank(ScaledInfGain$MinPVal)
    }
    else {
      ScaledInfGain$Rank    <- rank(-ScaledInfGain$WtMinPval) #Average target-weighted p-value of all significant targets
    }
  }
  else if (evalmethod == "mostPVal") {
    if (all(is.na(ScaledInfGain$AvgMostPVal))) {
      ScaledInfGain$Rank    <- rank(ScaledInfGain$MostIG)
    }
    # If no targets are significant --> avg. sig = NA --> default to most targets with IG >= IGcutoff
    else {
      ScaledInfGain$Rank    <- rank(-ScaledInfGain$MostPVal) #Average target-weighted p-value of all significant targets
    }
  }
  else if (evalmethod == "partErr") {
    ScaledInfGain$Rank    <- rank(ScaledInfGain$AvgIG) #Most information gain (average) across targets
  }
  else {
    ScaledInfGain$Rank    <- rank(-ScaledInfGain$AvgIG) #Most information gain (average) across targets
  }
  return(list(Eval, ScaledInfGain))
}

#' @usage Cat_Eval(X, m, Z, k, data, continuous,
#'          quantseq = seq(0,1,0.05), splitmin=10)
#' @export
#' @aliases Cat_Eval
#' @rdname MultiEval
Cat_Eval <- function(X, m, Z, k, data, continuous, quantseq = seq(0,1,0.05), splitmin=10) {
  contFeature <- list("Zk"=colnames(Z)[k], "Ztype"="Categorical", "Xm"=colnames(X)[m], "Xtype"=NA, "Eval"=NA)
  znum <- which( colnames(data)==names(Z[k]) )                           # Extract column number of outcome k
  xnum <- which( colnames(data)==names(X[m]) )                           # Extract column number of feature m
  Xm <- data[complete.cases(data[xnum]), c(xnum)]                     # Extract feature m from the sample or subsample and eliminate missing values
  Xtype <- class(Xm)                                                        # Extract feature attributes
  # Find root entropy
  df <- data[complete.cases(data[,c(xnum,znum)]), c(1:ncol(data))]
  rootProb <- table(Z[,k] )/length(Z[,k] )
  rootGini <- (1-sum(rootProb**2))
  # Case of Continuous X Feature
  if (Xtype == "numeric" || Xtype =="integer") {
    contFeature$Xtype <- "Continuous"
    l <- sort(unique(Xm) )                                                  # Extract sequence of unique m values (split threshold candidates)
    # Check that split candidates result in partitions with sufficient opservations (n >= 30)
    valcheck <- data.frame(matrix(nrow = length(l), ncol = 3))
    valcheck[,1] <- l
    for (c in 1:length(l)) {
      valcheck[c,2] <- sum(table(Xm[which(Xm < l[c])]))
      valcheck[c,3] <- sum(table(Xm[which(Xm >= l[c])]))
    }
    ln <- valcheck[which(valcheck[,2] >= splitmin & valcheck[,3] >= splitmin),1]  # Define new split candidate list based on above check
    if (length(ln)>0) {
      # If there are less than or equal to 100 split candidates, assess all
      if (continuous == "all") {
        ln <- ln
      }
      # If there are more than 100 split candidates, assess quantiles (default: each 5% seq(0,1,0.05))
      else {
        ln <- (quantile(ln, probs = quantseq))
      }
      contFeature$Eval <- matrix(nrow = length(ln), ncol = 12)                   # Create matrix to store evaluations on each unique m value
      colnames(contFeature$Eval) <- c("Zk","Ztype","Xm","Xtype","Split","Eval","SplitEntropy", "Inf_Gain", "Inf_GainProp", "P_Value", "Rank", "Scale")
      contFeature$Eval[,1] <- contFeature$Zk
      contFeature$Eval[,2] <- contFeature$Ztype
      contFeature$Eval[,3] <- colnames(X)[m]
      contFeature$Eval[,4] <- contFeature$Xtype
      contFeature$Eval[,5] <- paste(colnames(X)[m],ln,sep=" < ",collapse=NULL)   # Set matrix row names to unique m values
      contFeature$Eval[,6] <- "Gini"
      for (i in 1:length(ln)) {
        split     <- ln[i]                                                                                 # Extract threshold value i
        splitvar  <- data$splitvar <- ifelse(data[xnum]<split, "Left", ifelse(data[xnum]>=split, "Right", NA))
        splitnum <- which( colnames(data)=="splitvar" )                           # Extract column number of feature m

        base_prob <- table(splitvar)/length(splitvar)
        crosstab  <- table(data[,znum],data[,splitnum])
        crossprob <- prop.table(crosstab,2)
        No_Node_Gini  <- 1-sum(crossprob[,1]**2)                                                          # Assess left child node
        Yes_Node_Gini <- 1-sum(crossprob[,2]**2)                                                          # Assess right child node
        contFeature$Eval[i,7] <- mean((base_prob[[1]]*min(crossprob[,1])), (base_prob[[2]]*min(crossprob[,2])))                          # Assess entropy of split
        contFeature$Eval[i,8] <- (rootGini - sum(base_prob * c(No_Node_Gini,Yes_Node_Gini)))              # Assess information gain of split
        contFeature$Eval[i,9] <- ((rootGini - sum(base_prob * c(No_Node_Gini,Yes_Node_Gini)))/rootGini)   # Assess information gain of split
        contFeature$Eval[i,10] <- NA                                                                      # Assess average p-value of split
      }
    }
    else {
      contFeature$Eval <- NA
    }
  }
  # Case of Categorical X Feature
  else {
    if (length(levels(Xm)) == 2) {
      contFeature$Xtype <- "Categorical"
      l <- levels(Xm)
      contFeature$Eval <- matrix(nrow = length(l)-1, ncol = 12)                     # Create matrix to store evaluations on each unique m value
      colnames(contFeature$Eval) <- c("Zk","Ztype","Xm","Xtype","Split","Eval","SplitEntropy", "Inf_Gain", "Inf_GainProp", "P_Value", "Rank", "Scale")
      contFeature$Eval[,1] <- contFeature$Zk
      contFeature$Eval[,2] <- contFeature$Ztype
      contFeature$Eval[,3] <- colnames(X)[m]
      contFeature$Eval[,4] <- contFeature$Xtype
      contFeature$Eval[,5] <- paste(colnames(X)[m],l[1],sep=" in ",collapse=NULL)   # Set matrix row names to unique m values
      contFeature$Eval[,6] <- "Gini"
      #For all levels i of feature m, where m has 2 levels:
      split     <- l[1]                                                                                 # Extract threshold value i
      splitvar  <- data$splitvar <- ifelse(data[xnum]==split, "Left", ifelse(data[xnum] != split, "Right", NA))
      splitnum <- which( colnames(data)=="splitvar" )                           # Extract column number of feature m

      base_prob <- table(splitvar)/length(splitvar)
      crosstab  <- table(data[,znum],data[,splitnum])
      crossprob <- prop.table(crosstab,2)
      No_Node_Gini  <- 1-sum(crossprob[,1]**2)                                                          # Assess left child node
      Yes_Node_Gini <- 1-sum(crossprob[,2]**2)                                                          # Assess right child node
      contFeature$Eval[,7] <- mean((base_prob[[1]]*min(crossprob[,1])), (base_prob[[2]]*min(crossprob[,2])))                          # Assess entropy of split
      contFeature$Eval[,8] <- (rootGini - sum(base_prob * c(No_Node_Gini,Yes_Node_Gini)))              # Assess information gain of split
      contFeature$Eval[,9] <- ((rootGini - sum(base_prob * c(No_Node_Gini,Yes_Node_Gini)))/rootGini)   # Assess information gain of split
      contFeature$Eval[,10] <- NA                                                                      # Assess average p-value of split
    }
    else if (length(levels(Xm)) > 2) {
      contFeature$Xtype <- "Categorical"
      l <- levels(Xm)
      contFeature$Eval <- matrix(nrow = length(l), ncol = 12)         # Create matrix to store evaluations on each unique m value
      colnames(contFeature$Eval) <- c("Zk","Ztype","Xm","Xtype","Split","Eval","SplitEntropy", "Inf_Gain", "Inf_GainProp", "P_Value", "Rank", "Scale")
      contFeature$Eval[,1] <- contFeature$Zk
      contFeature$Eval[,2] <- contFeature$Ztype
      contFeature$Eval[,3] <- colnames(X)[m]
      contFeature$Eval[,4] <- contFeature$Xtype
      contFeature$Eval[,5] <- paste(colnames(X)[m],l,sep=" in ",collapse=NULL)   # Set matrix row names to unique m values
      contFeature$Eval[,6] <- "Gini"
      #For all levels i of feature m, where m has 2 levels:
      for (i in 1:length(levels(Xm))) {
        split     <- l[i]                                                                                 # Extract threshold value i
        splitvar  <- ifelse(data[xnum]==split, "Left", ifelse(data[xnum] != split, "Right", NA))
        base_prob <- table(splitvar)/length(splitvar)
        crosstab  <- table(Z[,k],splitvar)
        crossprob <- prop.table(crosstab,2)
        No_Node_Gini  <- 1-sum(crossprob[,1]**2)                                                          # Assess left child node
        Yes_Node_Gini <- 1-sum(crossprob[,2]**2)                                                          # Assess right child node
        contFeature$Eval[i,7] <- mean((base_prob[[1]]*min(crossprob[,1])), (base_prob[[2]]*min(crossprob[,2])))                          # Assess entropy of split
        contFeature$Eval[i,8] <- (rootGini - sum(base_prob * c(No_Node_Gini,Yes_Node_Gini)))              # Assess information gain of split
        contFeature$Eval[i,9] <- ((rootGini - sum(base_prob * c(No_Node_Gini,Yes_Node_Gini)))/rootGini)   # Assess information gain of split
        contFeature$Eval[i,10] <- NA                                                                      # Assess average p-value of split
      }
    }
  }
  # Define feature/outcome ranking will be performed based on information gain
  if (!is.na(contFeature$Eval[1])) {
    forRank <- as.numeric(as.character(contFeature$Eval[,8]))
    contFeature$Eval[,11] <- (rank(-forRank))                                             # Define raw ranking (kRank)
    forScale <- as.numeric(as.character(contFeature$Eval[,9]))
    contFeature$Eval[,12] <- scale(forScale)         # Apply standard normal scaling N(0,1) to proportion information gain
  }
  else {
    contFeature$Eval <- NA
  }
  # Output standardized information gain for SplitVars on k to matrix of splitvars rows, and k columns
  return(contFeature$Eval)
}

#' @usage Cont_Eval(X, m, Z, k, data, continuous,
#'           quantseq = seq(0,1,0.05), splitmin=10)
#' @export
#' @aliases Cont_Eval
#' @rdname MultiEval
Cont_Eval <- function(X, m, Z, k, data, continuous, quantseq = seq(0,1,0.05), splitmin=10) {
  contFeature <- list("Zk"=colnames(Z)[k], "Ztype"="Continuous", "Xm"=colnames(X)[m], "Xtype"=NA, "Eval"=NA)
  znum <- which( colnames(data)==names(Z[k]) )                           # Extract column number of outcome k
  xnum <- which( colnames(data)==names(X[m]) )                           # Extract column number of feature m
  Xm <- data[complete.cases(data[xnum]), c(xnum)]                     # Extract feature m from the sample or subsample and eliminate missing values
  Xtype <- class(Xm)                                                        # Extract feature attributes
  # Find root entropy
  df <- data[complete.cases(data[,c(xnum,znum)]), c(1:ncol(data))]
  root_fit <- lm(df[,(znum)] ~ 1, df)
  Root_Deviance <- deviance(root_fit)
  # Case of Continuous X Feature
  if (Xtype == "numeric" || Xtype == "integer") {
    contFeature$Xtype <- "Continuous"
    l <- sort(unique(Xm) )                                                  # Extract sequence of unique m values (split threshold candidates)
    # Check that split candidates result in partitions with sufficient opservations (n >= 30)
    valcheck <- data.frame(matrix(nrow = length(l), ncol = 3))
    valcheck[,1] <- l
    for (c in 1:length(l)) {
      valcheck[c,2] <- sum(table(Xm[which(Xm < l[c])]))
      valcheck[c,3] <- sum(table(Xm[which(Xm >= l[c])]))
    }
    ln <- valcheck[which(valcheck[,2] >= splitmin & valcheck[,3] >= splitmin),1]  # Define new split candidate list based on above check
    # If there are less than or equal to 100 split candidates, assess all
    if (length(ln) <= 100 & length(ln) > 0 || continuous == "all") {
      ln <- ln
    }
    # If there are more than 100 split candidates, assess quantiles (default: each 5% seq(0,1,0.05))
    if (length(ln) > 100 || continuous == "quantile") {
      ln <- (quantile(ln, probs = quantseq))
    }
    contFeature$Eval <- matrix(nrow = length(ln), ncol = 12)         # Create matrix to store evaluations on each unique m value
    colnames(contFeature$Eval) <- c("Zk","Ztype","Xm","Xtype","Split","Eval","SplitEntropy", "Inf_Gain", "Inf_GainProp", "P_Value", "Rank", "Scale")
    contFeature$Eval[,1] <- contFeature$Zk
    contFeature$Eval[,2] <- contFeature$Ztype
    contFeature$Eval[,3] <- colnames(X)[m]
    contFeature$Eval[,4] <- contFeature$Xtype
    contFeature$Eval[,5] <- paste(colnames(X)[m],ln,sep=" < ",collapse=NULL)   # Set matrix row names to unique m values
    contFeature$Eval[,6] <- "Deviance"
    if (length(ln) > 0) {
      for (i in 1:length(ln)) {
        split <- ln[i]                                                                                 # Extract threshold value i
        df_left            <- data[which(data[xnum]<split),c(1:ncol(data))]                   # Assess left child node
        left_fit           <- lm(df_left[,(znum)] ~ df_left[,(xnum)], df_left)
        lf                 <- summary(left_fit)
        Left_Deviance      <- deviance(left_fit)
        df_right           <- data[which(data[xnum]>=split),c(1:ncol(data))]                  # Assess right child node
        right_fit          <- lm(df_right[,(znum)] ~ df_right[,(xnum)], df_right)
        rf                 <- summary(right_fit)
        Right_Deviance     <- deviance(right_fit)
        contFeature$Eval[i,7] <- mean(((nrow(df_left)/nrow(data))*Left_Deviance), ((nrow(df_right)/nrow(data))*Right_Deviance))                                 # Assess entropy of split
        contFeature$Eval[i,8] <- (Root_Deviance - sum(((nrow(df_left)/nrow(data))*Left_Deviance), ((nrow(df_right)/nrow(data))*Right_Deviance)))                  # Assess information gain of split
        contFeature$Eval[i,9] <- ((Root_Deviance - sum(((nrow(df_left)/nrow(data))*Left_Deviance), ((nrow(df_right)/nrow(data))*Right_Deviance)))/Root_Deviance)  # Assess proportion information gain of split

        x <- ifelse(dim(lf$coefficients)[1] >= 2 && length(lf$coefficients[2, 4]) > 0, lf$coefficients[2, 4], NaN)
        y <- ifelse(dim(rf$coefficients)[1] >= 2 && length(rf$coefficients[2, 4]) > 0, rf$coefficients[2, 4], NaN)
        contFeature$Eval[i,10] <- ifelse(is.nan(x) & is.nan(y), NaN,
                                         ifelse(is.nan(x), y,
                                                ifelse(is.nan(y), x,mean(x, y))))
        # Assess average p-value of split
      }
    }
  }
  # Case of Categorical X Feature
  else {
    if (length(levels(Xm)) == 2) {
      contFeature$Xtype <- "Categorical"
      l     <- levels(Xm)
      contFeature$Eval <- matrix(nrow = length(l)-1, ncol = 12)         # Create matrix to store evaluations on each unique m value
      colnames(contFeature$Eval) <- c("Zk","Ztype","Xm","Xtype","Split","Eval","SplitEntropy", "Inf_Gain", "Inf_GainProp", "P_Value", "Rank", "Scale")
      contFeature$Eval[,1] <- contFeature$Zk
      contFeature$Eval[,2] <- contFeature$Ztype
      contFeature$Eval[,3] <- colnames(X)[m]
      contFeature$Eval[,4] <- contFeature$Xtype
      contFeature$Eval[,5] <- paste(colnames(X)[m],l[1],sep=" in ",collapse=NULL)   # Set matrix row names to unique m values
      contFeature$Eval[,6] <- "Deviance"
      #For all levels i of feature m, where m has 2 levels:
      fit                <- lm(data[,(znum)] ~ data[,(xnum)], data)
      fsum               <- summary(fit)
      SplitDeviance      <- deviance(fit)
      contFeature$Eval[,7] <- mean(((nrow(data[which(data[xnum] != l[1]),c(1:ncol(data))])/nrow(data))*SplitDeviance),((nrow(data[which(data[xnum]== l[1]),c(1:ncol(data))])/nrow(data))*SplitDeviance))                                 # Assess entropy of split
      contFeature$Eval[,8] <- (Root_Deviance - ((nrow(data[which(data[xnum] != l[1]),c(1:ncol(data))])/nrow(data))*SplitDeviance))                           # Assess information gain of split
      contFeature$Eval[,9] <- ((Root_Deviance - ((nrow(data[which(data[xnum] != l[1]),c(1:ncol(data))])/nrow(data))*SplitDeviance))/Root_Deviance)           # Assess proportion information gain of split
      contFeature$Eval[,10] <- (fsum$coefficients[2,4])                                 # Assess average p-value of split
    }
    else if (length(levels(Xm)) > 2) {
      contFeature$Xtype <- "Categorical"
      l     <- levels(Xm)
      contFeature$Eval <- matrix(nrow = length(l), ncol = 12)         # Create matrix to store evaluations on each unique m value
      colnames(contFeature$Eval) <- c("Zk","Ztype","Xm","Xtype","Split","Eval","SplitEntropy", "Inf_Gain", "Inf_GainProp", "P_Value", "Rank", "Scale")
      contFeature$Eval[,1] <- contFeature$Zk
      contFeature$Eval[,2] <- contFeature$Ztype
      contFeature$Eval[,3] <- colnames(X)[m]
      contFeature$Eval[,4] <- contFeature$Xtype
      contFeature$Eval[,5] <- paste(colnames(X)[m],l,sep=" in ",collapse=NULL)   # Set matrix row names to unique m values
      contFeature$Eval[,6] <- "Deviance"
      #For all levels i of feature m, where m has 2 levels:
      for (i in 1:length(levels(Xm))) {
        split <- l[i]                                                                     # Extract threshold value i
        df_right           <- data[which(data[xnum] != split),c(1:ncol(data))]   # Assess right child node
        right_fit          <- lm(df_right[,(znum)] ~ df_right[,(xnum)], df_right)
        rf                 <- summary(right_fit)
        Right_Deviance     <- deviance(right_fit)
        contFeature$Eval[i,7] <- mean(((nrow(data[which(data[xnum] != split),c(1:ncol(data))])/nrow(data))*Right_Deviance),((nrow(data[which(data[xnum]==split),c(1:ncol(data))])/nrow(data))*Right_Deviance))                                 # Assess entropy of split
        contFeature$Eval[i,8] <- (Root_Deviance - ((nrow(df_right)/nrow(data))*Right_Deviance))                  # Assess information gain of split
        contFeature$Eval[i,9] <- ((Root_Deviance - ((nrow(df_right)/nrow(data))*Right_Deviance))/Root_Deviance)  # Assess proportion information gain of split
        contFeature$Eval[i,10] <- (rf$coefficients[2,4])                              # Assess average p-value of split
      }
    }
  }
  # Define feature/outcome ranking will be performed based on information gain
  if (!is.na(contFeature$Eval[1])) {
    forRank <- as.numeric(as.character(contFeature$Eval[,8]))
    contFeature$Eval[,11] <- (rank(-forRank))                                             # Define raw ranking (kRank)
    forScale <- as.numeric(as.character(contFeature$Eval[,9]))
    contFeature$Eval[,12] <- scale(forScale)         # Apply standard normal scaling N(0,1) to proportion information gain
  }
  else {
    contFeature$Eval <- NA
  }
  # Output standardized information gain for SplitVars on k to matrix of splitvars rows, and k columns
  return(contFeature$Eval)
}

#' @usage Count_Eval(X, m, Z, k, data, continuous,
#'            quantseq = seq(0,1,0.05), splitmin=10)
#' @export
#' @aliases Count_Eval
#' @rdname MultiEval
Count_Eval <- function(X, m, Z, k, data, continuous, quantseq = seq(0,1,0.05), splitmin=10) {
  contFeature <- list("Zk"=colnames(Z)[k], "Ztype"="Event Rate", "Xm"=colnames(X)[m], "Xtype"=NA, "Eval"=NA)
  znum <- which( colnames(data)==names(Z[k]) )                           # Extract column number of outcome k
  xnum <- which( colnames(data)==names(X[m]) )                           # Extract column number of feature m
  Xm <- data[complete.cases(data[xnum]), c(xnum)]                        # Extract feature m from the sample or subsample and eliminate missing values
  Xtype <- class(Xm)                                                     # Extract feature attributes
  # Find root entropy
  df            <- data[complete.cases(data[,c(xnum,znum)]), c(1:ncol(data))]
  m1            <- glm(df[,(znum)] ~ 1, family="poisson")
  Root_Deviance <- abs(mean(residuals(m1, type = "deviance")))
  # Case of Continuous X Feature
  if (Xtype == "numeric" || Xtype =="integer") {
    contFeature$Xtype <- "Continuous"
    l <- sort(unique(Xm) )                                                  # Extract sequence of unique m values (split threshold candidates)
    # Check that split candidates result in partitions with sufficient opservations (n >= 30)
    valcheck <- data.frame(matrix(nrow = length(l), ncol = 3))
    valcheck[,1] <- l
    for (c in 1:length(l)) {
      valcheck[c,2] <- sum(table(Xm[which(Xm < l[c])]))
      valcheck[c,3] <- sum(table(Xm[which(Xm >= l[c])]))
    }
    ln <- valcheck[which(valcheck[,2] >= splitmin & valcheck[,3] >= splitmin),1]  # Define new split candidate list based on above check
    # If there are less than or equal to 100 split candidates, assess all
    if ((length(ln) <= 100 & length(ln) >0) || (length(ln) >0 & continuous == "all")) {
      ln <- ln
      contFeature$Eval <- matrix(nrow = length(ln), ncol = 12)                   # Create matrix to store evaluations on each unique m value
      colnames(contFeature$Eval) <- c("Zk","Ztype","Xm","Xtype","Split","Eval","SplitEntropy", "Inf_Gain", "Inf_GainProp", "P_Value", "Rank", "Scale")
      contFeature$Eval[,1] <- contFeature$Zk
      contFeature$Eval[,2] <- contFeature$Ztype
      contFeature$Eval[,3] <- colnames(X)[m]
      contFeature$Eval[,4] <- contFeature$Xtype
      contFeature$Eval[,5] <- paste(colnames(X)[m],ln,sep=" < ",collapse=NULL)   # Set matrix row names to unique m values
      contFeature$Eval[,6] <- "LR_Test"
      for (i in 1:length(ln)) {
        split     <- ln[i]                                                                                 # Extract threshold value i
        df_left            <- data[which(data[xnum]<split),c(1:ncol(data))]                   # Assess left child node
        if (nrow(df_left) > 0) {
          Left_m1            <- glm(df_left[,(znum)] ~ 1, family="poisson")
          Left_Deviance      <- abs(mean(residuals(Left_m1, type = "deviance")))
        }
        else {
          Left_Deviance <- NA
        }
        df_right           <- data[which(data[xnum]>=split),c(1:ncol(data))]                  # Assess right child node
        if (nrow(df_right) > 0) {
          Right_m1            <- glm(df_right[,(znum)] ~ 1, family="poisson")
          Right_Deviance      <- abs(mean(residuals(Right_m1, type = "deviance")))
        }
        else {
          Right_Deviance <- NA
        }
        contFeature$Eval[i,7]  <- mean(((nrow(df_left)/nrow(data))*Left_Deviance), ((nrow(df_right)/nrow(data))*Right_Deviance))                                    # Assess entropy of split
        contFeature$Eval[i,8]  <- (Root_Deviance - sum(((nrow(df_left)/nrow(data))*Left_Deviance), ((nrow(df_right)/nrow(data))*Right_Deviance)))                  # Assess information gain of split
        contFeature$Eval[i,9]  <- ((Root_Deviance - sum(((nrow(df_left)/nrow(data))*Left_Deviance), ((nrow(df_right)/nrow(data))*Right_Deviance)))/Root_Deviance)  # Assess proportion information gain of split
        contFeature$Eval[i,10] <- NA               # Assess average p-value of split
      }
    }
    # If there are more than 100 split candidates, assess quantiles (default: each 5% seq(0,1,0.05))
    else if (length(ln) > 100 || (length(ln) >0 & continuous == "quantile")) {
      ln <- (quantile(ln, probs = quantseq))
      contFeature$Eval <- matrix(nrow = length(ln), ncol = 12)                   # Create matrix to store evaluations on each unique m value
      colnames(contFeature$Eval) <- c("Zk","Ztype","Xm","Xtype","Split","Eval","SplitEntropy", "Inf_Gain", "Inf_GainProp", "P_Value", "Rank", "Scale")
      contFeature$Eval[,1] <- contFeature$Zk
      contFeature$Eval[,2] <- contFeature$Ztype
      contFeature$Eval[,3] <- colnames(X)[m]
      contFeature$Eval[,4] <- contFeature$Xtype
      contFeature$Eval[,5] <- paste(colnames(X)[m],ln,sep=" < ",collapse=NULL)   # Set matrix row names to unique m values
      contFeature$Eval[,6] <- "LR_Test"
      for (i in 1:length(ln)) {
        split     <- ln[i]                                                                                 # Extract threshold value i
        df_left            <- data[which(data[xnum]<split),c(1:ncol(data))]                   # Assess left child node
        if (nrow(df_left) > 0) {
          Left_m1            <- glm(df_left[,(znum)] ~ 1, family="poisson")
          Left_Deviance      <- abs(mean(residuals(Left_m1, type = "deviance")))
        }
        else {
          Left_Deviance <- NA
        }
        df_right           <- data[which(data[xnum]>=split),c(1:ncol(data))]                  # Assess right child node
        if (nrow(df_right) > 0) {
          Right_m1            <- glm(df_right[,(znum)] ~ 1, family="poisson")
          Right_Deviance      <- abs(mean(residuals(Right_m1, type = "deviance")))
        }
        else {
          Right_Deviance <- NA
        }
        contFeature$Eval[i,7]  <- mean(((nrow(df_left)/nrow(data))*Left_Deviance), ((nrow(df_right)/nrow(data))*Right_Deviance))                                    # Assess entropy of split
        contFeature$Eval[i,8]  <- (Root_Deviance - sum(((nrow(df_left)/nrow(data))*Left_Deviance), ((nrow(df_right)/nrow(data))*Right_Deviance)))                  # Assess information gain of split
        contFeature$Eval[i,9]  <- ((Root_Deviance - sum(((nrow(df_left)/nrow(data))*Left_Deviance), ((nrow(df_right)/nrow(data))*Right_Deviance)))/Root_Deviance)  # Assess proportion information gain of split
        contFeature$Eval[i,10] <- NA               # Assess average p-value of split
      }
    }
    else {
      contFeature$Eval <- NA
    }
  }
  # Case of Categorical X Feature
  else {
    if (length(levels(Xm)) == 2) {
      contFeature$Xtype <- "Categorical"
      l <- levels(Xm)
      contFeature$Eval <- matrix(nrow = length(l)-1, ncol = 12)                     # Create matrix to store evaluations on each unique m value
      colnames(contFeature$Eval) <- c("Zk","Ztype","Xm","Xtype","Split","Eval","SplitEntropy", "Inf_Gain", "Inf_GainProp", "P_Value", "Rank", "Scale")
      contFeature$Eval[,1] <- contFeature$Zk
      contFeature$Eval[,2] <- contFeature$Ztype
      contFeature$Eval[,3] <- colnames(X)[m]
      contFeature$Eval[,4] <- contFeature$Xtype
      contFeature$Eval[,5] <- paste(colnames(X)[m],l[1],sep=" in ",collapse=NULL)   # Set matrix row names to unique m values
      contFeature$Eval[,6] <- "LR_Test"
      #For all levels i of feature m, where m has 2 levels:
      split     <- l[1]                                                                                 # Extract threshold value i
      Right_m1            <- glm(data[,(znum)] ~ 1, family="poisson")
      Right_Deviance      <- abs(mean(residuals(Right_m1, type = "deviance")))

      contFeature$Eval[,7]  <- mean(((nrow(data[which(data[xnum] != split),c(1:ncol(data))])/nrow(data))*Right_Deviance), ((nrow(data[which(data[xnum]== split),c(1:ncol(data))])/nrow(data))*Right_Deviance))
      contFeature$Eval[,8]  <- (Root_Deviance - ((nrow(data[which(data[xnum] != split),c(1:ncol(data))])/nrow(data))*Right_Deviance))                  # Assess information gain of split
      contFeature$Eval[,9]  <- ((Root_Deviance - ((nrow(data[which(data[xnum] != split),c(1:ncol(data))])/nrow(data))*Right_Deviance))/Root_Deviance)  # Assess proportion information gain of split
      contFeature$Eval[,10] <- NA                     # Assess average p-value of split
    }
    else if (length(levels(Xm)) > 2) {
      contFeature$Xtype <- "Categorical"
      l <- levels(Xm)
      contFeature$Eval <- matrix(nrow = length(l), ncol = 12)         # Create matrix to store evaluations on each unique m value
      colnames(contFeature$Eval) <- c("Zk","Ztype","Xm","Xtype","Split","Eval","SplitEntropy", "Inf_Gain", "Inf_GainProp", "P_Value", "Rank", "Scale")
      contFeature$Eval[,1] <- contFeature$Zk
      contFeature$Eval[,2] <- contFeature$Ztype
      contFeature$Eval[,3] <- colnames(X)[m]
      contFeature$Eval[,4] <- contFeature$Xtype
      contFeature$Eval[,5] <- paste(colnames(X)[m],l,sep=" in ",collapse=NULL)   # Set matrix row names to unique m values
      contFeature$Eval[,6] <- "LR_Test"
      #For all levels i of feature m, where m has 2 levels:
      for (i in 1:length(levels(Xm))) {
        split     <- l[i]                                                                                 # Extract threshold value i
        df_right           <- data[which(data[xnum] != split),c(1:ncol(data))]                  # Assess right child node
        Right_m1            <- glm(df_right[,(znum)] ~ 1, family="poisson")
        Right_Deviance      <- abs(mean(residuals(Right_m1, type = "deviance")))

        contFeature$Eval[i,7]  <- mean(((nrow(data[which(data[xnum] != split),c(1:ncol(data))])/nrow(data))*Right_Deviance),((nrow(data[which(data[xnum]== split),c(1:ncol(data))])/nrow(data))*Right_Deviance))                                   # Assess entropy of split
        contFeature$Eval[i,8]  <- (Root_Deviance - ((nrow(df_right)/nrow(data))*Right_Deviance))                  # Assess information gain of split
        contFeature$Eval[i,9]  <- ((Root_Deviance - ((nrow(df_right)/nrow(data))*Right_Deviance))/Root_Deviance)  # Assess proportion information gain of split
        contFeature$Eval[i,10] <- NA                        # Assess average p-value of split
      }
    }
  }
  if (!is.na(contFeature$Eval[1])) {
    # Define feature/outcome ranking will be performed based on information gain
    forRank <- as.numeric(as.character(contFeature$Eval[,8]))
    contFeature$Eval[,11] <- (rank(-forRank))                                             # Define raw ranking (kRank)
    forScale <- as.numeric(as.character(contFeature$Eval[,9]))
    contFeature$Eval[,12] <- scale(forScale)         # Apply standard normal scaling N(0,1) to proportion information gain
  }
  else {
    contFeature$Eval <- NA
  }
  # Output standardized information gain for SplitVars on k to matrix of splitvars rows, and k columns
  return(contFeature$Eval)
}

#' @usage Surv_Eval(X, m, Z, k, data, continuous,
#'           quantseq = seq(0,1,0.05), splitmin=10)
#' @export
#' @aliases Surv_Eval
#' @rdname MultiEval
Surv_Eval <- function(X, m, Z, k, data, continuous, quantseq = seq(0,1,0.05), splitmin=10) {
  contFeature <- list("Zk"=colnames(Z)[k], "Ztype"="Event Rate", "Xm"=colnames(X)[m], "Xtype"=NA, "Eval"=NA)
  survvar <- colnames(Z)[k]
  eventvar  <- (unlist(strsplit(survvar, "_")))[4]
  timevar   <- (unlist(strsplit(survvar, "_")))[3]
  numE <- which(colnames(data)==eventvar)
  numT <- which(colnames(data)==timevar)
  E <- data[,numE]
  Tim <- data[,numT]
  xnum <- which( colnames(data)==names(X[m]) )                           # Extract column number of feature m
  znum <- which( colnames(data)==names(Z[k]) )                           # Extract column number of outcome k
  Xm <- data[complete.cases(data[xnum]), c(xnum)]                        # Extract feature m from the sample or subsample and eliminate missing values
  Xtype <- class(Xm)                                                     # Extract feature attributes
  # Find root entropy
  df            <- data[complete.cases(data[,c(xnum,znum)]), c(1:ncol(data))]
  m1            <- survival::survreg(Surv(Tim, E) ~ 1, df, dist="exponential")
  Root_Deviance <- abs(mean(residuals(m1, type="deviance")))
  # Case of Continuous X Feature
  if (Xtype == "numeric" || Xtype =="integer") {
    contFeature$Xtype <- "Continuous"
    l <- sort(unique(Xm) )                                                  # Extract sequence of unique m values (split threshold candidates)
    # Check that split candidates result in partitions with sufficient opservations (n >= 30)
    valcheck <- data.frame(matrix(nrow = length(l), ncol = 3))
    valcheck[,1] <- l
    for (c in 1:length(l)) {
      valcheck[c,2] <- sum(table(Xm[which(Xm < l[c])]))
      valcheck[c,3] <- sum(table(Xm[which(Xm >= l[c])]))
    }
    ln <- valcheck[which(valcheck[,2] >= splitmin & valcheck[,3] >= splitmin),1]  # Define new split candidate list based on above check
    # If there are less than or equal to 100 split candidates, assess all
    if ((length(ln) <= 100 & length(ln) >0) || (length(ln) >0 & continuous == "all")) {
      ln <- ln
      contFeature$Eval <- matrix(nrow = length(ln), ncol = 12)                   # Create matrix to store evaluations on each unique m value
      colnames(contFeature$Eval) <- c("Zk","Ztype","Xm","Xtype","Split","Eval","SplitEntropy", "Inf_Gain", "Inf_GainProp", "P_Value", "Rank", "Scale")
      contFeature$Eval[,1] <- contFeature$Zk
      contFeature$Eval[,2] <- contFeature$Ztype
      contFeature$Eval[,3] <- colnames(X)[m]
      contFeature$Eval[,4] <- contFeature$Xtype
      contFeature$Eval[,5] <- paste(colnames(X)[m],ln,sep=" < ",collapse=NULL)   # Set matrix row names to unique m values
      contFeature$Eval[,6] <- "LR_Test"

      for (i in 1:length(ln)) {
        split     <- ln[i]                                                                                 # Extract threshold value i
        df_left            <- data[which(data[xnum]<split),c(1:ncol(data))]                   # Assess left child node

        if (nrow(df_left) > 0) {
          Left_m1            <- survival::survreg(Surv(df_left[,(numT)], df_left[,(numE)]) ~ 1, df_left, dist="exponential")
          Left_Deviance      <- abs(mean(residuals(Left_m1, type = "deviance")))

        }
        else {
          Left_Deviance <- NA

        }
        df_right           <- data[which(data[xnum]>=split),c(1:ncol(data))]                  # Assess right child node
        if (nrow(df_right) > 0) {
          Right_m1            <- survival::survreg(Surv(df_right[,(numT)], df_right[,(numE)]) ~ 1, df_right, dist="exponential")
          Right_Deviance      <- abs(mean(residuals(Right_m1, type = "deviance")))
          m1       <- survival::survreg(Surv(Tim, E) ~ 1, data, dist="exponential")
          abs(mean(residuals(m1, type="deviance")))
        }
        else {
          Right_Deviance <- NA
        }
        contFeature$Eval[i,7]  <- mean(((nrow(df_left)/nrow(data))*Left_Deviance), ((nrow(df_right)/nrow(data))*Right_Deviance))                                    # Assess entropy of split
        contFeature$Eval[i,8]  <- (Root_Deviance - sum(((nrow(df_left)/nrow(data))*Left_Deviance), ((nrow(df_right)/nrow(data))*Right_Deviance)))                  # Assess information gain of split
        contFeature$Eval[i,9]  <- ((Root_Deviance - sum(((nrow(df_left)/nrow(data))*Left_Deviance), ((nrow(df_right)/nrow(data))*Right_Deviance)))/Root_Deviance)  # Assess proportion information gain of split
        contFeature$Eval[i,10] <- NA                # Assess average p-value of split

      }
    }
    # If there are more than 100 split candidates, assess quantiles (default: each 5% seq(0,1,0.05))
    else if (length(ln) > 100 || (length(ln) >0 & continuous == "quantile")) {
      ln <- (quantile(ln, probs = quantseq))
      contFeature$Eval <- matrix(nrow = length(ln), ncol = 12)                   # Create matrix to store evaluations on each unique m value
      colnames(contFeature$Eval) <- c("Zk","Ztype","Xm","Xtype","Split","Eval","SplitEntropy", "Inf_Gain", "Inf_GainProp", "P_Value", "Rank", "Scale")
      contFeature$Eval[,1] <- contFeature$Zk
      contFeature$Eval[,2] <- contFeature$Ztype
      contFeature$Eval[,3] <- colnames(X)[m]
      contFeature$Eval[,4] <- contFeature$Xtype
      contFeature$Eval[,5] <- paste(colnames(X)[m],ln,sep=" < ",collapse=NULL)   # Set matrix row names to unique m values
      contFeature$Eval[,6] <- "LR_Test"
      for (i in 1:length(ln)) {
        split     <- ln[i]                                                                                 # Extract threshold value i
        df_left            <- data[which(data[xnum]<split),c(1:ncol(data))]                   # Assess left child node
        if (nrow(df_left) > 0) {
          Left_m1            <- survival::survreg(Surv(df_left[,(numT)], df_left[,(numE)]) ~ 1, df_left, dist="exponential")
          Left_Deviance      <- abs(mean(residuals(Left_m1, type = "deviance")))
        }
        else {
          Left_Deviance <- NA
        }
        df_right           <- data[which(data[xnum]>=split),c(1:ncol(data))]                  # Assess right child node
        if (nrow(df_right) > 0) {
          Right_m1            <- survival::survreg(Surv(df_right[,(numT)], df_right[,(numE)]) ~ 1, df_right, dist="exponential")
          Right_Deviance      <- abs(mean(residuals(Right_m1, type = "deviance")))
        }
        else {
          Right_Deviance <- NA
        }
        contFeature$Eval[i,7]  <- mean(((nrow(df_left)/nrow(data))*Left_Deviance), ((nrow(df_right)/nrow(data))*Right_Deviance))                                    # Assess entropy of split
        contFeature$Eval[i,8]  <- (Root_Deviance - sum(((nrow(df_left)/nrow(data))*Left_Deviance), ((nrow(df_right)/nrow(data))*Right_Deviance)))                  # Assess information gain of split
        contFeature$Eval[i,9]  <- ((Root_Deviance - sum(((nrow(df_left)/nrow(data))*Left_Deviance), ((nrow(df_right)/nrow(data))*Right_Deviance)))/Root_Deviance)  # Assess proportion information gain of split
        contFeature$Eval[i,10] <- NA             # Assess average p-value of split
      }
    }
    else {
      contFeature$Eval <- NA
    }
  }
  # Case of Categorical X Feature
  else {
    if (length(levels(Xm)) == 2) {
      contFeature$Xtype <- "Categorical"
      l <- levels(Xm)
      contFeature$Eval <- matrix(nrow = length(l)-1, ncol = 12)                     # Create matrix to store evaluations on each unique m value
      colnames(contFeature$Eval) <- c("Zk","Ztype","Xm","Xtype","Split","Eval","SplitEntropy", "Inf_Gain", "Inf_GainProp", "P_Value", "Rank", "Scale")
      contFeature$Eval[,1] <- contFeature$Zk
      contFeature$Eval[,2] <- contFeature$Ztype
      contFeature$Eval[,3] <- colnames(X)[m]
      contFeature$Eval[,4] <- contFeature$Xtype
      contFeature$Eval[,5] <- paste(colnames(X)[m],l[1],sep=" in ",collapse=NULL)   # Set matrix row names to unique m values
      contFeature$Eval[,6] <- "LR_Test"
      #For all levels i of feature m, where m has 2 levels:
      split     <- l[1]                                                                                 # Extract threshold value i
      Right_m1            <- survival::survreg(Surv(data[,(numT)], data[,(numE)]) ~ 1, data, dist="exponential")
      Right_Deviance      <- abs(mean(residuals(Right_m1, type = "deviance")))

      contFeature$Eval[,7]  <- mean(((nrow(data[which(data[xnum] != split),c(1:ncol(data))])/nrow(data))*Right_Deviance),((nrow(data[which(data[xnum]== split),c(1:ncol(data))])/nrow(data))*Right_Deviance))
      contFeature$Eval[,8]  <- (Root_Deviance - ((nrow(data[which(data[xnum] != split),c(1:ncol(data))])/nrow(data))*Right_Deviance))                  # Assess information gain of split
      contFeature$Eval[,9]  <- ((Root_Deviance - ((nrow(data[which(data[xnum] != split),c(1:ncol(data))])/nrow(data))*Right_Deviance))/Root_Deviance)  # Assess proportion information gain of split
      contFeature$Eval[,10] <- NA                        # Assess average p-value of split
    }
    else if (length(levels(Xm)) > 2) {
      contFeature$Xtype <- "Categorical"
      l <- levels(Xm)
      contFeature$Eval <- matrix(nrow = length(l), ncol = 12)         # Create matrix to store evaluations on each unique m value
      colnames(contFeature$Eval) <- c("Zk","Ztype","Xm","Xtype","Split","Eval","SplitEntropy", "Inf_Gain", "Inf_GainProp", "P_Value", "Rank", "Scale")
      contFeature$Eval[,1] <- contFeature$Zk
      contFeature$Eval[,2] <- contFeature$Ztype
      contFeature$Eval[,3] <- colnames(X)[m]
      contFeature$Eval[,4] <- contFeature$Xtype
      contFeature$Eval[,5] <- paste(colnames(X)[m],l,sep=" in ",collapse=NULL)   # Set matrix row names to unique m values
      contFeature$Eval[,6] <- "LR_Test"
      #For all levels i of feature m, where m has 2 levels:
      for (i in 1:length(levels(Xm))) {
        split     <- l[i]                                                                                 # Extract threshold value i
        df_right           <- data[which(data[xnum] != split),c(1:ncol(data))]                  # Assess right child node
        Right_m1            <- survival::survreg(Surv(df_right[,(numT)], df_right[,(numE)]) ~ 1, df_right, dist="exponential")
        Right_Deviance      <- abs(mean(residuals(Right_m1, type = "deviance")))

        contFeature$Eval[i,7]  <- mean(((nrow(data[which(data[xnum] != split),c(1:ncol(data))])/nrow(data))*Right_Deviance),((nrow(data[which(data[xnum]== split),c(1:ncol(data))])/nrow(data))*Right_Deviance))
        contFeature$Eval[i,8]  <- (Root_Deviance - ((nrow(df_right)/nrow(data))*Right_Deviance))                  # Assess information gain of split
        contFeature$Eval[i,9]  <- ((Root_Deviance - ((nrow(df_right)/nrow(data))*Right_Deviance))/Root_Deviance)  # Assess proportion information gain of split
        contFeature$Eval[i,10] <- NA                          # Assess average p-value of split
      }
    }
  }
  if (!is.na(contFeature$Eval[1])) {
    # Define feature/outcome ranking will be performed based on information gain
    forRank <- as.numeric(as.character(contFeature$Eval[,8]))
    contFeature$Eval[,11] <- (rank(-forRank))                                             # Define raw ranking (kRank)
    forScale <- as.numeric(as.character(contFeature$Eval[,9]))
    contFeature$Eval[,12] <- scale(forScale)         # Apply standard normal scaling N(0,1) to proportion information gain
  }
  else {
    contFeature$Eval <- NA
  }
  # Output standardized information gain for SplitVars on k to matrix of splitvars rows, and k columns
  return(contFeature$Eval)
}


