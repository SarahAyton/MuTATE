#' @name CVsplitEval
#' @aliases CVsplitEval
#' @title Evaluate the performance of binary classifiers
#'
#' @import stringr
#' @import Hmisc
#' @import irr
#' @import Metrics
#'
# @import glm
# @import accuracy
# @import kappa2
# @import rmse
#'
#' @description A function that evaluates the performance of binary classifiers based on data split criteria
#'
#' @param splitvar A character string specifying the split criteria (e.g. "age < 60").
#' @param X A dataframe with predictor variables.
#' @param targets A data frame containing target variables.
#' @param data A data frame containing the predictor and response variables.
#' @return A list with the performance metrics for each target variable
#'
#' @details The function takes four inputs - \code{splitvar}, \code{X}, \code{targets}
#' and \code{data} - and outputs a list with measures of the performance of binary
#' classifiers. It first splits the data into two groups based on the specified split
#' criteria, and then fits binary classifiers to the resulting data. The performance
#' of the classifiers is then evaluated using various measures, such as accuracy and
#' kappa, and is stored in the output list.
#'
#' The function first splits the data into two groups based on the split criteria
#' specified in `splitvar`. The binary classifiers are then fit to the resulting data,
#' and their performance is evaluated using various measures such as accuracy, kappa, and
#' root mean squared error. The results are stored in a list and returned as the output
#' of the function.
#'
#' Note that the function uses the `glm` function in R to fit the binary classifiers, and
#' the accuracy, kappa, and root mean squared error are calculated using the `accuracy`,
#' `kappa2`, and `rmse` functions, respectively.
#'
#' @usage CVsplitEval(splitvar, X, targets, data)
#' @seealso Other relevant R functions include `glm`, `accuracy`, `kappa2`, and `rmse`.
#'

CVsplitEval <- function(splitvar, X, targets, data) {#enter dataframe and outcomes
  MT <- as.list(1:((length(colnames(targets[[2]])))))
  names(MT) <- c(colnames(targets[[2]]))
  var       <- (unlist(strsplit(splitvar, " ")))[1]
  varnum    <- which(colnames(data)==var)
  notation  <- (unlist(strsplit(splitvar, " ")))[2]
  thresh    <- str_remove(splitvar, var)
  ifelse(notation == "in", threshold <- str_remove(thresh, " in "),
         ifelse(notation == "<", threshold <- str_remove(thresh, " < "), NA))
  #if (notation == "in") {data$binrule <- as.factor(ifelse(data[,(varnum)] == threshold, "yes", "no"))}
  #else {data$binrule <- as.factor(ifelse(data[,(varnum)] < as.numeric(threshold), "yes", "no"))}
  if (notation == "in") { data$binrule <- as.factor(ifelse(!is.na(data[,(varnum)]) & data[,(varnum)] == threshold, "yes", "no")) }
  else { data$binrule <- as.factor(ifelse(!is.na(data[,(varnum)]) & data[,(varnum)] < threshold, "yes", "no")) }
  start <- Sys.time()
  for (i in 1:length(targets$Definitions)) {
    if (length(unique((targets$Z)[i])) > 1) {
      if (targets$Definitions[i] =="Cat") {
        Z <- targets$Z[,i]
        Zname    <- colnames(targets$Z)[i]
        Zvarnum    <- which(colnames(data)==Zname)
        propZ    <- prop.table(table(Z))

        dataX <- data[which(!is.na(data[,c(colnames(data)==Zname)])), c(1:ncol(data))]
        dataX <- data.frame(data[, c(Zvarnum, ncol(data))])
        dataX <- dataX[complete.cases(dataX),]
        varlev <- paste((names(propZ[(propZ)>=max(propZ)])))
        dataX$Zvar <- as.factor(dataX[,1])

        set.seed(123)
        model_b <- glm(Zvar ~ binrule, family="binomial", data = (dataX))
        Accuracy <- accuracy(dataX[complete.cases(dataX),]$Zvar, predict(model_b))

        set1 <- dataX[complete.cases(dataX),]$Zvar
        set2 <- predict(model_b)
        set  <- cbind(set1, set2)
        Kappa    <- (kappa2(set))$value

        vardata           <- as.list(1:4)

        vardata[(1)]      <- Accuracy
        vardata[(2)]      <- NA
        vardata[(3)]      <- Kappa
        vardata[(4)]      <- NA
        names(vardata)    <- c("Accuracy","AccuracySD","Kappa","KappaSD")

        MT[[i]] <- as.list(1:4)
        MT[[i]] <- vardata
      }
      else if (targets$Definitions[i] =="Surv") {

        Z <- targets$Z[,i]
        Zname    <- colnames(targets$Z)[i]
        Zvarnum    <- which(colnames(data)==Zname)

        dataX <- data[which(!is.na(data[,c(colnames(data)==Zname)])), c(1:ncol(data))]
        dataX <- data.frame(data[, c(Zvarnum, ncol(data))])
        dataX <- dataX[complete.cases(dataX),]
        dataX$Zvar <- (dataX[,1])

        set.seed(123)

        model_b <- glm(Zvar ~ binrule, family="poisson", data = (dataX[complete.cases(dataX),]))

        RMSE <- rmse(dataX[complete.cases(dataX),]$Zvar, predict(model_b))
        MAE  <- mae(dataX[complete.cases(dataX),]$Zvar, predict(model_b))
        Rsquared <- 1 - model_b$deviance/model_b$null.deviance

        vardata           <- as.list(1:6)

        vardata[(1)]      <- RMSE
        vardata[2]        <- NA
        vardata[3]        <- Rsquared
        vardata[4]        <- NA
        vardata[5]        <- MAE
        vardata[6]        <- NA
        names(vardata)    <- c("RMSE","RMSESD","Rsquared","RsquaredSD","MAE","MAESD")

        MT[[i]] <- as.list(1:6)
        MT[[i]] <- vardata
      }
      else if (targets$Definitions[i] =="Count") {
        Z <- targets$Z[,i]
        Zname    <- colnames(targets$Z)[i]
        Zvarnum    <- which(colnames(data)==Zname)

        dataX <- data[which(!is.na(data[,c(colnames(data)==Zname)])), c(1:ncol(data))]
        dataX <- data.frame(data[, c(Zvarnum, ncol(data))])
        dataX <- dataX[complete.cases(dataX),]
        dataX$Zvar <- (dataX[,1])

        set.seed(123)

        model_b <- glm(Zvar ~ binrule, family="poisson", data = (dataX[complete.cases(dataX),]))

        RMSE <- rmse(dataX[complete.cases(dataX),]$Zvar, predict(model_b))
        MAE  <- mae(dataX[complete.cases(dataX),]$Zvar, predict(model_b))
        Rsquared <- 1 - model_b$deviance/model_b$null.deviance

        vardata           <- as.list(1:6)

        vardata[(1)]      <- RMSE
        vardata[2]        <- NA
        vardata[3]        <- Rsquared
        vardata[4]        <- NA
        vardata[5]        <- MAE
        vardata[6]        <- NA
        names(vardata)    <- c("RMSE","RMSESD","Rsquared","RsquaredSD","MAE","MAESD")

        MT[[i]] <- as.list(1:6)
        MT[[i]] <- vardata
      }
      else {
        Zname    <- colnames(targets$Z)[i]

        dataX <- data[which(!is.na(data[,c(colnames(data)==Zname)])), c(1:ncol(data))]
        dataX <- dataX[ , colSums(is.na(dataX)) == 0]
        ZnumX <- which(colnames(dataX)==Zname)
        dataX$Zvar <- dataX[,c(ZnumX)]

        set.seed(123)

        model_b <- glm(Zvar ~ binrule, data = (dataX[complete.cases(dataX),]), family = gaussian)

        RMSE <- rmse(dataX[complete.cases(dataX),]$Zvar, predict(model_b))
        MAE  <- mae(dataX[complete.cases(dataX),]$Zvar, predict(model_b))
        Rsquared <- 1 - model_b$deviance/model_b$null.deviance

        vardata           <- as.list(1:6)

        vardata[(1)]      <- RMSE
        vardata[2]        <- NA
        vardata[3]        <- Rsquared
        vardata[4]        <- NA
        vardata[5]        <- MAE
        vardata[6]        <- NA
        names(vardata)    <- c("RMSE","RMSESD","Rsquared","RsquaredSD","MAE","MAESD")

        MT[[i]] <- as.list(1:6)
        MT[[i]] <- vardata
      }
    }
    else {
      if (targets$Definitions[i] =="Cat") {

        set.seed(123)
        Accuracy <- 0
        Kappa    <- 0

        vardata           <- as.list(1:4)

        vardata[(1)]      <- Accuracy
        vardata[(2)]      <- NA
        vardata[(3)]      <- Kappa
        vardata[(4)]      <- NA
        names(vardata)    <- c("Accuracy","AccuracySD","Kappa","KappaSD")

        MT[[i]] <- as.list(1:4)
        MT[[i]] <- vardata
      }
      else if (targets$Definitions[i] =="Surv") {

        set.seed(123)

        RMSE <- 0
        MAE  <- NA
        Rsquared <- 0

        vardata           <- as.list(1:6)

        vardata[(1)]      <- RMSE
        vardata[2]        <- NA
        vardata[3]        <- Rsquared
        vardata[4]        <- NA
        vardata[5]        <- MAE
        vardata[6]        <- NA
        names(vardata)    <- c("RMSE","RMSESD","Rsquared","RsquaredSD","MAE","MAESD")

        MT[[i]] <- as.list(1:6)
        MT[[i]] <- vardata
      }
      else if (targets$Definitions[i] =="Count") {

        set.seed(123)

        RMSE <- 0
        MAE  <- NA
        Rsquared <- 0

        vardata           <- as.list(1:6)

        vardata[(1)]      <- RMSE
        vardata[2]        <- NA
        vardata[3]        <- Rsquared
        vardata[4]        <- NA
        vardata[5]        <- MAE
        vardata[6]        <- NA
        names(vardata)    <- c("RMSE","RMSESD","Rsquared","RsquaredSD","MAE","MAESD")

        MT[[i]] <- as.list(1:6)
        MT[[i]] <- vardata
      }
      else {
        set.seed(123)

        RMSE <- 0
        MAE  <- NA
        Rsquared <- 0

        vardata           <- as.list(1:6)

        vardata[(1)]      <- RMSE
        vardata[2]        <- NA
        vardata[3]        <- Rsquared
        vardata[4]        <- NA
        vardata[5]        <- MAE
        vardata[6]        <- NA
        names(vardata)    <- c("RMSE","RMSESD","Rsquared","RsquaredSD","MAE","MAESD")

        MT[[i]] <- as.list(1:6)
        MT[[i]] <- vardata
      }
    }
  }
  end <- Sys.time()
  difftime(end, start)
  return(MT)
}
