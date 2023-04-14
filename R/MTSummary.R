#' @name MTSummary
#' @aliases MTSummary
#' @title Summarize multiple outcome variables
#'
#' @description This function takes in a data frame of outcomes and a target table with definitions
#' for each target, and outputs summary statistics for each target. It currently
#' supports continuous, categorical, survival, and count outcome variable types.

#' @param targets A data frame of targets with two columns, Definitions and Z. Definitions are defined as character
#' vectors indicating the type of variable: "Cont" for continuous, "Cat" for categorical, "Surv" for survival, and "Count" for count data.
#' Z is a matrix of outcomes with each row representing a patient and each column representing a target.
#' @param data A data frame of patient-level data with columns matching the names in targets.
#'
#' @details The `MTSummary()` function computes a summary of the outcomes in a dataset for multiple
#' types of variables, including categorical, survival and count variables. The targets
#' argument is a list containing a vector of outcome definitions and a matrix or data
#' frame of outcomes. The data argument is a matrix or data frame with the data to be summarized.
#'
#' For continuous targets, the function computes the mean, median, inter-quartile range, and MSE.
#'
#' For categorical targets, the function computes the counts and proportions of each
#' category, the most frequent category and its proportion, and the expected loss (i.e.,
#' proportion of observations that do not belong to the most frequent category).
#'
#' For survival targets, the function computes the number of missing values, the number
#' of events, the estimated rate of the event, and the deviance.
#'
#' For count targets, the function computes the number of missing values, the total count,
#' the estimated rate of the count, and the deviance.
#'
#' @return A list of lists with each element corresponding to an outcome in targets.
#' @export
#'
#' @importFrom stats complete.cases deviance gaussian glm lm pnorm predict quantile residuals weighted.mean
#' @importFrom utils head

MTSummary <- function(targets, data) {#enter dataframe and outcomes
  MT <- as.list(1:((length(colnames(targets[[2]])))))
  names(MT) <- c(colnames(targets[[2]]))
  for (i in 1:length(targets$Definitions)) {
    if (targets$Definitions[i] =="Cat") {
      Z <- targets$Z[,i]
      Zname    <- colnames(targets$Z)[i]
      num <- which(colnames(data)==Zname)

      missingZ <- sum(is.na (data[which(colnames(data)==Zname)]))
      countsZ  <- table(Z)
      propZ    <- prop.table(countsZ)

      vardata           <- (as.list(1:(3+2*(length(levels(Z))))))
      vardata_names     <- as.list(1:(3+2*(length(levels(Z)))))
      vardata_names[1]  <- c("absent")
      vardata[(1)]    <- ((missingZ))
      for (j in 1:length(levels(as.factor(Z)))) {
        vardata[(j*2)]          <- (countsZ[[j]])
        vardata[(j*2+1)]        <- (round(propZ[[j]], 4))
        vardata_names[(j*2)]   <- (paste(levels(Z)[j],"count",sep = "_"))
        vardata_names[(j*2+1)] <- (paste(levels(Z)[j],"prop",sep = "_"))
      }
      vardata[(2+2*(length(levels(Z))))] <- paste(names(propZ[(propZ)>=max(propZ)]))
      vardata_names[(2+2*(length(levels(Z))))]  <- c("predicted")
      vardata[(3+2*(length(levels(Z))))] <- (round(1-max(propZ),4))
      vardata_names[(3+2*(length(levels(Z))))]  <- c("expectedloss")
      names(vardata)  <- vardata_names

      MT[[i]] <- as.list(1:(3+2*(length(levels(as.factor(Z))))))
      MT[[i]] <- vardata
    }
    else if (targets$Definitions[i] =="Surv") {
      survvar <- colnames(targets$Z)[i]
      Zname   <- (unlist(strsplit(survvar, "_")))[4]
      eventvar  <- (unlist(strsplit(survvar, "_")))[4]
      timevar   <- (unlist(strsplit(survvar, "_")))[3]
      num <- which(colnames(data)==Zname)
      numE <- which(colnames(data)==eventvar)
      numT <- which(colnames(data)==timevar)
      Z <- as.numeric(as.character(data[,num]))
      E <- data[,numE]
      Tim <- data[,numT]

      missingZ <- sum(is.na (Z))
      sumcheck <- sum(Z)
      countsZ  <- table(Z)
      propZ    <- prop.table(countsZ)

      m1       <- survreg(Surv(Tim, E) ~ 1, data, dist="exponential")

      vardata           <- as.list(1:4)
      vardata_names     <- as.list(1:4)
      vardata_names[1]  <- c("absent")
      vardata[(1)]      <- ((missingZ))
      vardata_names[2]  <- c("events_count")
      vardata[2]        <- ifelse(sumcheck > 0 & length(countsZ) > 1, (countsZ[[2]]), sum(Z))
      vardata_names[3]  <- c("estimatedrate")
      vardata[3]        <- ifelse(sumcheck > 0 & length(propZ) > 1, (propZ[[2]]), 1)
      vardata_names[4]  <- c("deviance")
      vardata[4]        <- abs(mean(residuals(m1, type="deviance")))
      names(vardata)    <- vardata_names

      MT[[i]] <- as.list(1:4)
      MT[[i]] <- vardata
    }
    else if (targets$Definitions[i] =="Count") {
      Zname    <- colnames(targets$Z)[i]
      num      <- which(colnames(data)==Zname)
      Z <- data[,num]

      missingZ <- sum(is.na (data[which(colnames(data)==Zname)]))
      countsZ  <- sum(Z)
      propZ    <- mean(Z)
      m1       <- glm(Z ~ 1, family="poisson", data)


      vardata           <- as.list(1:4)
      vardata_names     <- as.list(1:4)
      vardata_names[1]  <- c("absent")
      vardata[(1)]      <- ((missingZ))
      vardata_names[2]  <- c("events_count")
      vardata[2]        <- countsZ
      vardata_names[3]  <- c("estimatedrate")
      vardata[3]        <- propZ
      vardata_names[4]  <- c("deviance")
      vardata[4]        <- abs(mean(residuals(m1, type = "deviance")))
      names(vardata)    <- vardata_names

      MT[[i]] <- as.list(1:4)
      MT[[i]] <- vardata
    }
    else {
      Z <- targets$Z[,i]
      Zname    <- colnames(targets$Z)[i]
      num      <- which(colnames(data)==Zname)

      missingZ <- sum(is.na (Z))
      summaryZ <- summary(Z)
      model <- lm(Z~1, data)
      model_summ <-summary(model)

      vardata           <- as.list(1:6)
      vardata_names     <- as.list(1:6)
      vardata_names[1]  <- c("absent")
      vardata[(1)]      <- ((missingZ))
      vardata_names[2]  <- c("mean")
      vardata[2]        <- round(summaryZ[[4]], 4)
      vardata_names[3]  <- c("median")
      vardata[3]        <- round(summaryZ[[3]], 4)
      vardata_names[4]  <- c("IQR25")
      vardata[4]        <- round(summaryZ[[2]],4)
      vardata_names[5]  <- c("IQR75")
      vardata[5]        <- round(summaryZ[[5]],4)
      vardata_names[6]  <- c("MSE")
      vardata[6]        <- round((mean(model_summ$residuals^2)),4)
      names(vardata)    <- vardata_names

      MT[[i]] <- as.list(1:6)
      MT[[i]] <- vardata
    }
  }
  return(MT)
}
