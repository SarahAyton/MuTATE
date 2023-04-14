#' @name SplitPrep
#' @title Prepare data for splitting in a decision tree
#'
#' @description \code{SplitPrep()} is a function that prepares data frames for use in tree-based
#' models by identifying and removing candidate predictor variables that are unsuitable for
#' splitting at a given node. This function takes in three data frames, \code{Xdf}, \code{df},
#' and \code{data_splt} and a character string \code{parentsplit}.
#' It returns a list of two data frames: \code{splitX_df}, which contains the predictor variables
#' that can be used for splitting at a given node, and \code{split_data_df}, which contains the response
#' variables and any other variables that should be included in the data frame for further analysis.
#'
#' @details \code{SplitPrep()} takes in two data frames, one containing the predictor variables
#' and another containing the response variable. It identifies which predictor
#' variables to use for splitting in the decision tree and returns a list
#' containing two data frames: the modified predictor variable data frame and the
#' modified response variable data frame.
#'
#' The function first drops any unused factor levels in both input data frames
#' (\code{Xdf} and \code{df}) using the \code{droplevels()} function. It then loops through the predictor
#' variables in \code{Xdf} to determine which ones should be kept for splitting at the current node.
#'
#' For each predictor variable, the function first checks if it was used for splitting at the
#' parent node (identified by the parentsplit argument). If so, it is added to a list of
#' variables to drop (\code{droplist}).
#'
#' If the predictor variable is categorical (i.e., a factor variable), the function checks if
#' it has less than 2 levels or a prevalence of less than 5% in the current node. If either of
#' these conditions is true, the variable is added to \code{droplist}.
#'
#' If the predictor variable is continuous, the function checks if it has only one observed
#' value. If so, it is added to \code{droplist}.
#'
#' Finally, if \code{droplist} exists, the function removes the corresponding columns from \code{splitX_df}
#' and renames the remaining columns to match the original column names. The function then
#' returns a list containing \code{splitX_df} and \code{split_data_df}.
#'
#' This function is used to pre-process data before fitting tree-based models to ensure
#' that only suitable predictor variables are used for splitting.
#'
#' @param Xdf A data frame containing the predictor variables.
#' @param df A data frame containing the response variable.
#' @param data_splt The data used for the split.
#' @param parentsplit The parent split.
#'
#' @return A list containing two data frames: the modified predictor variable
#' data frame and the modified response variable data frame.
#' @export

SplitPrep <- function(Xdf, df, data_splt, parentsplit) {
  df <- droplevels(data.frame(df))
  splitX_df     <- droplevels(data.frame(Xdf), c(1:ncol(data.frame(Xdf))))
  split_data_df <- droplevels(data.frame(df),  c(1:ncol(data.frame(df))))
  for (n in 1:ncol(data.frame(splitX_df))) {
    # If the X-var was used in the parent split:
    if (!is.null(data.frame(splitX_df)[,c(n)]) & colnames(data.frame(splitX_df))[(n)] == parentsplit) {
      # if the droplist does exist, append additional variables with only one level
      if (exists("droplist")){
        droplist <- c(droplist, n)
      }
      # if the droplist doesn't exist, create it
      else {
        droplist <- c(n)
      }
    }
    # If the X-var is categorical:
    else if (!is.null(data.frame(splitX_df)[,c(n)]) & as.character(class(data.frame(splitX_df)[,c(n)])) == "factor") {
      # If the X-var has < 2 levels, remove it from split candidates
      if (length(levels(data.frame(splitX_df)[,c(n)])) < 2) {
        # if the droplist does exist, append additional variables with only one level
        if (exists("droplist")){
          droplist <- c(droplist, n)
        }
        # if the droplist doesn't exist, create it
        else {
          droplist <- c(n)
        }
      }
      # If the X-var has a prevalence in < 5% of observations in that node, remove it from split candidates
      else if (!is.null(table(data.frame(splitX_df)[,c(n)])) & (min(table(data.frame(splitX_df)[,c(n)])) < (0.05*nrow(data_splt)))) {
        # if the droplist does exist, append additional variables with less than 5% level prevalence
        if (exists("droplist")){
          droplist <- c(droplist, n)
        }
        # if the droplist doesn't exist, create it
        else {
          droplist <- c(n)
        }
      }
      # Otherwise, keep all split candidates
      else {
        # if the droplist does exist, append additional variables with less than 5% level prevalence
        if (exists("droplist")){
          droplist <- (droplist)
        }
        # if the droplist doesn't exist, create it
        else {
          droplist <- (NA)
          rm(droplist)
        }
      }
    }
    # If the X-var is continuous: only keep those with at least two observed values as split candidates
    else {
      if (!is.null(data.frame(splitX_df)[,c(n)]) & length(unique(data.frame(splitX_df)[,c(n)])) < 2) {
        # if the droplist does exist, append additional variables with only one value
        if (exists("droplist")){
          droplist <- c(droplist, n)
        }
        # if the droplist doesn't exist, create it
        else {
          droplist <- c(n)
        }
      }
      else {
        # if the droplist does exist, append additional variables with less than 5% level prevalence
        if (exists("droplist")){
          droplist <- (droplist)
        }
        # if the droplist doesn't exist, create it
        else {
          droplist <- (NA)
          rm(droplist)
        }
      }
    }
  }
  remnames <- names(splitX_df)
  if (exists("droplist")) {
    splitX_df  <- (splitX_df[,-(droplist)])
    splitX_df <- data.frame(splitX_df)
    colnames(splitX_df) <- remnames[-(droplist)]
  }
  else {
    splitX_df  <- splitX_df
  }
  return(list(splitX_df, split_data_df))
}
