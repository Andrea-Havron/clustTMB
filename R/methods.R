
#' Extract the log likelihood of a clustTMB model
#'
#' @param object The fitted clustTMB model
#' @importFrom stats logLik
#' @method logLik clustTMB
#' @exportS3Method logLik clustTMB
#' @return object of class \code{logLik} with attributes
logLik.clustTMB <- function(object, ...) {
  val <- -object$opt$objective
  df <- length(object$opt$par) # fixed effects only
 
  structure(val,
            df = df,
            class = "logLik"
  )
}

#' Extract the AIC of a clustTMB model
#'
#' @param object The fitted clustTMB model
#' @param scale The scale (not used)
#' @param k Penalization parameter, defaults to 2
#' @param ... Anything else
#' 
#' @method extractAIC clustTMB
#' @exportS3Method extractAIC clustTMB
#' @return numeric value
extractAIC.clustTMB <- function(object, scale, k = 2, ...) {
  L <- logLik(object)
  edf <- attr(L, "df")
  return(c(edf, c(-2 * L + k * edf)))
}

#' Get fixed-effect coefficients
#'
#' @param object The fitted clustTMB model
#' @param complete Currently ignored
#' @param ... Currently ignored
#' @importFrom stats coef
#' 
#' @method coef clustTMB
#' @exportS3Method coef clustTMB
#' @return names numeric vector
coef.clustTMB <- function(object, complete = FALSE, ...) {
  out <- object$opt$par
  out
}


#' Invoke TMB's summary.sdreport function
#'
#' @title summary tables of model parameters
#' @param object The fitted clustTMB model
#' @param select Parameter classes to select. Can be any subset of
#' \code{"fixed"} (\eqn{\hat\theta}), \code{"random"} (\eqn{\hat u}) or
#' \code{"report"} (\eqn{\phi(\hat u,\hat\theta)}) using notation as
#' \code{\link{sdreport}}.
#' @param p.value Add column with approximate p-values
#' @param ... Not used
#' @return matrix
#' @method summary clustTMB
#' @exportS3Method summary clustTMB
#' 
#' @return numeric matrix of parameter estimate and standard errors
summary.clustTMB <- function(object, select = c("all", "fixed", "random", "report"),
                             p.value=FALSE, ...)
{
  ans <- summary(object$sdr, select, p.value, ...)
  ans
}


#' Invoke TMB's print.report function
#'
#' @title Print brief model summary
#' @param object The fitted clustTMB model
#' @return NULL
#' @method print clustTMB
#' @exportS3Method print clustTMB
#' 
#' @return numeric matrix of parameter estimate and standard errors
print.clustTMB <- function(object){
  print(object$sdr)
  invisible(object)
}