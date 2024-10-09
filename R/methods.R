#' Extract the log likelihood of a clustTMB model
#'
#' @param object The fitted clustTMB model
#' @param ... Currently ignored
#' @importFrom stats logLik
#' @return object of class \code{logLik} with attributes
#' @exportS3Method logLik clustTMB
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
#' @param fit The fitted clustTMB model
#' @param scale The scale, currently ignored
#' @param k Penalization parameter, defaults to 2
#' @param ... Currently ignored
#' @return numeric value
#' @exportS3Method extractAIC clustTMB
extractAIC.clustTMB <- function(fit, scale, k = 2, ...) {
  L <- logLik(fit)
  edf <- attr(L, "df")
  return(c(edf, c(-2 * L + k * edf)))
}

#' Get fixed-effect coefficients
#'
#' @param object The fitted clustTMB model
#' @param complete Currently ignored
#' @param ... Currently ignored
#' @importFrom stats coef
#' @return names numeric vector
#' @exportS3Method coef clustTMB
coef.clustTMB <- function(object, complete = FALSE, ...) {
  out <- object$opt$par
  out
}


#' Invoke TMB's summary.sdreport function
#'
#' @title summary tables of model parameters
#' @param object The fitted clustTMB model
#' @param select Parameter classes to select. Can be any subset of
#' \code{"fixed"} (\eqn{\hat{\theta}}), \code{"random"} (\eqn{\hat{u}}) or
#' \code{"report"} (\eqn{\phi(\hat{u},\hat{\theta)}}) using notation as
#' [TMB::sdreport()].
#' @param p.value Add column with approximate p-values
#' @param ... Currently ignored
#' @return numeric matrix of parameter estimate and standard errors
#' @exportS3Method summary clustTMB
summary.clustTMB <- function(object, select = c("all", "fixed", "random", "report"),
                             p.value = FALSE, ...) {
  ans <- summary(object$sdr, select, p.value, ...)
  ans
}


#' Invoke TMB's print.report function
#'
#' @title Print brief model summary
#' @param x The fitted clustTMB model
#' @param ... Not used
#' @return numeric matrix of parameter estimate and standard errors
#' @exportS3Method print clustTMB
print.clustTMB <- function(x, ...) {
  print(x$sdr)
  invisible(x)
}
