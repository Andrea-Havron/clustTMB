
#' Extract the log likelihood of a clustTMB model
#'
#' @param object The fitted clustTMB model
#' @importFrom stats logLik
#' @method logLik clustTMB
#' @S3method logLik clustTMB
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
#' @S3method extractAIC clustTMB
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
#' @method extractAIC clustTMB
#' @S3method extractAIC clustTMB
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
#' @method summary sdreport
#' @S3method summary sdreport
summary.sdreport <- function(object, select = c("all", "fixed", "random", "report"),
                             p.value=FALSE, ...)
{
  select <- match.arg(select, several.ok = TRUE)# *several* : e.g. c("fixed", "report")
  ## check if 'meth' (or "all") is among the 'select'ed ones :
  s.has <- function(meth) any(match(c(meth, "all"), select, nomatch=0L)) > 0L
  ans1 <- ans2 <- ans3 <- NULL
  if(s.has("fixed"))  ans1 <- summary(object$sdr, "fixed")
  if(s.has("random")) ans2 <- summary(object$sdr, "random")
  if(s.has("report")) ans3 <- summary(object$sdr, "report")
  ans <- rbind(ans1, ans2, ans3)
  
  ans
}


#' Invoke TMB's summary.report function
#'
#' @title Print brief model summary
#' @param object The fitted clustTMB model
#' @return NULL
#' @method print sdreport
#' @S3method print sdreport
print.sdrpeort <- function(object){
  print(object$sdr)
  invisible(object)
}