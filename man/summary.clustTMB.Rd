% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{summary.clustTMB}
\alias{summary.clustTMB}
\title{summary tables of model parameters}
\usage{
\method{summary}{clustTMB}(
  object,
  select = c("all", "fixed", "random", "report"),
  p.value = FALSE,
  ...
)
}
\arguments{
\item{object}{The fitted clustTMB model}

\item{select}{Parameter classes to select. Can be any subset of
\code{"fixed"} (\eqn{\hat{\theta}}), \code{"random"} (\eqn{\hat{u}}) or
\code{"report"} (\eqn{\phi(\hat{u},\hat{\theta)}}) using notation as
\code{\link[TMB:sdreport]{TMB::sdreport()}}.}

\item{p.value}{Add column with approximate p-values}

\item{...}{Currently ignored}
}
\value{
numeric matrix of parameter estimate and standard errors
}
\description{
Invoke TMB's summary.sdreport function
}
