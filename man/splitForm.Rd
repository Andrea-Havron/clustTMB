% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lme4_utils.R
\name{splitForm}
\alias{splitForm}
\title{Split formula containing special random effect terms}
\usage{
splitForm(
  formula,
  defaultTerm = "norm",
  allowFixedOnly = TRUE,
  allowNoSpecials = TRUE,
  debug = FALSE
)
}
\arguments{
\item{formula}{a formula containing special random effect terms}

\item{defaultTerm}{default type for non-special RE terms}

\item{allowFixedOnly}{(logical) are formulas with no RE terms OK?}

\item{allowNoSpecials}{(logical) are formulas with only standard RE terms OK?}

\item{debug}{(logical) debug?}
}
\value{
a list containing elements \code{fixedFormula};
\code{reTrmFormulas} list of \code{x | g} formulas for each term;
\code{reTrmAddArgs} list of function+additional arguments, i.e. \code{list()} (non-special), \code{foo()} (no additional arguments), \code{foo(addArgs)} (additional arguments); \code{reTrmClasses} (vector of special functions/classes, as character)
}
\description{
Parse a formula into fixed formula and random effect terms,
treating 'special' terms appropriately
}
\details{
Taken from Steve Walker's lme4ord,
ultimately from the flexLambda branch of lme4
\url{https://github.com/stevencarlislewalker/lme4ord/blob/master/R/formulaParsing.R}.  Mostly for internal use.
}
\examples{
splitForm(~x+y)                     ## no specials or RE
splitForm(~x+y+(f|g))               ## no specials
splitForm(~x+y+diag(f|g))           ## one special
splitForm(~x+y+(diag(f|g)))         ## 'hidden' special
splitForm(~x+y+(f|g)+cs(1|g))       ## combination
splitForm(~x+y+(1|f/g))             ## 'slash'; term
splitForm(~x+y+(1|f/g/h))             ## 'slash'; term
splitForm(~x+y+(1|(f/g)/h))             ## 'slash'; term
splitForm(~x+y+(f|g)+cs(1|g)+cs(a|b,stuff))  ## complex special
splitForm(~(((x+y))))               ## lots of parentheses
splitForm(~1+rr(f|g,n=2))

}
\author{
Steve Walker
}
