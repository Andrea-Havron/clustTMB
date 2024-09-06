##' list of specials -- taken from enum.R
##' @noRd
findReTrmClasses <- function() {
  names(.valid_reStruct) # modified from covstruct
}



##' Parse a formula into fixed formula and random effect terms,
##' treating 'special' terms appropriately
##'
##' Taken from Steve Walker's lme4ord,
##' ultimately from the flexLambda branch of lme4
##' <https://github.com/stevencarlislewalker/lme4ord/blob/master/R/formulaParsing.R>.  Mostly for internal use.
##' @title Split formula containing special random effect terms
##' @param formula a formula containing special random effect terms
##' @param defaultTerm default type for non-special RE terms
##' @param allowFixedOnly (logical) are formulas with no RE terms OK?
##' @param allowNoSpecials (logical) are formulas with only standard RE terms OK?
##' @param debug (logical) debug?
##' @return a list containing elements \code{fixedFormula};
##' \code{reTrmFormulas} list of \code{x | g} formulas for each term;
##' \code{reTrmAddArgs} list of function+additional arguments, i.e. \code{list()} (non-special), \code{foo()} (no additional arguments), \code{foo(addArgs)} (additional arguments); \code{reTrmClasses} (vector of special functions/classes, as character)
##' @examples
##' splitForm(~x+y)                     ## no specials or RE
##' splitForm(~x+y+(f|g))               ## no specials
##' splitForm(~x+y+diag(f|g))           ## one special
##' splitForm(~x+y+(diag(f|g)))         ## 'hidden' special
##' splitForm(~x+y+(f|g)+cs(1|g))       ## combination
##' splitForm(~x+y+(1|f/g))             ## 'slash'; term
##' splitForm(~x+y+(1|f/g/h))             ## 'slash'; term
##' splitForm(~x+y+(1|(f/g)/h))             ## 'slash'; term
##' splitForm(~x+y+(f|g)+cs(1|g)+cs(a|b,stuff))  ## complex special
##' splitForm(~(((x+y))))               ## lots of parentheses
##' splitForm(~1+rr(f|g,n=2))
##'
##' @importFrom reformulas makeOp
##'
##' @author Steve Walker
##' @export
splitForm <- function(formula,
                      defaultTerm = "norm",
                      allowFixedOnly = TRUE,
                      allowNoSpecials = TRUE,
                      debug = FALSE) {
  ## logic:

  ## string for error message *if* specials not allowed
  ## (probably package-specific)
  noSpecialsAlt <- "lmer or glmer"

  specials <- findReTrmClasses()

  ## formula <- expandDoubleVerts(formula)
  ## split formula into separate
  ## random effects terms
  ## (including special terms)

  fbxx <- reformulas::findbars_x(formula, debug, specials)
  formSplits <- glmmTMB::expandAllGrpVar(fbxx)

  if (length(formSplits) > 0) {
    formSplitID <- vapply(
      lapply(formSplits, "[[", 1),
      as.character, rep(" ", 1)
    )
    # warn about terms without a
    # setReTrm method

    parenTerm <- formSplitID == "("
    # capture additional arguments
    reTrmAddArgs <- lapply(formSplits, "[", -2)[!parenTerm]
    # remove these additional
    # arguments
    formSplits <- lapply(formSplits, "[", 1:2)
    # standard RE terms
    formSplitStan <- formSplits[parenTerm]
    # structured RE terms
    formSplitSpec <- formSplits[!parenTerm]

    if (!allowNoSpecials) {
      if (length(formSplitSpec) == 0) {
        stop(
          "no special covariance structures. ",
          "please use ", noSpecialsAlt,
          " or use findReTrmClasses() for available structures."
        )
      }
    }

    reTrmFormulas <- c(
      lapply(formSplitStan, "[[", 2),
      lapply(formSplitSpec, "[[", 2)
    )
    reTrmFormulas <- unlist(reTrmFormulas)
    reTrmClasses <- c(
      # Fix me:: added for rr structure when it has n = 2,
      # gives a list of list... quick fix
      rep(defaultTerm, length(formSplitStan)),
      vapply(
        lapply(formSplitSpec, "[[", 1),
        as.character, rep(" ", 1)
      )
    )
  } else {
    reTrmFormulas <- reTrmAddArgs <- reTrmClasses <- NULL
  }
  fixedFormula <- reformulas::noSpecials(lme4::nobars(formula))

  list(
    fixedFormula = fixedFormula,
    reTrmFormulas = reTrmFormulas,
    reTrmAddArgs = reTrmAddArgs,
    reTrmClasses = reTrmClasses
  )
}
