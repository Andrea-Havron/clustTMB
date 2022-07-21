## combine unary or binary operator + arguments (sugar for 'substitute')
## FIXME: would be nice to have multiple dispatch, so
## (arg,op) gave unary, (arg,arg,op) gave binary operator
makeOp <- function(x, y, op = NULL) {
  if (is.null(op) || missing(y)) { ## unary
    if (is.null(op)) {
      substitute(OP(X), list(X = x, OP = y))
    } else {
      substitute(OP(X), list(X = x, OP = op))
    }
  } else {
    substitute(OP(X, Y), list(X = x, OP = op, Y = y))
  }
}

##' list of specials -- taken from enum.R
##' @noRd
findReTrmClasses <- function() {
  names(.valid_reStruct) # modified from covstruct
}

# toLang <- function(x) parse(text=x)[[1]]
#
# ## expandGrpVar(quote(x*y))
# ## expandGrpVar(quote(x/y))
# expandGrpVar <- function (f) {
#   form <- as.formula(makeOp(f, quote(`~`)))
#   mm <- terms(form)
#   tl <- attr(mm, "term.labels")
#   switch_order <- function(x) paste(rev(unlist(strsplit(x,
#                                                         ":"))), collapse = ":")
#   if (inForm(f, quote(`/`))) {
#     tl <- unname(vapply(tl, switch_order, character(1)))
#     tl <- rev(tl)
#   }
#   res <- lapply(tl, toLang)
#   return(res)
# }
#
# ##' expand interactions/combinations of grouping variables
# ##'
# ##' Modeled after lme4:::expandSlash, by Doug Bates
# ##' @param bb a list of naked grouping variables, i.e. 1 | f
# ##' @examples
# ##' ff <- clustTMB:::(y~1+(x|f/g))
# ##' clustTMB:::expandAllGrpVar(ff)
# ##' clustTMB:::expandAllGrpVar(quote(1|(f/g)/h))
# ##' clustTMB:::expandAllGrpVar(quote(1|f/g/h))
# ##' clustTMB:::expandAllGrpVar(quote(1|f*g))
# ##' @importFrom utils head
# ##' @keywords internal
# expandAllGrpVar <- function(bb) {
#   ## Return the list of '/'-separated terms
#   if (!is.list(bb)) {
#     expandAllGrpVar(list(bb))
#   } else {
#     for (i in seq_along(bb)) {
#       esfun <- function(x) {
#         if (length(x) == 1) {
#           return(x)
#         }
#         if (length(x) == 2) {
#           ## unary operator such as diag(1|f/g)
#           ## return diag(...) + diag(...) + ...
#           return(lapply(esfun(x[[2]]),
#             makeOp,
#             y = head(x)
#           ))
#         }
#         if (length(x) == 3) {
#           ## binary operator
#           if (x[[1]] == quote(`|`)) {
#             return(lapply(expandGrpVar(x[[3]]),
#               makeOp,
#               x = x[[2]], op = quote(`|`)
#             ))
#           } else {
#             return(setNames(makeOp(esfun(x[[2]]), esfun(x[[3]]),
#               op = x[[1]]
#             ), names(x)))
#           }
#         }
#       } ## esfun def.
#       return(unlist(lapply(bb, esfun)))
#     } ## loop over bb
#   }
# }
#
# ##' test whether a formula contains a particular element?
# ##' @rdname formfuns
# ##' @examples
# ##' inForm(z~.,quote(.))
# ##' inForm(z~y,quote(.))
# ##' inForm(z~a+b+c,quote(c))
# ##' inForm(z~a+b+(d+e),quote(c))
# ##' f <- ~ a + offset(x)
# ##' f2 <- z ~ a
# ##' inForm(f,quote(offset))
# ##' inForm(f2,quote(offset))
# ##' @export
# ##' @keywords internal
# inForm <- function(form, value) {
#   if (any(sapply(form,identical,value))) return(TRUE)
#   if (all(sapply(form,length)==1)) return(FALSE)
#   return(any(vapply(form,inForm,value,FUN.VALUE=logical(1))))
# }

##' find and process random effects terms
##'
##' @param term a formula or piece of a formula
##' @param debug (logical) debug?
##' @param specials list of special terms
##' @param default.special character: special to use for parenthesized terms - i.e. random effects terms with unspecified structure
##' 1. atom (not a call or an expression): NULL
##' 2. special, i.e. foo(...) where "foo" is in specials: return term
##' 3. parenthesized term: \emph{if} the head of the head is | (i.e.
##'    it is of the form (xx|gg), then convert it to the default
##'    special type; we won't allow pathological cases like
##'    ((xx|gg))
##'
##' @importFrom glmmTMB expandAllGrpVar
##' @importFrom utils head
##' @examples
##' splitForm(quote(us(x,n=2)))
##' findbars_x(~ 1 + (1 | f) + (1 | g))
##' findbars_x(~ 1 + (1 | f) + (1 | g))
##' findbars_x(~ 1 + (1|Subject))
##' findbars_x(~ (1||Subject))
##' findbars_x(~ (1|Subject))
##' findbars_x(~ 1 + x)
##' @export
findbars_x <- function(term,
                       debug = FALSE,
                       specials = character(0),
                       default.special = "norm") { # ,
  # expand_doublevert_method = c("diag_special", "split")) {

  # expand_doublevert_method <- match.arg(expand_doublevert_method)

  ds <- if (is.null(default.special)) {
    NULL
  } else {
    eval(substitute(as.name(foo), list(foo = default.special)))
  }

  ## base function
  ## defining internally in this way makes debugging slightly
  ## harder, but (1) allows easy propagation of the top-level
  ## arguments down the recursive chain; (2) allows the top-level
  ## expandAllGrpVar() operation (which also handles cases where
  ## a naked term rather than a list is returned)

  fbx <- function(term) {
    if (is.name(term) || !is.language(term)) {
      return(NULL)
    }
    if (list(term[[1]]) %in% lapply(specials, as.name)) {
      if (debug) cat("special: ", deparse(term), "\n")
      return(term)
    }
    if (head(term) == as.name("|")) { ## found x | g
      if (debug) cat("bar term:", deparse(term), "\n")
      if (is.null(ds)) {
        return(term)
      }
      return(makeOp(term, ds))
    }
    ## TODO: functionality not available in clustTMB yet
    # if (head(term) == as.name("||")) {
    #   if (expand_doublevert_method == "diag_special") {
    #     return(makeOp(makeOp(term[[2]], term[[3]],
    #                          op = quote(`|`)),
    #                   as.name("diag")))
    #   }
    #   if (expand_doublevert_method == "split") {
    #     ## need to return *multiple* elements
    #     return(lapply(expandDoubleVert(term), fbx))
    #   }
    #   stop("unknown doublevert method ", expand_doublevert_method)
    # }
    if (head(term) == as.name("(")) { ## found (...)
      if (debug) cat("paren term:", deparse(term), "\n")
      return(fbx(term[[2]]))
    }
    stopifnot(is.call(term))
    if (length(term) == 2) {
      ## unary operator, decompose argument
      if (debug) cat("unary operator:", deparse(term[[2]]), "\n")
      return(fbx(term[[2]]))
    }
    ## binary operator, decompose both arguments
    if (debug) {
      cat(
        "binary operator:", deparse(term[[2]]), ",",
        deparse(term[[3]]), "\n"
      )
    }
    c(fbx(term[[2]]), fbx(term[[3]]))
  }

  expandAllGrpVar(fbx(term))
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
##' @importFrom glmmTMB expandAllGrpVar noSpecials
##'
##' @author Steve Walker
##' @importFrom lme4 nobars
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

  fbxx <- findbars_x(formula, debug, specials)
  formSplits <- expandAllGrpVar(fbxx)

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
    reTrmFormulas <- unlist(reTrmFormulas) # Fix me:: added for rr structure when it has n = 2, gives a list of list... quick fix
    reTrmClasses <- c(
      rep(defaultTerm, length(formSplitStan)),
      vapply(
        lapply(formSplitSpec, "[[", 1),
        as.character, rep(" ", 1)
      )
    )
  } else {
    reTrmFormulas <- reTrmAddArgs <- reTrmClasses <- NULL
  }
  fixedFormula <- noSpecials(nobars(formula))

  list(
    fixedFormula = fixedFormula,
    reTrmFormulas = reTrmFormulas,
    reTrmAddArgs = reTrmAddArgs,
    reTrmClasses = reTrmClasses
  )
}
