#' Fit a finite mixture model using TMB
#'
#'
#' @param response A numeric vector, matrix, or data frame of observations. When data are multivariate, rows correspond to observations and columns correspond to the multivariate response.
#' @param expertformula Formula defining expert model. This formula corresponds to the covariates included in the response densities. Defaults to intercept only (~1) when no covariates are used.
#' @param gatingformula Formula defining gating model.  This formula corresponds to the covariates included in the mixing proportions (logistic regression). Defaults to intercept only (~1) when no covariates are used. When a random effects term is included in the gating network, this formula will be updated so that the intercept term is removed.
#' @param expertdata Data frame containing expert model covariates.
#' @param gatingdata Data frame containing gating model covariates.
#' @param family Statistical distribution and link function of observations.
#' @param Offset Constant in expertformula only used to offset density expectation.
#' @param G Integer specifying the number of clusters.
#' @param rr List specifying dimension of rank reduction in spatial, temporal, and/or random effects. Dimension must be smaller than the total dimension of the response. Rank reduction is applied only to the expertformula random effects. The rank reduction reduces the dimensionality of a correlated multivariate response to a smaller dimension independent response. When used, the covariance structure of the response is switched to 'Diagonal.' Defaults to NULL, no rank reduction. If rank reduction is used in conjunction with a random effect, that random effect must also be specified in the expert formula. Currently, rank reduction on temporal random effects is disabled.
#' @param covariance.structure A character string specifying the covariance structure of the response using mclust naming scheme. See description of modelNames under ?Mclust for details.
#' @param Start Set initial values for random effects parameters (fixed and random terms)
#' @param Map Vector indicating parameter maps, see ?[TMB::MakeADFun()] for details. Defaults in clustTMB control this map argument and user input is limited
#' @param initialization.args A list consisting of initialization settings used to generate initial values.
#' control Calls [init.options()] to generate settings for initial values. Arguments of [init.options()] can be specified by the user.
#' 1. init.method - Single character string indicating initial clustering method. Methods include: hc, quantile, random, mclust, kmeans, mixed, user. Defaults to 'hc'. In the case where data are univariate and there are no covariates in the gating/expert formula, this defaults to 'quantile'
#' 2. hc.options - Named list of two character strings specifying hc modelName and hcUse when init.method = 'hc'. The default modelName is 'VVV' and the default use is 'SVD' unless gating/expert covariates specified, in which case default in VARS. See ?mclust::mclust.options for complete list of options.
#' 3. mix.method - String stating initialization method for mixed-type data (init.method = 'mixed'). Current default when Tweedie family specified. Options include: Gower kmeans (default), Gower hclust, and kproto.
#' 4. user - Numeric or character vector defining user specified initial classification. init.method must be set to 'user' when using this option.
#' @param spatial.list List of data objects needed when fitting a spatial GMRF model
#' @param projection.dat Spatial Points class of projection coordinates or Spatial Points Dataframe containing projection coordinates and projection covariates
#' @param control List controlling whether models are run and whether standard errors are calculated.
#' @importFrom stats model.matrix model.frame as.formula setNames terms update vcov
#' @importFrom lme4 nobars
#' @importFrom utils data
#' @importFrom Matrix t
#' @useDynLib clustTMB
#'
#' @return list of objects from fitted model
#' @export
#'
#' @examples
#' data("faithful")
#' m1 <- clustTMB(response = faithful, covariance.structure = "VVV")
#' plot(faithful$eruptions, faithful$waiting, pch = 16, col = m1$report$classification + 1)
clustTMB <- function(response = NULL,
                     expertformula = ~1,
                     gatingformula = ~1,
                     expertdata = NULL,
                     gatingdata = NULL,
                     family = gaussian(link = "identity"),
                     Offset = NULL,
                     G = 2,
                     rr = list(spatial = NULL, temporal = NULL, random = NULL),
                     covariance.structure = NULL,
                     Start = list(),
                     Map = list(),
                     initialization.args = list(control = init.options()),
                     spatial.list = list(
                       loc = NULL,
                       mesh = NULL,
                       init.range = list(
                         gating.range = NULL,
                         expert.range = NULL
                       )
                     ),
                     projection.dat = NULL,
                     control = run.options()) {
  # set str of response to matrix
  response <- as.matrix(response)

  # initialize dim.list to keep track of data/parameter dimensions
  dim.list <- list(
    n.i = nrow(response),
    n.j = ncol(response),
    n.t = NULL,
    n.g = G,
    n.f.sp = NULL,
    n.f.rand = NULL,
    n.r.g = spatial.list$init.range$gating.range,
    n.r.e = spatial.list$init.range$expert.range,
    n.v = NULL
  )
  dim.list$nl.fix <- ifelse(dim.list$n.j > 1,
    (dim.list$n.j^2 - dim.list$n.j) / 2,
    1
  )

  # exception-handling
  if (is.null(covariance.structure)) {
    if (dim.list$n.j == 1) {
      covariance.structure <- "E"
      message("Setting covariance structure to univariate E")
    } else {
      covariance.structure <- "EII"
      message("Setting covariance structure to multivariate EII")
    }
  }

  if (!(covariance.structure %in% fixStruct.names())) {
    stop("covariance structure not supported")
  }

  if (gatingformula == ~1) {
    ll.method <- 0
  } else {
    ll.method <- 2
  }

  if ((covariance.structure == "E" | covariance.structure == "V") &
    dim.list$n.j > 1) {
    stop("Need to specify multivariate covariance structure")
  }
  if (!(covariance.structure == "E" | covariance.structure == "V") &
    dim.list$n.j == 1) {
    stop("Need to specify univariate covariance structure")
  }
  if (!(family[[1]] %in% names(.valid_family))) {
    stop("family not supported")
  }
  if (!(family[[2]] %in% names(.valid_link))) {
    stop("link not supported in specified family")
  }
  if (!is.null(spatial.list$loc)) {
    if (!inherits(spatial.list$loc, "sf")) {
      stop("Locations need to be defined as sf data.frames")
    }
  }

  # set up random component of model
  reOut <- mkRandom(
    expertformula,
    gatingformula,
    expertdata,
    gatingdata,
    spatial.list,
    dim.list
  )
  reStruct <- reOut$reStruct
  random.names <- reOut$random.names
  expert.time <- reOut$expert.time

  # update dim.list
  ## TODO: develop time component of .cpp to distinguish between gating/expert.
  dim.list$n.t <- length(unique(expert.time))

  # remove intercept from gatingformula when random effects specified
  if (sum(reStruct[1, ] > 0) &
    attributes(terms(gatingformula))$intercept == 1) {
    gatingformula <- update(gatingformula, ~ . - 1)
    message("intercept removed from gatingformula
            when random effects specified")
  }



  ## Rank reduction settings
  # clustTMB currently allows for spatial and random rank reduction
  rrStruct <- c(0, 0)
  if (!is.null(rr$random)) {
    if (rr$random >= dim.list$n.j) {
      stop("random rank reduction must be smaller
           than the number of columns in the response")
    }
    rrStruct[1] <- 1
    dim.list$n.f.rand <- rr$random
    random.names <- c(random.names, "v_ifg")
  } else {
    dim.list$n.f.rand <- dim.list$n.j
  }
  dim.list$nl.rand <- dim.list$n.j * dim.list$n.f.rand -
    (dim.list$n.f.rand * (dim.list$n.f.rand - 1)) / 2
  if (!is.null(rr$spatial)) {
    if (rr$spatial >= dim.list$n.j) {
      stop("spatial rank reduction must be smaller than
           the number of columns in the response")
    }
    rrStruct[2] <- 1
    dim.list$n.f.sp <- rr$spatial
  } else {
    dim.list$n.f.sp <- dim.list$n.j
  }
  dim.list$nl.sp <- dim.list$n.j * dim.list$n.f.sp -
    (dim.list$n.f.sp * (dim.list$n.f.sp - 1)) / 2
  if (sum(rrStruct) > 0) {
    if (dim.list$n.j == 1) {
      stop("cannot implement rank reduction on univariate models")
    }
    if (covariance.structure == "EEE" | covariance.structure == "VVV") {
      stop("Need to specify diagonal covariance structure
           when implementing rank reduction")
    }
  }

  # if((rrStruct[1] == 1) & reStruct[2,3] == 0){
  #   stop ("You have specified an random error rank reduction in the expert model and therefore need to specify a normal overdispersion model in the expertformula")
  # }
  if (rrStruct[2] == 1 & reStruct[2, 1] == 0) {
    stop("You have specified a spatal rank reduction in the
         expert model and therefore need to specify a
         spatial model in the expertformula")
  }
  ## ! fix time component of .cpp to distinguish between gating/expert. Implement something similar to glmmTMB?
  Dat <- mkDat(
    response,
    time.vector = expert.time - 1,
    expert.data = expertdata,
    gating.data = gatingdata,
    expert.formula = expertformula,
    gating.formula = gatingformula,
    family = family,
    ll.method = ll.method,
    fixStruct = covariance.structure,
    rrStruct = rrStruct,
    reStruct = reStruct,
    dim.list = dim.list,
    offset = Offset,
    spatial.list = spatial.list,
    projection.dat = projection.dat
  )
  dim.list$n.v <- Dat$spde$n_s
  initialization.args$Data <- Dat
  initialization.args$dim.list <- dim.list
  initialization.args$family <- family
  init.parm <- do.call(genInit, initialization.args)
  arg.map <- mkMap(
    Dat$family,
    covariance.structure,
    Dat$rrStruct,
    Dat$reStruct,
    dim.list
  )
  # update starting values
  for (p in names(Start)) {
    if (!(p %in% names(init.parm$parms))) {
      stop(sprintf("unrecognized parameter name in Start: %s", p))
    }
    if (!(p %in% start.names())) {
      stop(sprintf("setting initial value unsupported for this
                   parameter as initial value controlled by user
                   specified cluster id: %s", p))
    }
    if (p %in% names(arg.map)) {
      warning(sprintf("parameter is not estimated in specified model,
                      setting starting value may affect inference: %s", p))
    }
    Ds <- dim(Start[[p]])
    Dp <- dim(init.parm$parms[[p]])
    if (is.null(Ds)) Ds <- length(Start[[p]])
    if (is.null(Dp)) Dp <- length(init.parm$parms[[p]])
    if (!all(Ds == Dp)) {
      stop(sprintf("parameter dimension mismatch,
                   see parm.lookup() for dimensions: %s", p))
    }
    init.parm$parms[[p]] <- Start[[p]]
  }
  for (m in names(Map)) {
    if (!(m %in% names(init.parm$parms))) {
      stop(sprintf("unrecognized parameter name in Map: %s", m))
    }
    if (!is.factor(Map[[m]])) {
      stop(sprintf("map values need to be specified as factors,
                   see ?TMB::MakeADFun() for details: %s", m))
    }
    if (m %in% names(arg.map)) {
      if (all(is.na(arg.map[[m]])) & !(all(is.na(Map[[m]])))) {
        stop(sprintf("this parameter is turned off for this model
                     specification and should not be estimated: %s", m))
      }
    }
    Dm <- dim(Map[[m]])
    Dp <- dim(init.parm$parms[[m]])
    if (is.null(Dm)) Dm <- length(Map[[m]])
    if (is.null(Dp)) Dp <- length(init.parm$parms[[m]])
    if (!all(Dm == Dp)) {
      stop("parameter dimension mismatch, see parm.lookup() for dimensions")
    }
    arg.map[[m]] <- Map[[m]]
  }
  if (control$check.input == TRUE) {
    clustTMB.mod <- list(
      Dat = Dat,
      inits = init.parm,
      map = arg.map,
      random = random.names
    )
  } else {
    clustTMB.mod <- fit.tmb(
      obj.args = list(
        data = Dat,
        parameters = init.parm$parms,
        map = arg.map,
        random = random.names,
        DLL = "clustTMB",
        silent = TRUE
      ),
      control = control
    )
  }
  return(clustTMB.mod)
}

#' Run Options
#'
#' @param check.input TRUE: Return initial values before running TMB
#' @param run.model FALSE: Return TMB object before optimizing model
#' @param do.sdreport TRUE: Run delta method to obtain standard errors
#'
#' @return list
#' @export
#' @examples
#' run.options()
run.options <- function(check.input = NULL, run.model = NULL,
                        do.sdreport = NULL) {
  if (missing(check.input)) check.input <- FALSE
  if (missing(run.model)) run.model <- TRUE
  if (missing(do.sdreport)) do.sdreport <- TRUE
  return(list(
    check.input = check.input,
    run.model = run.model,
    do.sdreport = do.sdreport
  ))
}
