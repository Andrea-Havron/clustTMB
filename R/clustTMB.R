#' Fit a finite mixture model using TMB
#'
#'
#' @param response A numeric vector, matrix, or data frame of observations. When data are multivariate, rows correspond to observations and columns correspond to the multivariate response.
#' @param expertformula Formula defining expert model. This formula corresponds to the covariates included in the response densities. Defaults to intercept only (~1) when no covariates are used.
#' @param gatingformula formula defining gating model.  This formula corresponds to the covariates included in the mixing proportions (logistic regression). Defaults to intercept only (~1) when no covariates are used. When a random effects term is included in the gating network, this formula will be updated so that the intercept term is removed.
#' @param expertdata Data frame containing expert model covariates.
#' @param gatingdata Data frame containing gating model covariates.
#' @param family Statistical distribution and link function of observations.
#' @param Offset Constant in expertformula only used to offset density expectation.
#' @param G Integer specifying the number of clusters.
#' @param rr List specifying dimension of rank reduction in spatial, temporal, and/or random effects. Dimension must be smaller than the total dimension of the reponse. Rank reduction is applied only to the expertformula random effects. The rank reduction reduces the dimensionality of a correlated multivariate ressponse to a smaller dimension independent reponse. When used, the covariance structure of the response is swicthed to 'Diagonal.' Defaults to NULL, no rank reduction. If rank reduction is used in conjunction with a random effect, that random effect must also be specified in the expert formula. Currently, rank reduction on temporal random effects is disabled.
#' @param covariance.structure A character string specifying the covariance structure of the response using mclust naming scheme. See description of modelNames under ?Mclust for details.
#' @param initialization.args A list consisting of initalization settings used to generate initial values.
##' \item{init.method}{methods include: qclass (univariate), hc, kmeans, mixed}
##' \item{data.trans}{data transformation options include: log+1, log, yeojin boxcox }
##' \item{hc.options}{list specifying (hcModelName and hcUse) from mclust, see ?mclust.options.}
##' \item{mix.method}{When fitting Tweedie data, this method is used to determine a better intial clustering. Options include kproto, Gower kmeans, and Gower hclust.}
#' @param spatial.list List of data objects needed when fitting a spatial GMRF model
#' @param projection.list List of data objects needed for projection
#' @param control List controlling whether models are run and whether standard errors are calculated.
#' @importFrom stats model.matrix model.frame as.formula setNames terms update vcov
#' @importFrom lme4 subbars nobars findbars mkReTrms
#' @importFrom glmmTMB splitForm
#' @importFrom utils data
#' @importFrom Matrix t
#' @useDynLib clustTMB
#'
#' @return
#' @export
#'
#' @examples
#' data("faithful")
#' m1 <- clustTMB(response = faithful, covariance.structure = "VVV",
#' initialization.args = list(init.method = "hc", hc.options = c("VVV", "VARS")))
#' plot(faithful$eruptions, faithful$waiting, pch = 16, col = m1$report$classification)
clustTMB <- function(response = NULL,
                     expertformula = ~ 1,
                     gatingformula = ~ 1,
                     expertdata = NULL,
                     gatingdata = NULL,
                     family = gaussian(link = "identity"),
                     Offset = NULL,
                     G = 2,
                     rr = list(spatial = NULL, temporal = NULL, random = NULL),
                     covariance.structure = NULL,
                     ##! more advanced arguments can be passed to list. How do I specify this, using ...?
                     initialization.args = list(init.method = NULL,
                                                data.trans = NA,
                                                ##! default needs to be 'E' when response univariate
                                                hc.options = c(),
                                                mix.method = NULL),
                     spatial.list = list(loc = NULL, mesh = NULL),
                     projection.list = list(grid.df = NULL, ##!Need more rules about grid.df spatial structure
                                            ##!need to ensure order of covariates preserved somehow
                                            ##?match names from expert.dat/gating.dat to grid.df?
                                            expert.pred.names = NULL,
                                            gating.pred.names = NULL),
                     control = list(run.model = TRUE, do.sdreport = TRUE)){

  response <- as.matrix(response)
  dim.list <- list(n.i = nrow(response), n.j = ncol(response), n.t = NULL,
                   n.g = G, n.f.sp = NULL, n.f.rand = NULL, n.v = NULL)
  if(is.null(covariance.structure)){
    if(dim.list$n.j == 1){
      covariance.structure = 'E'
      warning("Setting covariance structure to univariate E")
    } else {
      covariance.structure = 'EII'
      warning("Setting covariance structure to univariate EII")
    }
  }
  if(gatingformula == ~1){
    ll.method <- 0
  } else {
    ll.method <- 2
  }
  if(is.null(expertdata)) expertdata <- data.frame(x=rep(1,dim.list$n.i))
  if(is.null(gatingdata)) gatingdata <- data.frame(x=rep(1,dim.list$n.i))
  if(!is.data.frame(expertdata)){
    stop ("expert data needs to be a data.frame")
  }
  if(!is.data.frame(gatingdata)){
    stop ("gating data needs to be a data.frame")
  }
  if(!is.null(expertdata)){
    if(nrow(expertdata) != dim.list$n.i){
      stop ("expert data and response need to have the same number of observations")
    }
  }
  if(!is.null(gatingdata)){
    if(nrow(gatingdata) != dim.list$n.i){
      stop ("gating data and response need to have the same number of observations")
    }
  }
  if((covariance.structure == 'E' | covariance.structure == 'V' ) & dim.list$n.j > 1){
    stop("Need to specify multivariate covariance structure")
  }
  if(!(covariance.structure == 'E' | covariance.structure == 'V' ) & dim.list$n.j == 1){
    stop("Need to specify univariate covariance structure")
  }
  if(!(family[[1]] %in% names(.valid_family))){
    stop ("family not supported")
  }
  if(!(family[[2]] %in% names(.valid_link))){
    stop ("link not supported in specified family")
  }
  if(!is.null(spatial.list$loc)){
    if( (class(spatial.list$loc) != 'SpatialPoints') &
        (class(spatial.list$loc) != 'SpatialPointsDataFrame') ){
    stop ("Locations need to be SpatialPoints of SpatialPointsDataFrame class")
    }
  }


  #coding this way might increase risk of mismatched data and coordinates due to reordering
  #when data are spatial, it might be more interesting to include everything in as a spatial data frame
  #or is it better to keep structure simple using glmmTMB's numFactor method?
  if(!is.null(spatial.list$loc)){
    expertdata$loc <- factor(row.names(spatial.list$loc@coords))
    gatingdata$loc <- factor(row.names(spatial.list$loc@coords))
  }

  expert.split <- splitForm(expertformula)
  expert.re.names <- expert.split$reTrmClasses
  gating.split <- splitForm(gatingformula)
  gating.re.names <- gating.split$reTrmClasses

  if(('gmrf' %in% expert.re.names | 'gmrfSpeedup' %in% expert.re.names) &
     (is.null(spatial.list$loc) & is.null(spatial.list$mesh))){
    stop ("You have specified a spatal model and need to provide location or
          mesh data in the spatial.list argument")
  }
  if(('gmrf' %in% gating.re.names | 'gmrfSpeedup' %in% gating.re.names) &
     (is.null(spatial.list$loc) & is.null(spatial.list$mesh))){
    stop ("You have specified a spatal model and need to provide location or
          mesh data in the spatial.list argument")
  }

  specials <- c("ar1", "gmrf", "gmrfSpeedup")
  if(length(expert.re.names)>0){
    for(i in 1:length(expert.re.names)){
      if(!(expert.re.names[i] %in% specials)){
        stop ("Currently clustTMB only works with spatio-temporal random effects")
      }
    }
  }
  if(length(gating.re.names)>0){
    for(i in 1:length(gating.re.names)){
      if(!(gating.re.names[i] %in% specials)){
        stop ("Currently clustTMB only works with spatio-temporal random effects")
      }
    }
  }

  if('ar1' %in% expert.re.names){
    idx <- which(expert.re.names == 'ar1')
    ar1.form <- as.formula(paste("~", deparse(expert.split$reTrmFormulas[[idx]])))
    ##! make sure this enters TMB as numeric - although consider using glmmTMB factor approach
    expert.time <- model.frame(subbars(ar1.form), expertdata)
  } else {
    expert.time <- rep(1,dim.list$n.i)
  }
  ##! fix time component of .cpp to distinguish between gating/expert. Implement something similar to glmmTMB?
  dim.list$n.t <- length(unique(expert.time))

  if('ar1' %in% gating.re.names){
    idx <- which(gating.re.names == 'ar1')
    ar1.form <- as.formula(paste("~", deparse(gating.split$reTrmFormulas[[idx]])))
    gating.time <- model.frame(subbars(ar1.form), gatingdata)
  } else {
    gating.time <- rep(1,dim.list$n.i)
  }

  if(('gmrf' %in% expert.re.names)|('gmrfSpeedup' %in% expert.re.names)){
    idx <- which((expert.re.names == 'gmrf') | (expert.re.names == 'gmrfSpeedup'))
    gmrf.form <- as.formula(paste("~", deparse(expert.split$reTrmFormulas[[idx]])))
    expert.gmrf <- model.frame(subbars(gmrf.form), expertdata)
  } else {
    expert.gmrf <- NA
  }

  if(('gmrf' %in% gating.re.names)|('gmrfSpeedup' %in% gating.re.names)){
    idx <- which((gating.re.names == 'gmrf') | (gating.re.names == 'gmrfSpeedup'))
    gmrf.form <- as.formula(paste("~", deparse(gating.split$reTrmFormulas[[idx]])))
    gating.gmrf <- model.frame(subbars(gmrf.form), gatingdata)
  } else {
    gating.gmrf <- NA
  }

  ##! expert/gating .ar1 and .gmrf provide information about interactions between space/time
  ##!! Spatio-temporal interactions not implemented yet!!
  reStruct <- matrix(0, 2, 3)
  random.names <- c()

  for(i in seq_along(gating.re.names)){
    if(gating.re.names[i] == "gmrf"){
      reStruct[1,1] <- 3
      random.names <- c(random.names, "Gamma_vg")
    }
    if(gating.re.names[i] == "gmrfSpeedup"){
      reStruct[1,1] <- 4
      random.names <- c(random.names, "Gamma_vg")
    }
    if(gating.re.names[i] == "ar1"){
      reStruct[1,2] <- 2
      random.names <- c(random.names, "upsilon_tg")
    }
  }
  for(i in seq_along(expert.re.names)){
    if(expert.re.names[i] == "gmrf"){
      reStruct[2,1] <- 3
      random.names <- c(random.names, "Omega_vfg")
    }
    if(expert.re.names[i] == "gmrfSpeedup"){
      reStruct[2,1] <- 4
      random.names <- c(random.names, "Omega_vfg")
    }
    if(expert.re.names[i] == "ar1"){
      reStruct[2,2] <- 2
      random.names <- c(random.names, "epsilon_tjg")
    }
  }
  if(sum(reStruct[1,]>0) & attributes(terms(gatingformula))$intercept == 1){
    gatingformula <- update(gatingformula, ~ . -1)
    warning('intercept removed from gatingformula when random effects specified')
  }

  expert.fix.dat <- model.matrix(nobars(expertformula))
  gating.fix.dat <- model.matrix(nobars(gatingformula))
  if( (length(dimnames(expert.fix.dat)[[2]]) == 1) &
      dimnames(expert.fix.dat)[[2]][1] == "(Intercept)" ){
    expert.fix.dat <- matrix(1, dim.list$n.i, 1, dimnames = list(NULL, '(Intercept)'))
  }
  if( (length(dimnames(gating.fix.dat)[[2]]) == 1) &
      dimnames(gating.fix.dat)[[2]][1] == "(Intercept)" ){
    gating.fix.dat <- matrix(1, dim.list$n.i, 1, dimnames = list(NULL, '(Intercept)'))
  }

  ## Rank reduction settings
  #clustTMB currently allows for spatial and random rank reduction
  rrStruct = c(0,0)
  if(!is.null(rr$random)){
    if(rr$random == n.j){
      stop("random rank reduction must be smaller than the number of columns in the response")
    }
    rrStruct[1] <- 1
    dim.list$n.f.rand = rr$random
  } else {
    dim.list$n.f.rand = dim.list$n.j
  }
  if(!is.null(rr$spatial)){
    if(rr$spatial == n.j){
      stop("spatial rank reduction must be smaller than the number of columns in the response")
    }
    rrStruct[2] <- 1
    dim.list$n.f.sp = rr$spatial
  } else {
    dim.list$n.f.sp = dim.list$n.j
  }
  if(sum(rrStruct)>0){
    if(n.j == 1){
      stop('cannot implement rank reduction on univariate models')
    }
    if(covariance.structure != 'RR'){
      covariance.structure = 'Diag'
      warning('Covariance structure changed to Diag. Multivariate models will be reduced to conditionally independent.')
    }
  }
  if(sum(rrStruct)==0 & covariance.structure == 'RR'){
    stop('Need to specify dimensions of rank reduction using rr')
  }

  if((rrStruct[1] == 1) & reStruct[2,3] == 0){
    stop ("You have specified an random error rank reduction in the expert model and therefore need to specify a normal overdispersion model in the expertformula")
  }
  if(rrStruct[2] == 1 & reStruct[2,1] == 0){
    stop ("You have specified a spatal rank reduction in th expert model and therefore need to specify a spatial model in the expertformula")
  }

  ##! fix time component of .cpp to distinguish between gating/expert. Implement something similar to glmmTMB?
  Dat <- mkDat(response, time.vector = expert.time-1,
               expert.dat = expert.fix.dat,
               gating.dat = gating.fix.dat,
               family = family,
               ll.method = ll.method,
               fixStruct = covariance.structure,
               rrStruct = rrStruct,
               reStruct = reStruct,
               dim.list =  dim.list,
               offset = Offset,
               spatial.list = spatial.list,
               projection.list = projection.list)
  dim.list$n.v <- Dat$spde$n_s
  initialization.args$Data <- Dat
  initialization.args$dim.list <- dim.list
  initialization.args$family <- family
    clustTMB.mod <- fit.tmb(
      obj.args = list(data = Dat,
           parameters = do.call(genInit, initialization.args)$parms,
           map = mkMap(Dat$family, Dat$fixStruct, Dat$rrStruct, Dat$reStruct, dim.list),
           random = random.names,
           DLL="clustTMB",
           silent = TRUE),
      control = control
    )

  return(clustTMB.mod)
}


