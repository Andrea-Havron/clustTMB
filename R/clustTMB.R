#' Fit a finite mixture model using TMB
#'
#'
#' @param response as vector or matrix
#' @param expertformula formula defining expert model (response)
#' @param gatingformula formula defining gating model (clusters)
#' @param expertdata data frame containing expert model covariates
#' @param gatingdata data frame containing gating model covariates
#' @param family statistical distribution and link functin
#' @param G number of cluster
#' @param rr factor analysis list identifying random effect rank reduction
#' @param structure fixed effects covariance structure using mclust naming scheme
#' @param initialization how to generate initial values
#' @param spatial.list a list of data objects needed when fitting a spatial GMRF model
#' @param projection.list a list of data objects needed for projection
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
#' m1 <- clustTMB(response = faithful, structure = "VVV",
#' initialization.args = list(init.method = "hc", hc.options = c("VVV", "VARS")))
#' plot(faithful$eruptions, faithful$waiting, pch = 16, col = m1$report$classification)
clustTMB <- function(response = NULL,
                     expertformula = ~ 1,
                     gatingformula = ~ 1,
                     expertdata = NULL,
                     gatingdata = NULL,
                     family = gaussian(link = "identity"),
                     G = 2,
                     rr = list(spatial = NULL, temporal = NULL, dispersion = NULL),
                     structure = NULL,
                     ##! more advanced arguments can be passed to liist. How do I specify this, using ...?
                     initialization.args = list(init.method = "hc",
                                                data.trans = NA,
                                                ##! default need to be 'E' when response univariate
                                                hc.options = c("VVV", "VARS"),
                                                mix.method = NULL),
                     spatial.list = list(loc = NULL, mesh = NULL),
                     projection.list = list(grid.df = NULL, ##!Need more rules about grid.df spatial structure
                                            ##!need to ensure order of covariates preserved somehow
                                            ##?match names from expert.dat/gating.dat to grid.df?
                                            expert.pred.names = NULL,
                                            gating.pred.names = NULL)){

  response <- as.matrix(response)
  dim.list <- list(n.i = nrow(response), n.j = ncol(response), n.t = NULL,
                   n.g = G, n.f = NULL, n.v = NULL)
  #rr is only spatial at the moment
  if(!is.null(rr$spatial)){
    dim.list$n.f = rr$spatial
  } else {
    dim.list$n.f = dim.list$n.j
  }
  if(is.null(expertdata)) expertdata <- data.frame(x=rep(1,dim.list$n.i))
  if(is.null(gatingdata)) gatingdata <- data.frame(x=rep(1,dim.list$n.i))
  if(!is.data.frame(expertdata)){
    stop ("expert data needs to be a data.frame")
  }
  if(!is.data.frame(gatingdata)){
    stop ("gating data needs to be a data.frame ")
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
  if(length(expert.re.names >0)){
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
  covstruct <- matrix(0, 2, 3)
  random.names <- c()

  for(i in seq_along(gating.re.names)){
    if(gating.re.names[i] == "gmrf"){
      covstruct[1,1] <- 3
      random.names <- c(random.names, "Gamma_vg")
    }
    if(gating.re.names[i] == "gmrfSpeedup"){
      covstruct[1,1] <- 4
      random.names <- c(random.names, "Gamma_vg")
    }
    if(gating.re.names[i] == "ar1"){
      covstruct[1,2] <- 2
      random.names <- c(random.names, "upsilon_tg")
    }
  }
  for(i in seq_along(expert.re.names)){
    if(expert.re.names[i] == "gmrf"){
      covstruct[2,1] <- 3
      random.names <- c(random.names, "Omega_vfg")
    }
    if(expert.re.names[i] == "gmrfSpeedup"){
      covstruct[2,1] <- 4
      random.names <- c(random.names, "Omega_vfg")
    }
    if(expert.re.names[i] == "ar1"){
      covstruct[2,2] <- 2
      random.names <- c(random.names, "epsilon_tjg")
    }
  }

  expert.fix.dat <- model.matrix(nobars(expertformula))
  gating.fix.dat <- model.matrix(nobars(gatingformula))
  if( (length(dimnames(expert.fix.dat)[[2]]) == 1) &
      dimnames(expert.fix.dat)[[2]] == "(Intercept)" ){
    expert.fix.dat <- matrix(1, dim.list$n.i, 1)
  }
  if( (length(dimnames(gating.fix.dat)[[2]]) == 1) &
      dimnames(gating.fix.dat)[[2]] == "(Intercept)" ){
    gating.fix.dat <- matrix(1, dim.list$n.i, 1)
  }

  ##! fix time component of .cpp to distinguish between gating/expert. Implement something similar to glmmTMB?
  Dat <- mkDat(response, time.vector = expert.time-1,
               expert.dat = expert.fix.dat,
               gating.dat = gating.fix.dat,
               family = family,
               fixStruct = structure, covstruct = covstruct,
               dim.list =  dim.list,
               spatial.list = spatial.list,
               projection.list = projection.list)
  dim.list$n.v <- Dat$spde$n_s
  initialization.args$Data <- Dat
  initialization.args$dim.list <- dim.list
  clustTMB.mod <- fit.tmb(
    ##! Passing makeADFun args to to fit.tmb. Is there are more elegant way to do this? Lots of other options can be passed, need to include this somehow
         list(data = Dat,
         parameters = do.call(genInit, initialization.args)$parms,
         map = mkMap(Dat$covstruct, Dat$family, Dat$fixStruct, dim.list),
         random = random.names,
         DLL="clustTMB",
         silent = TRUE)
  )

  return(clustTMB.mod)
}


