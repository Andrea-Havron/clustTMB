#' Helper function to run TMB code and return output
#' @importFrom TMB openmp MakeADFun sdreport
#'
#' @param obj.args List of arguments for TMB MakeADFun() function
#' @param opt.args List of arguments for nlminb() function
#' @param control List controling model runs and standard error reporting
#'
#' @return Fitted objective function, nlminb output, reported values from model, sdreport if true
#'
#' @noRd
fit.tmb <- function(obj.args, opt.args = list(control = list(iter = 800, eval = 800),
                                              hessian = NULL, scale = 1,
                                              lower = -Inf, upper = Inf ),
                    control = list(run.model = TRUE, do.sdreport = TRUE)){
  obj <- do.call(MakeADFun, obj.args)
  if(control$run.model){
    opt <- with(obj, do.call(nlminb,  c(list(par, fn, gr), opt.args) ))
    report <- obj$report(obj$env$last.par.best)
    if(control$do.sdreport){
      sdr <- sdreport(obj)
    } else {
      sdr <- NULL
    }
  } else {
    opt <- NULL
    report <- NULL
    sdr <- NULL
  }
  fit.results <- list(obj = obj, opt = opt, report = report, sdr = sdr)
  return(fit.results)
}

#' Returns factors defining how to optionally collect and fix parameters to be passed to map list in MakeADFun. See ?MakeADFun for more details.
#'
#' @param d Dimension of factors for mapping
#' @param f Vector of values to be converted to a factor
#'
#' @return Integer or vector of factors
#' @keywords internal
#' @noRd
mkFac <- function(d, f = NA){
  ans <- factor(f)
  if(length(d)>1){
    dim(ans) <- d
  }
  return(ans)
}


#' Sets up map list for TMB model. The map list defines how to collect and fix parameters. See ?MakeADFun for more details.
#'
#' @param Family Distribution family
#' @param covstruct Covariance structure
#' @param rrStruct Rank reduction structure
#' @param reStruct Random effects structure
#' @param dim.list List of parameter dimensions
#' @param map.ops In future, allow users to alter map?
#'
#' @return List defining how to optionally collect and fix parameters
#' @keywords internal
#' @noRd
mkMap <- function(Family, covstruct, rrStruct, reStruct, dim.list, map.ops = NULL){
  #list2env(dim.list, environment(mkMap))
  n.i <- dim.list$n.i
  n.j <- dim.list$n.j
  n.t <- dim.list$n.t
  n.g <- dim.list$n.g
  n.f.sp <- dim.list$n.f.sp
  n.f.rand <- dim.list$n.f.rand
  n.v <- dim.list$n.v
  nl.rand <- n.j*n.f.rand - (n.f.rand*(n.f.rand-1))/2
  nl.sp <- n.j*n.f.sp - (n.f.sp*(n.f.sp-1))/2
  nl.fix <- ifelse(n.j>1, (n.j^2-n.j)/2, 1)

  #Map out parameters based on model structure
  Map <- list()
  # if(sum(reStruct[1,]>0)){ ##! Rather than map out intercept, change formula to remove intercept and include warning
  #   Map$betag = mkFac(d = c(1,n.g-1), f = rep(NA, n.g-1)) #map out betag if random effects terms
  # }
  if(Family != 700) Map$thetaf <- mkFac(d = c(n.j,n.g), f = rep(NA, n.j*n.g))
  if( covstruct != 30  ){ # if not a general covariance structure
    Map$logit_corr_fix <- mkFac(d = c(nl.fix, n.g), f = rep(NA, nl.fix*n.g))
  }
  if(rrStruct[1] == 0){
    Map$ld_rand <- mkFac(d = c(nl.rand, n.g), f = rep(NA, nl.rand*n.g))
  }
  if(rrStruct[2] == 0){
    Map$ld_sp <- mkFac(d = c(nl.sp, n.g), f = rep(NA, nl.sp*n.g))
  }
  if(covstruct == 'EII'){ #Diagonal, equal variance within and between cluster
    Map$theta <- mkFac(d = c(n.j, n.g), f = rep(1, n.j*n.g))
  }
  if(covstruct == 'VII'){ #Diagonal, equal variance within cluster, variable variance between clusters
    Map$theta <- mkFac(d = c(n.j, n.g),
                       f = as.vector(
                         matrix(rep(1:n.g, each = n.j), n.j, n.g)
                       ))
  }
  if(covstruct == 'EEI'){ #Diagonal, variable variance within cluster, equal variance between clusters
    Map$theta <- mkFac(d = c(n.j, n.g),
                       f = as.vector(
                         matrix(rep(1:n.j, each = n.g), n.j, n.g, byrow = TRUE)
                       ))
  }
  if(covstruct == 'EEE'){ #General, variable variance and covariance within cluster, equal variance and covariance between clusters
    Map$logit_corr_fix <- mkFac(d = c(nl.fix, n.g),
                                f = as.vector(
                                  matrix(rep(1:nl.fix, each = n.g),
                                         nl.fix, n.g, byrow = TRUE)
                                  ))
    Map$theta <- mkFac(d = c(n.j, n.g),
                       f = as.vector(
                         matrix(rep(1:n.j, each = n.g),
                                n.j, n.g, byrow = TRUE)
                       ))
  }
  if(reStruct[1,1] == 0){
    Map$Hg_input <- mkFac(d = c(2,n.g-1), f = rep(NA, 2*(n.g-1)))
    Map$ln_kappag <- mkFac(d = n.g-1, f = rep(NA, n.g-1))
    Map$Gamma_vg <- mkFac(d = c(n.v,n.g-1), f = rep(NA, n.v*(n.g-1)))
  }
  if(reStruct[1,1] > 0){
    Map$Hg_input = mkFac(d = c(2,n.g-1), f = rep(NA, 2*(n.g-1)))
    Map$ln_kappag <- mkFac(d = n.g-1, f = rep(1, n.g-1))
  }
  if(reStruct[2,1] == 0){
    Map$Hd_input <- mkFac(d = c(2,n.j, n.g), f = rep(NA,2*n.j*n.g))
    Map$ln_kappad <- mkFac(d = c(n.j, n.g), f = rep(NA, n.j*n.g))
    Map$ln_taud <- mkFac(d = c(n.j, n.g), f = rep(NA, n.j*n.g))
    Map$Omega_vfg <- mkFac(d = c(n.v, n.j, n.g), f = rep(NA, n.v*n.j*n.g))
  }
  if(reStruct[2,1] > 0){
   # use when anisotropy turned on Map$Hd_input <- mkFac(d = c(2,n.j, n.g), f = rep(rep(c(1,2),n.j),n.g))
    Map$Hd_input <- mkFac(d = c(2,n.j, n.g), f = rep(NA,2*n.j*n.g))
    Map$ln_kappad <- mkFac(d = c(n.j, n.g), f = rep(1, n.j*n.g))
    Map$ln_taud <- mkFac(d = c(n.j, n.g), f = rep(1, n.j*n.g))
  }
  if(rrStruct[2] == 1){
    Map$ln_taud <- mkFac(d = c(n.j,n.g), f = rep(NA, n.j*n.g))
  }
  if(reStruct[1,2] == 0){
    Map$logit_rhog <- mkFac(d = c(n.g-1), f = rep(NA, n.g-1))
    Map$ln_sigmaup <- mkFac(d = c(n.g-1), f = rep(NA, n.g-1))
    Map$upsilon_tg <- mkFac(d = c(n.t,n.g-1), f = rep(NA, n.t*(n.g-1)))
  }
  if(reStruct[2,2] == 0){
    Map$logit_rhod <- mkFac(d = c(n.j, n.g), f = rep(NA, n.j*n.g))
    Map$ln_sigmaep <-  mkFac(d = c(n.j, n.g), f = rep(NA, n.j*n.g))
    Map$epsilon_tjg <- mkFac(d = c(n.t,n.j,n.g), f = rep(NA, n.t*n.j*n.g))
  }
  if(reStruct[1,3] == 0){
    Map$ln_sigmau <- mkFac(d = c(n.g-1), f = rep(NA, n.g-1))
    Map$u_ig <- mkFac(d = c(n.i,n.g-1), f = rep(NA, n.i*(n.g-1)))
  }
  if(reStruct[2,3] == 0){
    Map$ln_sigmav <- mkFac(d = c(n.j, n.g), f = rep(NA, n.j*n.g))
    Map$v_ifg <- mkFac(d = c(n.i,n.j,n.g), f = rep(NA, n.i*n.j*n.g))
  }
  if(rrStruct[1] == 1){
    Map$ln_sigmav <- mkFac(d = c(n.j, n.g), f = rep(NA, n.j*n.g))
  }
  return(Map)
}



#' @importFrom stats make.link
#' @keywords internal
#' @noRd
Tweedie <- function(f = 'Tweedie', link = "log"){
  f <- c(f,list(link=link),make.link(link))
  class(f) <- "family"
  return(f)
}

#' Sets up data list for TMB model.
#'
#' @param response A numeric matrix of observations.
#'
#' @param time.vector A numeric indicator vector specifying time classification.
#' @param expert.dat A dataframe specifying covariates in the expert model. If no expret model, defaults to a unit matrix.
#' @param gating.dat A dataframe specifying covariates in the gating model. If no gating model, defaults to a unit matrix.
#' @param family Distribution family
#' @param ll.method Integer specifying whether to use marginal (0) or conditional (1) lieklihood during optimization
#' @param fixStruct Integer specigying whether covariance structure is Univariate, Diagonal, or General
#' @param rrStruct Integer vector specifying whether or not rank reduction is turned on for each random effect
#' @param reStruct Integer matrix specifying random effects structure
#' @param dim.list List of parameter dimensions
#' @param offset Expertformula offset term
#' @param spatial.list List of spatial objects when fitting a spatial model
#' @param projection.list List of spatial objects used when returning spatial projection results
#'
#' @importFrom INLA inla.mesh.create inla.spde.make.A
#' @importFrom sp coordinates
#'
#' @return Data list for input into TMB::MakeADFun
#' @keywords internal
#' @noRd

mkDat <- function(response, time.vector, expert.dat, gating.dat,
                  family = gaussian(link = "identity"), ll.method,
                  fixStruct, rrStruct, reStruct, dim.list,
                  offset = NULL,
                  spatial.list = list(loc = NULL, mesh = NULL),
                  projection.list = list(grid.df = NULL, ##!Need more rules about grid.df spatial structure
                                         ##!need to ensure order of covariates preserved somehow
                                         ##?match names from expert.dat/gating.dat to grid.df?
                                         expert.pred.names = NULL,
                                         gating.pred.names = NULL)){
  n.i <- dim.list$n.i
  n.j <- dim.list$n.j
  n.t <- dim.list$n.t
  n.g <- dim.list$n.g
  n.f <- dim.list$n.f
  n.v <- dim.list$n.v
  #list2env(dim.list, environment(mkDat)) ##! environment locked. try new.env?
  loc <- spatial.list$loc
  mesh <- spatial.list$mesh
  #list2env(spatial.list, environment(mkDat))
  grid.df <- projection.list$grid.df
  expert.pred.names <- projection.list$expert.pred.names
  gating.pred.names <- projection.list$gating.pred.names
  #list2env(projection.list)
  if((is.null(mesh)&is.null(loc)) & !is.null(grid.df)){
    warning("loc and mesh are null. Need to provide locations or mesh in spatial.list to initiate spatial model for spatial predictions")
  }
  if((!is.null(mesh)|!is.null(loc)) & is.null(grid.df)){
    warning("grid.df is null. Need to provide locations in projection.list$grid.df for spatial predictions")
  }
  if(is.null(loc)&is.null(mesh)){
    #create dummy mesh for model - not used in inference or projection
    loc <- data.frame(x=runif(n.i,0,1), y=runif(n.i,0,1))
    mesh <- inla.mesh.create(as.matrix(loc))
    coordinates(loc) <- ~x*y
  }
  if(!is.null(loc)& is.null(mesh)){
    #default mesh - in future add options to include arguments for inla.mesh.2d
    #for now, user can supply mesh if a more complex mesh is needed
    mesh <- inla.mesh.create(loc@coords)
    warning("Building simple spatial mesh. If using the SPDE-FEM GMRF method,
            the simple mesh may result in spatial bias. Consider bulding a
            more appropriate mesh using INLA::meshbuilder()")
  }
  if(is.null(loc)&!is.null(mesh)){
    if(is.null(mesh$idx$loc)){
      #if user-supplied mesh built without observation locations, the user must also provide observation lovations
      stop ("Need to provide locations for observations in spatial.list$loc")
    } else {
      #if user-supplied mesh built using observation locations, these can be obtained from the mesh if not provided
      loc <- mesh$loc[mesh$idx$loc,1:2]
    }
  }
  n.v <- mesh$n
  A <- inla.spde.make.A(mesh, loc@coords)

  if(is.null(grid.df)){
    #create dummy projection grid - ##! be sure to turn off projection when reporting out
    grid.loc <- as.matrix(expand.grid(x = seq(0,1,.1), y = seq(0,1,.1)))
    Xd_proj = matrix(1, nrow(grid.loc), 1)
    Xg_proj = matrix(1, nrow(grid.loc), 1)
  } else {
    grid.loc <- as.matrix(grid.df@coords)
    if(class(grid.df) == "SpatialPoints" ){ ##!is this the best way to do this?
      Xd_proj = matrix(1, nrow(grid.loc), 1)
      Xg_proj = matrix(1, nrow(grid.loc), 1)
    } else {
      grid.data <- grid.df@data
      Xd_proj <- as.matrix(grid.data[expert.pred.names])
      Xg_proj <- as.matrix(grid.data[gating.pred.names])
    }
  }
  A.proj <- inla.spde.make.A(mesh, grid.loc)


  if(is.null(offset)) offset = rep(1, dim.list$n.i)

  Dat <- list(
    Y = as.array(response),
    t = time.vector,
    Xd = expert.dat,
    Xg = gating.dat,
    Xpz = matrix(1, n.i, 1),
    Offset = offset,
    A = A,
    A_proj = A.proj,
    Xd_proj = matrix(1, nrow(grid.loc), 1),
    Xg_proj = matrix(1, nrow(grid.loc), 1)
  )

  Dat$spde <- spdeStruct(mesh)
  if(family[[1]] == "gaussian") Dat$family <- 0
  if(family[[1]] == "Tweedie") Dat$family <- 700
  if(family[[2]] == "log") Dat$link <- 0
  if(family[[2]] == "identity") Dat$link <- 5
  Dat$loglike <- ll.method
  if(fixStruct == 'E' | fixStruct == 'V'){
    Dat$fixStruct <- 10
  }
  if(fixStruct == 'EII' | fixStruct == 'VII' | fixStruct == 'EEI' | fixStruct == 'VVI'){
    Dat$fixStruct <- 20
  }
  if(fixStruct == "VVV" | fixStruct == 'EEE'){
    Dat$fixStruct <- 30
  }
  Dat$rrStruct <- rrStruct
  Dat$reStruct <- reStruct


  return(Dat)
}


fixStruct.names <- function(...){
  return(c('E', 'V', 'EII', 'VII', 'EEI', 'VVI', 'VVV', 'EEE', 'RR'))
}
