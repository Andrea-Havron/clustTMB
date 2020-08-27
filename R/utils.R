#' @importFrom TMB openmp MakeADFun sdreport
#' @param obj.args list of arguments for TMB MakeADFun() function
#' @return fitted objective function, nlminb output, reported values from model
fit.tmb <- function(obj.args){
  obj <- do.call(MakeADFun, obj.args)
  opt <- with(obj, do.call(nlminb,  c(list(par, fn, gr, control = list(iter = 800, eval = 800))) ))
  report <- obj$report(obj$env$last.par.best)
  fit.results <- list(obj = obj, opt = opt, report = report)
  ##! add sdreport here or in a separate function?
  return(fit.results)
}

mkFac <- function(dims, f = NA){
  ans <- rep(factor(f), prod(dims))
  if(length(dims)>1){
    dim(ans) <- dims
  }
  return(ans)
}

mkMap <- function(covstruct, Family, fixStruct, dim.list, map.ops = NULL){
  #list2env(dim.list, environment(mkMap))
  n.i <- dim.list$n.i
  n.j <- dim.list$n.j
  n.t <- dim.list$n.t
  n.g <- dim.list$n.g
  n.f <- dim.list$n.f
  n.v <- dim.list$n.v
  if(fixStruct == 36){
    nl <- n.j*n.f - (n.f*(n.f-1))/2
  } else {
    nl <- (n.j^2-n.j)/2
  }
  Map <- list()
  if(sum(covstruct[1,]>0))  Map$betag = mkFac(c(1,n.g-1)) #map out betag if random effects terms
  if(Family != "Tweedie") Map$thetaf <- mkFac(c(n.j,n.g))
  if( fixStruct == 20 ) Map$ld <- mkFac(c(nl, n.g))
  if(covstruct[1,1] == 0){
    Map$Hg_input <- mkFac(c(2,n.g-1))
    Map$ln_kappag <- mkFac(n.g-1)
    Map$ln_taug <- mkFac(n.g-1)
    Map$Gamma_vg <- mkFac(c(n.v,n.g-1))
  }
  if(covstruct[1,1] > 0){
    Map$Hg_input = mkFac(c(2,n.g-1))
    Map$ln_kappag <- mkFac(n.g-1, f=1)
    Map$ln_taug <- mkFac(n.g-1, f=1)
  }
  if(covstruct[2,1] == 0){
    Map$Hd_input = mkFac(c(2,n.j,n.g))
    Map$ln_kappad <- mkFac(c(n.j,n.g))
    Map$ln_taud <- mkFac(c(n.j,n.g))
    Map$Omega_vfg <- mkFac(c(n.v, n.j, n.g))
  }
  if(covstruct[2,1] > 0){
    Map$Hd_input =  mkFac(c(2,n.j,n.g))
    Map$ln_kappad <- mkFac(c(n.j,n.g), f=1)
    Map$ln_taud <- mkFac(c(n.j,n.g), f=1)
  }
  if(fixStruct == 36){
    Map$ln_taud <- mkFac(c(n.j,n.g))
  }
  if(covstruct[1,2] == 0){
    Map$logit_rhog = mkFac(c(n.g-1))
    Map$ln_sigmaup = mkFac(c(n.g-1))
    Map$upsilon_tg = mkFac(c(n.t,n.g-1))
  }
  if(covstruct[2,2] == 0){
    Map$logit_rhod = mkFac(c(n.j,n.g))
    Map$ln_sigmaep = mkFac(c(n.j,n.g))
    Map$epsilon_tjg = mkFac(c(n.t,n.j,n.g))
  }
  if(covstruct[1,3] == 0){
    Map$ln_sigmau = mkFac(c(n.g-1))
    Map$u_ig = mkFac(c(n.i,n.g-1))
  }
  if(covstruct[2,3] == 0){
    Map$ln_sigmav = mkFac(c(n.j,n.g))
    Map$v_ifg = mkFac(c(n.i,n.j,n.g))
  }
  return(Map)
}

#' @importFrom stats make.link
Tweedie <- function(f = 'Tweedie', link = "log"){
  f <- c(f,list(link=link),make.link(link))
  class(f) <- "family"
  return(f)
}

#' @importFrom INLA inla.mesh.create inla.spde.make.A
mkDat <- function(response, time.vector, expert.dat, gating.dat,
                  family = gaussian(link = "identity"),
                  fixStruct, covstruct, dim.list,
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
  if(is.null(loc)&is.null(mesh)){
    if(is.null(grid.df)){
      #create dummy mesh and spde for model - not used in inference or projection
      loc <- matrix(runif(n.i*2,0,1), ncol = 2)
      mesh <- inla.mesh.create(loc)
      n.v <- mesh$n
      A <- inla.spde.make.A(mesh, loc)
      #create dummy projection grid - ##! be sure to turn off projection in reporting out
      grid.loc <- as.matrix(expand.grid(x = seq(0,1,.1), y = seq(0,1,.1)))
      A.proj <- inla.spde.make.A(mesh, grid.loc)
      Xd_proj = matrix(1, nrow(grid.loc), 1)
      Xg_proj = matrix(1, nrow(grid.loc), 1)
    } else {
      stop ("Need location information for a spatial model to make projections")
    }
  }
  if(!is.null(loc)& is.null(mesh)){
    mesh <- inla.mesh.create(loc)
    n.v <- mesh$n
    warning("Building simple spatial mesh. If using the SPDE-FEM GMRF method,
            the simple mesh may result in spatial bias. Consider bulding a
            more appropriate mesh using INLA::meshbuilder()")
    if(is.null(grid.df)){
      #create dummy projection grid - ##! be sure to turn off projection in reporting out
      grid.loc <- as.matrix(expand.grid(x = seq(0,1,.1), y = seq(0,1,.1)))
      A.proj <- inla.spde.make.A(mesh, grid.loc)
      Xd_proj = matrix(1, nrow(grid.loc), 1)
      Xg_proj = matrix(1, nrow(grid.loc), 1)
    } else {
      grid.loc <- as.matrix(grid.df@coords)
      if(class(grid.df) == "SpatialPoints"){ ##!is this the best way to do this?
        Xd_proj = matrix(1, nrow(grid.loc), 1)
        Xg_proj = matrix(1, nrow(grid.loc), 1)
      } else {
        grid.data <- grid.df@data
        Xd_proj <- as.matrix(grid.data[expert.pred.names])
        Xg_proj <- as.matrix(grid.data[gating.pred.names])
      }
    }
  }
  if( (!is.null(loc)) & (!is.null(mesh)) ){
    n.v <- mesh$n
    A <- inla.spde.make.A(mesh, loc)
    if(is.null(grid.df)){
      grid.loc <- as.matrix(expand.grid(x = seq(0,1,.1), y = seq(0,1,.1)))
      A.proj <- inla.spde.make.A(mesh, grid.loc)
      Xd_proj = matrix(1, nrow(grid.loc), 1)
      Xg_proj = matrix(1, nrow(grid.loc), 1)
    } else {
      grid.loc <- as.matrix(grid.df@coords)
      if(class(grid.df) == "SpatialPoints"){ ##!is this the best way to do this?
        Xd_proj = matrix(1, nrow(grid.loc), 1)
        Xg_proj = matrix(1, nrow(grid.loc), 1)
      } else {
        grid.data <- grid.df@data
        Xd_proj <- as.matrix(grid.data[expert.pred.names])
        Xg_proj <- as.matrix(grid.data[gating.pred.names])
      }
    }
  }
  if(is.null(loc)&!is.null(mesh)){
    stop ("Need to provide locations for observations in spatial list")
  }

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
  if(family[[1]] == "gaussian") Dat$family <- 0
  if(family[[1]] == "Tweedie") Dat$family <- 700
  if(family[[2]] == "log") Dat$link <- 0
  if(family[[2]] == "identity") Dat$link <- 5
  #  if( ncol(response) == 1 | fixStruct == "SFA") Dat$method <- 10
  # if( ncol(response) > 1 & fixStruct != "SFA" ) Dat$method <- 11
  #if( ncol(response) == 1 & fixStruct == "SFA") stop("SFA is a multivariate method and cannot be used in a univariate analysis")
  if(fixStruct == "VVV"){
    Dat$fixStruct <- 35
    n.f <- n.j
  }
  if(fixStruct == "SFA"){
    Dat$fixStruct <- 36
    Dat$method <- 10
    n.f <- n.f
  }
  if(fixStruct == "E"){
    Dat$fixStruct <- 20
    Dat$method <- 10
  }

  Dat$covstruct = covstruct

  Dat$spde <- spdeStruct(mesh)
  return(Dat)
}
