#' Helper function to run TMB code and return output
#' @importFrom TMB openmp MakeADFun sdreport
#'
#' @param obj.args List of arguments for TMB MakeADFun() function
#' @param opt.args List of arguments for nlminb() function
#' @param control List controlling model runs and standard error reporting
#'
#' @return Fitted objective function, nlminb output, reported values from model, sdreport if true
#'
#' @noRd
fit.tmb <- function(obj.args, opt.args = list(
                      control = list(iter = 800, eval = 800),
                      hessian = NULL, scale = 1,
                      lower = -Inf, upper = Inf
                    ),
                    control = list(run.model = TRUE, do.sdreport = TRUE)) {
  obj <- do.call(MakeADFun, obj.args)
  if (control$run.model) {
    opt <- with(obj, do.call(nlminb, c(list(par, fn, gr), opt.args)))
    report <- obj$report(obj$env$last.par.best)
    if (control$do.sdreport) {
      sdr <- sdreport(obj)
      fit.results <- list(obj = obj, opt = opt, report = report, sdr = sdr)
    } else {
      fit.results <- list(obj = obj, opt = opt, report = report)
    }
  } else {
    fit.results <- list(obj = obj, inits = obj$par, init.report = obj$report())
  }
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
mkFac <- function(d, f = NA) {
  ans <- factor(f)
  if (length(d) > 1) {
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
mkMap <- function(Family, covstruct, rrStruct, reStruct, dim.list, map.ops = NULL) {
  # list2env(dim.list, environment(mkMap))
  n.i <- dim.list$n.i
  n.j <- dim.list$n.j
  n.t <- dim.list$n.t
  n.g <- dim.list$n.g
  n.f.sp <- dim.list$n.f.sp
  n.f.rand <- dim.list$n.f.rand
  n.v <- dim.list$n.v
  nl.rand <- n.j * n.f.rand - (n.f.rand * (n.f.rand - 1)) / 2
  nl.sp <- n.j * n.f.sp - (n.f.sp * (n.f.sp - 1)) / 2
  nl.fix <- ifelse(n.j > 1, (n.j^2 - n.j) / 2, 1)

  # Map out parameters based on model structure
  Map <- list()
  # if(sum(reStruct[1,]>0)){ ##! Rather than map out intercept, change formula to remove intercept and include warning
  #   Map$betag = mkFac(d = c(1,n.g-1), f = rep(NA, n.g-1)) #map out betag if random effects terms
  # }
  if (Family != 700) Map$thetaf <- mkFac(d = c(n.j, n.g), f = rep(NA, n.j * n.g))
  if (covstruct != 30) { # if not a general covariance structure
    Map$logit_corr_fix <- mkFac(d = c(nl.fix, n.g), f = rep(NA, nl.fix * n.g))
  }
  if (rrStruct[1] == 0) {
    Map$ld_rand <- mkFac(d = c(nl.rand, n.g), f = rep(NA, nl.rand * n.g))
  }
  if (rrStruct[2] == 0) {
    Map$ld_sp <- mkFac(d = c(nl.sp, n.g), f = rep(NA, nl.sp * n.g))
  }
  if (covstruct == "EII") { # Diagonal, equal variance within and between cluster
    Map$theta <- mkFac(d = c(n.j, n.g), f = rep(1, n.j * n.g))
  }
  if (covstruct == "VII") { # Diagonal, equal variance within cluster, variable variance between clusters
    Map$theta <- mkFac(
      d = c(n.j, n.g),
      f = as.vector(
        matrix(rep(1:n.g, each = n.j), n.j, n.g)
      )
    )
  }
  if (covstruct == "EEI") { # Diagonal, variable variance within cluster, equal variance between clusters
    Map$theta <- mkFac(
      d = c(n.j, n.g),
      f = as.vector(
        matrix(rep(1:n.j, each = n.g), n.j, n.g, byrow = TRUE)
      )
    )
  }
  if (covstruct == "EEE") { # General, variable variance and covariance within cluster, equal variance and covariance between clusters
    Map$logit_corr_fix <- mkFac(
      d = c(nl.fix, n.g),
      f = as.vector(
        matrix(rep(1:nl.fix, each = n.g),
          nl.fix, n.g,
          byrow = TRUE
        )
      )
    )
    Map$theta <- mkFac(
      d = c(n.j, n.g),
      f = as.vector(
        matrix(rep(1:n.j, each = n.g),
          n.j, n.g,
          byrow = TRUE
        )
      )
    )
  }
  if (reStruct[1, 1] == 0) {
    Map$Hg_input <- mkFac(d = c(2, n.g - 1), f = rep(NA, 2 * (n.g - 1)))
    Map$ln_kappag <- mkFac(d = n.g - 1, f = rep(NA, n.g - 1))
    Map$Gamma_vg <- mkFac(d = c(n.v, n.g - 1), f = rep(NA, n.v * (n.g - 1)))
  }
  if (reStruct[1, 1] > 0) {
    Map$Hg_input <- mkFac(d = c(2, n.g - 1), f = rep(NA, 2 * (n.g - 1)))
    Map$ln_kappag <- mkFac(d = n.g - 1, f = rep(1, n.g - 1))
  }
  if (reStruct[2, 1] == 0) {
    Map$Hd_input <- mkFac(d = c(2, n.f.sp, n.g), f = rep(NA, 2 * n.f.sp * n.g))
    Map$ln_kappad <- mkFac(d = c(n.f.sp, n.g), f = rep(NA, n.f.sp * n.g))
    Map$ln_taud <- mkFac(d = c(n.f.sp, n.g), f = rep(NA, n.f.sp * n.g))
    Map$Omega_vfg <- mkFac(d = c(n.v, n.f.sp, n.g), f = rep(NA, n.v * n.f.sp * n.g))
  }
  if (reStruct[2, 1] > 0) {
    # use when anisotropy turned on Map$Hd_input <- mkFac(d = c(2,n.j, n.g), f = rep(rep(c(1,2),n.j),n.g))
    Map$Hd_input <- mkFac(d = c(2, n.f.sp, n.g), f = rep(NA, 2 * n.f.sp * n.g))
    Map$ln_kappad <- mkFac(d = c(n.f.sp, n.g), f = rep(1, n.f.sp * n.g))
    Map$ln_taud <- mkFac(d = c(n.f.sp, n.g), f = rep(1, n.f.sp * n.g))
  }
  if (rrStruct[2] == 1) {
    Map$ln_taud <- mkFac(d = c(n.f.sp, n.g), f = rep(NA, n.f.sp * n.g))
  }
  if (reStruct[1, 2] == 0) {
    Map$logit_rhog <- mkFac(d = c(n.g - 1), f = rep(NA, n.g - 1))
    Map$ln_sigmaup <- mkFac(d = c(n.g - 1), f = rep(NA, n.g - 1))
    Map$upsilon_tg <- mkFac(d = c(n.t, n.g - 1), f = rep(NA, n.t * (n.g - 1)))
  }
  if (reStruct[2, 2] == 0) {
    Map$logit_rhod <- mkFac(d = c(n.j, n.g), f = rep(NA, n.j * n.g))
    Map$ln_sigmaep <- mkFac(d = c(n.j, n.g), f = rep(NA, n.j * n.g))
    Map$epsilon_tjg <- mkFac(d = c(n.t, n.j, n.g), f = rep(NA, n.t * n.j * n.g))
  }
  # if(reStruct[1,3] == 0){
  #   Map$ln_sigmau <- mkFac(d = c(n.g-1), f = rep(NA, n.g-1))
  #   Map$u_ig <- mkFac(d = c(n.i,n.g-1), f = rep(NA, n.i*(n.g-1)))
  # }
  if (reStruct[2, 3] == 0 | rrStruct[1] == 1) {
    Map$ln_sigmav <- mkFac(d = c(n.f.rand, n.g), f = rep(NA, n.f.rand * n.g))
    if (rrStruct[1] == 0) {
      Map$v_ifg <- mkFac(d = c(n.i, n.f.rand, n.g), f = rep(NA, n.i * n.f.rand * n.g))
    }
  }
  # if(rrStruct[1] == 1){
  #   Map$ln_sigmav <- mkFac(d = c(n.j, n.g), f = rep(NA, n.j*n.g))
  # }
  return(Map)
}


#' Tweedie family and link specification
#'
#' @param link link function association with family
#' @importFrom stats make.link
#' @export
#' @examples
#' fam <- tweedie()
#' fam$family
#' fam$link
tweedie <- function(link = "log") {
  r <- list(family = "tweedie")
  f <- c(r, list(link = link), make.link(link))
  class(f) <- "family"
  return(f)
}

#' Lognormal family and link specification
#'
#' @param link link function association with family
#' @importFrom stats make.link
#' @export
#' @examples
#' fam <- lognormal()
#' fam$family
#' fam$link
lognormal <- function(link = "identity") {
  r <- list(family = "lognormal")
  f <- c(r, list(link = link), make.link(link))
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
#' @return Data list for input into TMB::MakeADFun
#' @keywords internal
#' @noRd

mkDat <- function(response, time.vector, expert.dat, gating.dat,
                  family = gaussian(link = "identity"), ll.method,
                  fixStruct, rrStruct, reStruct, dim.list,
                  offset = NULL,
                  spatial.list = list(loc = NULL, mesh = NULL),
                  projection.list = list(
                    grid.df = NULL, ## !Need more rules about grid.df spatial structure
                    ## !need to ensure order of covariates preserved somehow
                    ## ?match names from expert.dat/gating.dat to grid.df?
                    expert.pred.names = NULL,
                    gating.pred.names = NULL
                  )) {
  n.i <- dim.list$n.i
  n.j <- dim.list$n.j
  n.t <- dim.list$n.t
  n.g <- dim.list$n.g
  n.f.rand <- dim.list$n.f.rand
  n.f.sp <- dim.list$n.f.sp
  n.v <- dim.list$n.v
  # list2env(dim.list, environment(mkDat)) ##! environment locked. try new.env?
  loc <- spatial.list$loc
  if (!is.null(loc)) {
    loc <- loc@coords
  }
  mesh <- spatial.list$mesh
  # list2env(spatial.list, environment(mkDat))
  grid.df <- projection.list$grid.df
  expert.pred.names <- projection.list$expert.pred.names
  gating.pred.names <- projection.list$gating.pred.names
  # list2env(projection.list)
  if ((is.null(mesh) & is.null(loc)) & !is.null(grid.df)) {
    warning("loc and mesh are null. Need to provide locations or mesh in spatial.list to initiate spatial model for spatial predictions")
  }
  if ((!is.null(mesh) | !is.null(loc)) & is.null(grid.df)) {
    warning("spatial projection is turned off. Need to provide locations in projection.list$grid.df for spatial predictions")
  }
  if (!is.null(loc) & is.null(mesh)) {
    # default mesh - in future add options to include arguments for inla.mesh.2d
    # for now, user can supply mesh if a more complex mesh is needed
    if (!requireNamespace("INLA", quietly = TRUE)) {
      stop("INLA must be installed to build a spatial mesh.")
    }
    mesh <- INLA::inla.mesh.create(loc@coords)
    warning("Building simple spatial mesh. If using the SPDE-FEM GMRF method,
            the simple mesh may result in spatial bias. Consider bulding a
            more appropriate mesh using [INLA::meshbuilder()]")
  }
  if (is.null(loc) & !is.null(mesh)) {
    if (is.null(mesh$idx$loc)) {
      # if user-supplied mesh built without observation locations, the user must also provide observation lovations
      stop("Need to provide locations for observations in spatial.list$loc")
    } else {
      # if user-supplied mesh built using observation locations, these can be obtained from the mesh if not provided
      loc <- mesh$loc[mesh$idx$loc, 1:2]
    }
  }
  if (is.null(mesh)) {
    A <- as(matrix(0, n.i, 1), "dgCMatrix")
    n.v <- 1
  } else {
    if (!requireNamespace("INLA", quietly = TRUE)) {
      stop("INLA must be installed to build a spatial mesh.")
    }
    A <- INLA::inla.spde.make.A(mesh, loc)
    n.v <- mesh$n
  }

  if (is.null(grid.df)) {
    Xd_proj <- matrix(1)
    Xg_proj <- matrix(1)
    doProj <- FALSE
    A.proj <- as(matrix(0), "dgCMatrix")
  } else {
    grid.loc <- as.matrix(grid.df@coords)
    if (class(grid.df) == "SpatialPoints") { ## !is this the best way to do this?
      Xd_proj <- matrix(1, nrow(grid.loc), 1)
      Xg_proj <- matrix(1, nrow(grid.loc), 1)
    } else {
      grid.data <- grid.df@data
      if (is.null(expert.pred.names)) {
        Xd_proj <- matrix(1, nrow(grid.loc), 1)
      } else {
        Xd_proj <- as.matrix(grid.data[expert.pred.names])
      }
      if (is.null(gating.pred.names)) {
        Xg_proj <- matrix(1, nrow(grid.loc), 1)
      } else {
        Xg_proj <- as.matrix(grid.data[gating.pred.names])
      }
    }
    doProj <- TRUE
    if (!requireNamespace("INLA", quietly = TRUE)) {
      stop("INLA must be installed to build a projection grid.")
    }
    A.proj <- INLA::inla.spde.make.A(mesh, grid.loc)
  }


  if (is.null(offset)) offset <- rep(1, dim.list$n.i)

  Dat <- list(
    Y = as.array(response),
    t = time.vector,
    Xd = expert.dat,
    Xg = gating.dat,
    Xpz = matrix(1, n.i, 1),
    Offset = offset,
    A = A,
    A_proj = A.proj,
    doProj = doProj,
    Xd_proj = Xd_proj,
    Xg_proj = Xg_proj
  )

  Dat$spde <- spdeStruct(mesh)
  Dat$family <- .valid_family[family[[1]]]
  Dat$link <- .valid_link[family[[2]]]
  Dat$loglike <- ll.method
  Dat$fixStruct <- fixStruct.lookup(fixStruct)

  Dat$rrStruct <- rrStruct
  Dat$reStruct <- reStruct

  return(Dat)
}

#' mkRandom: set up random effects component of the model
#'
#' @param expertformula Formula defining expert model.
#' @param gatingformula Formula defining gating model.
#' @param expertdata  Data frame containing expert model covariates.
#' @param gatingdata  Data frame containing gating model covariates.
#' @param spatial.list List of data objects needed when fitting a spatial GMRF model
#' @param dim.list Class object containing model dimensions
#'
#' @importFrom lme4 subbars
#'
#' @return list vector containing random effects components of the model
#' @noRd
mkRandom <- function(expertformula, gatingformula, expertdata, gatingdata, spatial.list, dim.list) {
  expert.split <- splitForm(expertformula)
  expert.re.names <- expert.split$reTrmClasses
  gating.split <- splitForm(gatingformula)
  gating.re.names <- gating.split$reTrmClasses

  if (("gmrf" %in% expert.re.names | "gmrfSpeedup" %in% expert.re.names) &
    (is.null(spatial.list$loc) & is.null(spatial.list$mesh))) {
    stop("You have specified a spatal model and need to provide location or
          mesh data in the spatial.list argument")
  }
  if (("gmrf" %in% gating.re.names | "gmrfSpeedup" %in% gating.re.names) &
    (is.null(spatial.list$loc) & is.null(spatial.list$mesh))) {
    stop("You have specified a spatal model and need to provide location or
          mesh data in the spatial.list argument")
  }

  specials <- c("ar1", "gmrf", "gmrfSpeedup")
  if (length(expert.re.names) > 0) {
    for (i in seq_along(expert.re.names)) {
      if (!(expert.re.names[i] %in% specials)) {
        stop("Currently clustTMB only works with spatio-temporal random effects")
      }
    }
  }
  if (length(gating.re.names) > 0) {
    for (i in seq_along(gating.re.names)) {
      if (!(gating.re.names[i] %in% specials)) {
        stop("Currently clustTMB only works with spatio-temporal random effects")
      }
    }
  }

  if ("ar1" %in% expert.re.names) {
    idx <- which(expert.re.names == "ar1")
    ar1.form <- as.formula(paste("~", deparse(expert.split$reTrmFormulas[[idx]])))
    ## ! make sure this enters TMB as numeric - although consider using glmmTMB factor approach
    expert.time <- model.frame(subbars(ar1.form), expertdata)
  } else {
    expert.time <- rep(1, dim.list$n.i)
  }

  ## TODO: section below unused
  # if ("ar1" %in% gating.re.names) {
  #   idx <- which(gating.re.names == "ar1")
  #   ar1.form <- as.formula(paste("~", deparse(gating.split$reTrmFormulas[[idx]])))
  #   gating.time <- model.frame(subbars(ar1.form), gatingdata)
  # } else {
  #   gating.time <- rep(1, dim.list$n.i)
  # }
  #
  # if (("gmrf" %in% expert.re.names) | ("gmrfSpeedup" %in% expert.re.names)) {
  #   loc.id <- spatial.list$loc@coords[,1] + spatial.list$loc@coords[,2]
  #   expertdata$loc <- as.numeric(factor())
  #   idx <- which((expert.re.names == "gmrf") | (expert.re.names == "gmrfSpeedup"))
  #   gmrf.form <- as.formula(paste("~", deparse(expert.split$reTrmFormulas[[idx]])))
  #   expert.gmrf <- model.frame(subbars(gmrf.form), expertdata)
  # } else {
  #   expert.gmrf <- NA
  # }
  #
  # if (("gmrf" %in% gating.re.names) | ("gmrfSpeedup" %in% gating.re.names)) {
  #   loc.id <- spatial.list$loc@coords[,1] + spatial.list$loc@coords[,2]
  #   gatingdata$loc <- as.numeric(factor(loc.id))
  #   idx <- which((gating.re.names == "gmrf") | (gating.re.names == "gmrfSpeedup"))
  #   gmrf.form <- as.formula(paste("~", deparse(gating.split$reTrmFormulas[[idx]])))
  #   gating.gmrf <- model.frame(subbars(gmrf.form), gatingdata)
  # } else {
  #   gating.gmrf <- NA
  # }
  ### end of unused section

  ## TODO: expert/gating .ar1 and .gmrf provide information about interactions between space/time
  ## TODO: Spatio-temporal interactions not implemented yet!!
  reStruct <- matrix(0, 2, 3)
  random.names <- c()

  for (i in seq_along(gating.re.names)) {
    if (gating.re.names[i] == "gmrf") {
      reStruct[1, 1] <- 3
      random.names <- c(random.names, "Gamma_vg")
    }
    if (gating.re.names[i] == "gmrfSpeedup") {
      reStruct[1, 1] <- 4
      random.names <- c(random.names, "Gamma_vg")
    }
    if (gating.re.names[i] == "ar1") {
      reStruct[1, 2] <- 2
      random.names <- c(random.names, "upsilon_tg")
    }
    # if(gating.re.names[i] == "norm"){
    #   reStruct[1,3] <- 1
    #   random.names <- c(random.names, "u_ig")
    # }
  }
  for (i in seq_along(expert.re.names)) {
    if (expert.re.names[i] == "gmrf") {
      reStruct[2, 1] <- 3
      random.names <- c(random.names, "Omega_vfg")
    }
    if (expert.re.names[i] == "gmrfSpeedup") {
      reStruct[2, 1] <- 4
      random.names <- c(random.names, "Omega_vfg")
    }
    if (expert.re.names[i] == "ar1") {
      reStruct[2, 2] <- 2
      random.names <- c(random.names, "epsilon_tjg")
    }
    # if(expert.re.names[i] == "norm"){
    #   reStruct[2,3] <- 1
    #   random.names <- c(random.names, "v_ifg")
    # }
  }

  out <- list(
    reStruct = reStruct, 
    random.names = random.names, 
    expert.time = expert.time
    )
  return(out)
}


#' Fixed Covariance Structure names
#' @export
#' @examples
#' fixStruct.names()
fixStruct.names <- function() {
  return(c("E", "V", "EII", "VII", "EEI", "VVI", "VVV", "EEE"))
}

#' Names of parameters with initial values that can be modified
#' @noRd
start.names <- function() {
  return(c(
    "thetaf", "ln_kappa_g", "ln_kappa_d", "ln_tau_d", "logit_rhog",
    "logit_rhod", "ln_sigmaup", "ln_sigma_ep", "ln_sigmav",
    "upsilon_tg", "epsilon_tjg", "v_ifg", "Gamma_vg", "Omega_vfg"
  ))
} # removed ln_sigmau, u_ig

#' Parameter Information
#'
#' @return Description of parameters, including dimension and structure
#' @export
#'
#' @examples
#' parm.lookup()
parm.lookup <- function() {
  df <- data.frame(
    parm = c(
      "betag", "betad", "betapz", "theta", "thetaf", "logit_corr_fix",
      "ld_rand", "ld_sp", "Hg_input", "Hd_input", "ln_kappa_g", "ln_kappa_d",
      "ln_tau_d", "logit_rhog", "logit_rhod", "ln_sigmaup", "ln_sigmaep",
      "ln_sigmav", "upsilon_tg", "epsilon_tjg", "v_ifg",
      "Gamma_vg", "Omega_vfg"
    ), # removed ln_sigmau, u_ig
    # Fixed 19 -> 18; Random 6 -> 5
    type = c(rep("Fixed", 18), rep("Random", 5)), 
    str = c(
      "Matrix", "Array", "Array", rep("Matrix", 6), "Array",
      "Vector", "Matrix", "Matrix", "Vector", "Matrix", "Vector",
      # 'Matrix', 'Vector', 'Matrix', rep('Array', 6)
      "Matrix", "Matrix", rep("Array", 5)
    ),
    dim = c(
      "Kg,G-1", "Kd,J,G", "M,J,G", "J,G", "J,G", "(J^2-J)/2", "Fr,G", "Fs,G",
      "2,G-1", "2,J,G", "G-1", "J,G", "J,G", "G-1", "J,G", "G-1", "J,G", #' G-1'
      # 'J,G', 'T,G-1', 'T,J,G',  'N,G-1', 'N,J/Fr,G', 'N,G-1', 'V,J/Fs,G'
      "Fr,G", "T,G-1", "T,J,G", "N,J/Fr,G", "N,G-1", "V,J/Fs,G"
    ),
    descr = c(
      "Gating covariate coefficients",
      "Expert covariate coefficients",
      "Zero-inflation covariate coefficients",
      "Observation variance terms (natural log)",
      "Tweedie power parameter (adjusted logit)",
      "Observation multivariate correlation (logit)",
      "Overdispersion rank reduction loadings",
      "Spatial rank reduction loadings",
      "Anisotropy terms in spatial gating",
      "Anisotropy terms in spatial expert",
      "Spatial gating decay (natural log)",
      "Spatial expert decay (natural log)",
      "Spatial expert precision (natural log)",
      "Temporal gating correlation (logit)",
      "Temporal expert correlation (logit)",
      "Temporal gating variance (natural log)",
      "Temporal expert variance (natural log)",
      # 'Overdispersion gating standard deviation (natural log)',
      "Overdispersion expert standard deviation (natural log)",
      "Temporal gating random effect",
      "Temporal expert random effect",
      # 'Overdispersion gating random effect',
      "Overdispersion expert random effect",
      "Spatial gating random effect",
      "Spatial expert random effect"
    )
  )
  key <- data.frame(
    dim = c(
      "Fr", "Fs", "G", "J", "Kd", "Kg", "M", "N", "T", "V"
    ),
    descr = c(
      "Rank reduction on random error",
      "Rank reduction on spatial effects",
      "Number of clusters",
      "Number of columns in the response",
      "Number of covariates in gating",
      "Number of covariates in expert",
      "Number of covariates in zero inflation",
      "Number of observations in the response",
      "Number of unique temporal units",
      "Number of vertices in INLA mesh"
    )
  )
  note <- list("ln_sigma_v not estimated - fixed to ln(1)")
  out <- list(parm = df, key = key, note = note)
  return(out)
}


#' Calculates skewness
#'
#' @param x numeric vector of values for which skewness is calculated
#'
#' @return skewness value of x
#' @export
#' @examples
#' skewness(rgamma(100, 1, 1))
skewness <- function(x) {
  n <- length(x)
  x <- x - mean(x)
  y <- sqrt(n) * sum(x^3) / (sum(x^2)^(3 / 2))
  y <- y * ((1 - 1 / n))^(3 / 2)
  return(y)
}

#' Check if INLA installed (i.e., not on CRAN)
#'
#' @export
#' @examples
#' inla_installed()
inla_installed <- function() {
  requireNamespace("INLA", quietly = TRUE)
}


#' fixStruct look-up table
#'
#' @param fixStruct user input character vector
#'
#' @return numeric value associated with character vector
#' @examples
#' fixStruct.lookup("E")
fixStruct.lookup <- function(fixStruct) {
  df <- data.frame(
    fixStruct =
      c(
        "E", "V", "EII",
        "VII", "EEI", "VVI",
        "VVV", "EEE"
      ),
    value =
      c(
        10, 10, 20, 20,
        20, 20, 30, 30
      )
  )
  out <- df[df$fixStruct == fixStruct, ]$value
  return(out)
}
