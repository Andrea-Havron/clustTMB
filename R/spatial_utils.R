## Functions used to setup spatial structure

#' Generate SPDE list needed as input to FEM-SPDE spatial TMB models
#'
#' @param mesh mesh object generated from [fmesher::fm_rcdt_2d()] or [fmesher::fm_mesh_2d()]
#'
#' @return spde List used as input into TMB anisotropic model
#'
#' @importFrom Matrix diag
#' @keywords internal
spdeStruct <- function(mesh) {
  if (is.null(mesh)) {
    spde.list <- list(
      "n_s"      = 1,
      "n_tri"    = 0,
      "Tri_Area" = matrix(0, 1, 2),
      "E0"       = matrix(0, 1, 2),
      "E1"       = matrix(0, 1, 2),
      "E2"       = matrix(0, 1, 2),
      "TV"       = matrix(0, 1, 3),
      "G0"       = fmesher::fm_as_dgTMatrix(matrix(0, 2, 2)),
      "G0_inv"   = fmesher::fm_as_dgTMatrix(matrix(0, 2, 2))
    )
  } else {
    spde <- fmesher::fm_fem(mesh)
    # ---------- Begin code that prepares object for anisotropy.
    Dset <- 1:2
    # Triangle info
    TV <- mesh$graph$tv # Triangle to vertex indexing
    V0 <- mesh$loc[TV[, 1], Dset] # V = vertices for each triangle
    V1 <- mesh$loc[TV[, 2], Dset]
    V2 <- mesh$loc[TV[, 3], Dset]
    E0 <- V2 - V1 # E = edge for each triangle
    E1 <- V0 - V2
    E2 <- V1 - V0
    # Calculate Areas
    TmpFn <- function(Vec1, Vec2) abs(det(rbind(Vec1, Vec2)))
    Tri_Area <- rep(NA, nrow(E0))
    # T = area of each triangle
    for (t in seq_along(Tri_Area)) Tri_Area[t] <- TmpFn(E0[t, ], E1[t, ]) / 2
    # ---------- End code that prepare objects for anisotropy.
    spde.list <- list(
      "n_s"      = mesh$n,
      "n_tri"    = nrow(TV),
      "Tri_Area" = Tri_Area,
      "E0"       = E0,
      "E1"       = E1,
      "E2"       = E2,
      "TV"       = TV - 1,
      "G0"       = spde$c0,
      "G0_inv"   = fmesher::fm_as_dgTMatrix(diag(1 / diag(spde$c0)))
    )
  }
  return(spde.list)
}


#' Setup spatial data for mkDat
#'
#' @param n.i number of observations
#' @param spatial.list list of spatial locations and mesh
#' @param projection.dat Points class of projection coordinates or Spatial Points Dataframe containing projection coordinates and projection covariates
#'
#' @importFrom methods is
#'
#' @return list of spatial mesh and sparse A matrix
#' @keywords internal
setup.spatialDat <- function(n.i, spatial.list, projection.dat) {
  loc <- spatial.list$loc
  if (!is.null(loc)) {
    # convert from sp to sf if sp type
    if (!is(loc, "sf")) {
      loc <- sf::st_as_sf(loc)
    }
    Loc <- sf::st_coordinates(loc)
  }
  mesh <- spatial.list$mesh

  if ((is.null(mesh) & is.null(loc)) & !is.null(projection.dat)) {
    warning("loc and mesh are null. Need to provide locations or mesh in spatial.list to initiate spatial model for spatial predictions")
  }
  if ((!is.null(mesh) | !is.null(loc)) & is.null(projection.dat)) {
    warning("spatial projection is turned off. Need to provide locations in projection.list$grid.df for spatial predictions")
  }
  if (!is.null(loc) & is.null(mesh)) {
    # default mesh - in future add options to include arguments for fmesher::fm_mesh_2d
    # for now, user can supply mesh if a more complex mesh is needed
    mesh <- fmesher::fm_rcdt_2d(Loc)
    warning("Building simple spatial mesh. If using the SPDE-FEM GMRF method,
            the simple mesh may result in spatial bias. Consider bulding a
            more appropriate mesh")
  }
  if (is.null(loc) & !is.null(mesh)) {
    if (is.null(mesh$idx$loc)) {
      # if user-supplied mesh built without observation locations, the user must also provide observation locations to the loc argument of spatial.list
      stop("Need to provide locations for observations in spatial.list$loc")
    } else {
      # if user-supplied mesh built using observation locations, these can be obtained from the mesh if not provided
      Loc <- mesh$loc[mesh$idx$loc, 1:2]
    }
  }
  if (is.null(mesh)) {
    A <- fmesher::fm_as_dgCMatrix(matrix(0, n.i, 1))
  }
  if (!is.null(mesh)) {
    A <- fmesher::fm_basis(mesh, Loc)
  }

  out <- list(A = A, mesh = mesh)
  return(out)
}

#' Setup projection data for mkDat
#'
#' @param mesh spatial constrained Delaunay triangulation derived from the fmesher R package
#' @param projection.dat Spatial Points class of projection coordinates or Spatial Points Dataframe containing projection coordinates and projection covariates
#' @param expert.formula Formula defining expert model. This formula corresponds to the covariates included in the response densities. Defaults to intercept only (~1) when no covariates are used.
#' @param gating.formula Formula defining gating model. This formula corresponds to the covariates included in the mixing proportions (logistic regression). Defaults to intercept only (~1) when no covariates are used. When a random effects term is included in the gating network, this formula will be updated so that the intercept term is removed.
#'
#' @return list of projection data
#' @keywords internal
setup.projDat <- function(mesh, projection.dat,
                          expert.formula,
                          gating.formula) {
  grid.df <- projection.dat

  if (is.null(grid.df)) {
    Xd_proj <- matrix(1)
    Xg_proj <- matrix(1)
    doProj <- FALSE
    A_proj <- fmesher::fm_as_dgCMatrix(matrix(0))
  } else {
    if (!is(grid.df, "sf")) {
      grid.df <- sf::st_as_sf(grid.df)
    }
    grid.loc <- sf::st_coordinates(grid.df)
    df <- sf::st_drop_geometry(grid.df)
    if (ncol(grid.df) == 0) {
      df <- data.frame(x = rep(1, nrow(grid.loc)))
    }

    Xd_proj <- model.matrix(lme4::nobars(expert.formula), df)
    Xg_proj <- model.matrix(lme4::nobars(gating.formula), df)
    doProj <- TRUE
    A_proj <- fmesher::fm_basis(mesh, grid.loc)
  }

  out <- list(
    doProj = doProj,
    Xd_proj = Xd_proj,
    Xg_proj = Xg_proj,
    A_proj = A_proj
  )

  return(out)
}
