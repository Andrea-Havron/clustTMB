## Functions used to setup spatial structure

#' Generate SPDE list needed as imput to FEM-SPDE spatial TMB models
#'
#' @param mesh INLA mesh object generated from inla.mesh.create or inla.mesh.2d
#'
#' @return spde List used as input into TMB anisotropic model
#'
#' @importFrom INLA inla.spde2.matern inla.mesh.create
#' @importFrom methods as
#' @importFrom Matrix diag
#' @keywords internal
spdeStruct <- function(mesh){
  spde <-inla.spde2.matern(mesh)
  # ---------- Begin code that prepares object for anisotropy.
  Dset <- 1:2
  # Triangle info
  TV <- mesh$graph$tv           # Triangle to vertex indexing
  V0 <- mesh$loc[TV[,1],Dset]   # V = vertices for each triangle
  V1 <- mesh$loc[TV[,2],Dset]
  V2 <- mesh$loc[TV[,3],Dset]
  E0 <- V2 - V1                      # E = edge for each triangle
  E1 <- V0 - V2
  E2 <- V1 - V0
  # Calculate Areas
  TmpFn <- function(Vec1, Vec2) abs(det( rbind(Vec1, Vec2) ))
  Tri_Area <- rep(NA, nrow(E0))
  for(t in 1:length(Tri_Area)) Tri_Area[t] <- TmpFn( E0[t,],E1[t,] )/2   # T = area of each triangle
  # ---------- End code that prepare objects for anisotropy.
  spde.list <- list(
    "n_s"      = spde$n.spde,
    "n_tri"    = nrow(TV),
    "Tri_Area" = Tri_Area,
    "E0"       = E0,
    "E1"       = E1,
    "E2"       = E2,
    "TV"       = TV - 1,
    "G0"       = spde$param.inla$M0,
    "G0_inv"   = as(diag(1/diag(spde$param.inla$M0)), "dgTMatrix"))

  return(spde.list)
}
