## Functions used to setup spatial structure

#' Generate SPDE list needed as imput to FEM-SPDE spatial TMB models
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

