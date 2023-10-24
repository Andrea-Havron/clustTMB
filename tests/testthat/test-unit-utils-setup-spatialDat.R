context("test utils-setup-spatialDat")


test_that("loc and mesh - sp object", {
  spatial.list <- list(
    loc = NULL,
    mesh = NULL,
    init.range = list(
      gating.range = NULL,
      expert.range = NULL
    )
  )

  n.i <- 100
  Loc <- matrix(runif(n.i * 2), n.i, 2)
  mesh1 <- fmesher::fm_rcdt_2d(Loc)
  bnd <- fmesher::fm_nonconvex_hull(Loc)
  mesh2 <- fmesher::fm_rcdt_2d(boundary = bnd)
  A1 <- fmesher::fm_basis(mesh1, Loc)
  A2 <- fmesher::fm_basis(mesh2, Loc)

  spatial.list$loc <- Loc
  spatial.list$mesh <- mesh1

  expect_error(setup.spatialDat(
    n.i,
    spatial.list,
    NULL
  ))
  loc <- data.frame(x = Loc[, 1], y = Loc[, 2])
  sp::coordinates(loc) <- ~ x * y
  spatial.list$loc <- loc

  expect_warning(setup.spatialDat(
    n.i,
    spatial.list,
    NULL
  ))
  grid.df <- 1

  spDat <- setup.spatialDat(
    n.i,
    spatial.list,
    grid.df
  )
  expect_equal(spDat$mesh, mesh1)
  expect_equal(spDat$A, A1)

  spatial.list$mesh <- mesh2
  spDat <- setup.spatialDat(
    n.i,
    spatial.list,
    grid.df
  )
  expect_equal(spDat$mesh, mesh2)
  expect_equal(spDat$A, A2)
})

test_that("loc and mesh - sf object", {
  spatial.list <- list(
    loc = NULL,
    mesh = NULL,
    init.range = list(
      gating.range = NULL,
      expert.range = NULL
    )
  )

  n.i <- 100
  Loc <- matrix(runif(n.i * 2), n.i, 2)
  mesh1 <- fmesher::fm_rcdt_2d(Loc)
  bnd <- fmesher::fm_nonconvex_hull(Loc)
  mesh2 <- fmesher::fm_rcdt_2d(boundary = bnd)
  A1 <- fmesher::fm_basis(mesh1, Loc)
  A2 <- fmesher::fm_basis(mesh2, Loc)


  loc <- data.frame(x = Loc[, 1], y = Loc[, 2])
  spatial.list$loc <- sf::st_as_sf(loc, coords = c("x", "y"))

  expect_warning(setup.spatialDat(
    n.i,
    spatial.list,
    NULL
  ))
  grid.df <- 1
  spatial.list$mesh <- mesh1

  spDat <- setup.spatialDat(
    n.i,
    spatial.list,
    grid.df
  )
  expect_equal(spDat$mesh, mesh1)
  expect_equal(spDat$A, A1)

  spatial.list$mesh <- mesh2
  spDat <- setup.spatialDat(
    n.i,
    spatial.list,
    grid.df
  )
  expect_equal(spDat$mesh, mesh2)
  expect_equal(spDat$A, A2)
})

test_that("loc, no mesh", {
  spatial.list <- list(
    loc = NULL,
    mesh = NULL,
    init.range = list(
      gating.range = NULL,
      expert.range = NULL
    )
  )

  n.i <- 100
  Loc <- matrix(runif(n.i * 2), n.i, 2)
  mesh <- fmesher::fm_rcdt_2d(Loc)
  A <- fmesher::fm_basis(mesh, Loc)
  loc <- data.frame(x = Loc[, 1], y = Loc[, 2])
  spatial.list$loc <- sf::st_as_sf(loc, coords = c("x", "y"))

  grid.df <- 1

  expect_warning(setup.spatialDat(
    n.i,
    spatial.list,
    grid.df
  ))

  spDat <- suppressWarnings(
    setup.spatialDat(
      n.i,
      spatial.list,
      grid.df
    )
  )
  # mesh meta does not match because fmesher::fm_rcdt_2d() called separately
  expect_equal(spDat$mesh[2:8], mesh[2:8])
  expect_equal(spDat$A, A)
})

test_that("no loc, mesh", {
  spatial.list <- list(
    loc = NULL,
    mesh = NULL,
    init.range = list(
      gating.range = NULL,
      expert.range = NULL
    )
  )

  Loc <- matrix(runif(100), 50, 2)
  mesh1 <- fmesher::fm_rcdt_2d(Loc)
  bnd <- fmesher::fm_nonconvex_hull(Loc)
  mesh2 <- fmesher::fm_rcdt_2d(boundary = bnd)
  A1 <- fmesher::fm_basis(mesh1, Loc)
  A2 <- fmesher::fm_basis(mesh2, Loc)

  spatial.list$mesh <- mesh1
  grid.df <- 1

  spDat <- setup.spatialDat(
    n.i,
    spatial.list,
    grid.df
  )
  expect_equal(mesh1, spDat$mesh)
  expect_equal(A1, spDat$A)

  spatial.list$mesh <- mesh2
  expect_error(setup.spatialDat(
    n.i,
    spatial.list,
    grid.df
  ))

  loc <- data.frame(x = Loc[, 1], y = Loc[, 2])
  spatial.list$loc <- sf::st_as_sf(loc, coords = c("x", "y"))
  spDat <- setup.spatialDat(
    n.i,
    spatial.list,
    grid.df
  )

  expect_equal(spDat$mesh[2:8], mesh2[2:8])
  expect_equal(spDat$A, A2)
})

test_that("no loc, no mesh", {
  spatial.list <- list(
    loc = NULL,
    mesh = NULL,
    init.range = list(
      gating.range = NULL,
      expert.range = NULL
    )
  )

  n.i <- 100
  mesh <- NULL
  A <- as(matrix(0, 100, 1), "dgCMatrix")

  spDat <- setup.spatialDat(
    n.i,
    spatial.list,
    NULL
  )

  grid.df <- 1
  expect_equal(mesh, spDat$mesh)
  expect_equal(A, spDat$A)

  expect_warning(setup.spatialDat(
    n.i,
    spatial.list,
    grid.df
  ))
})
