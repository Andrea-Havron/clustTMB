stopifnot(
  require("testthat"),
  require("clustTMB")
)

context("test utils-setup-projDat")


test_that("grid.loc, proj data - sp object", {
  n.i <- 100
  Loc <- matrix(runif(n.i * 2), n.i, 2)
  mesh <- INLA::inla.mesh.create(Loc)
  gating.formula <- ~Xg
  expert.formula <- ~Xd

  proj.grid <- expand.grid(
    x = seq(0, 1, 0.1),
    y = seq(0, 1, 0.1)
  )
  n.proj <- nrow(proj.grid)
  df <- data.frame(
    Xd = rnorm(n.proj),
    Xg = rnorm(n.proj)
  )
  dat <- sp::SpatialPointsDataFrame(
    coords = proj.grid,
    data = df
  )
  Aproj <- INLA::inla.spde.make.A(mesh, dat@coords)

  projDat <- setup.projDat(
    mesh, dat,
    expert.formula,
    gating.formula
  )

  expect_equal(TRUE, projDat$doProj)
  expect_equal(df$Xd, unname(projDat$Xd_proj[, "Xd"]))
  expect_equal(df$Xg, unname(projDat$Xg_proj[, "Xg"]))
  expect_equal(Aproj, projDat$A_proj)
})

test_that("grid.loc, proj data - sf object", {
  n.i <- 100
  Loc <- matrix(runif(n.i * 2), n.i, 2)
  mesh <- INLA::inla.mesh.create(Loc)
  gating.formula <- ~Xg
  expert.formula <- ~Xd

  proj.grid <- expand.grid(
    x = seq(0, 1, 0.1),
    y = seq(0, 1, 0.1)
  )
  Aproj <- INLA::inla.spde.make.A(mesh, as.matrix(proj.grid))
  n.proj <- nrow(proj.grid)

  df <- data.frame(
    x = proj.grid$x,
    y = proj.grid$y,
    Xd = rnorm(n.proj),
    Xg = rnorm(n.proj)
  )

  dat <- sf::st_as_sf(df, coords = c("x", "y"))
  projDat <- setup.projDat(
    mesh, dat,
    expert.formula,
    gating.formula
  )

  expect_equal(TRUE, projDat$doProj)
  expect_equal(df$Xd, unname(projDat$Xd_proj[, "Xd"]))
  expect_equal(df$Xg, unname(projDat$Xg_proj[, "Xg"]))
  expect_equal(Aproj, projDat$A_proj)
})


test_that("grid.loc, no proj data - sp object", {
  n.i <- 100
  Loc <- matrix(runif(n.i * 2), n.i, 2)
  mesh <- INLA::inla.mesh.create(Loc)
  gating.formula <- ~1
  expert.formula <- ~1

  proj.grid <- expand.grid(
    x = seq(0, 1, 0.1),
    y = seq(0, 1, 0.1)
  )
  n.proj <- nrow(proj.grid)
  sp::coordinates(proj.grid) <- ~ x * y
  Aproj <- INLA::inla.spde.make.A(mesh, proj.grid@coords)

  projDat <- setup.projDat(
    mesh, proj.grid,
    expert.formula,
    gating.formula
  )

  expect_equal(TRUE, projDat$doProj)
  expect_equal(rep(1, n.proj), unname(as.vector(projDat$Xd_proj)))
  expect_equal(rep(1, n.proj), unname(as.vector(projDat$Xg_proj)))
  expect_equal(Aproj, projDat$A_proj)
})

test_that("no grid.loc, no proj data", {
  projection.dat <- NULL
  mesh <- NULL

  Xd_proj <- Xg_proj <- matrix(1)
  doProj <- FALSE
  A_proj <- as(matrix(0), "dgCMatrix")

  projDat <- setup.projDat(mesh, projection.dat)

  expect_equal(Xd_proj, projDat$Xd_proj)
  expect_equal(Xg_proj, projDat$Xg_proj)
  expect_equal(doProj, projDat$doProj)
  expect_equal(A_proj, projDat$A_proj)
})

test_that("grid.loc, no proj data - sf object", {
  n.i <- 100
  Loc <- matrix(runif(n.i * 2), n.i, 2)
  mesh <- INLA::inla.mesh.create(Loc)
  gating.formula <- ~1
  expert.formula <- ~1

  proj.grid <- expand.grid(
    x = seq(0, 1, 0.1),
    y = seq(0, 1, 0.1)
  )
  Aproj <- INLA::inla.spde.make.A(mesh, as.matrix(proj.grid))
  n.proj <- nrow(proj.grid)
  proj.grid <- sf::st_as_sf(as.data.frame(proj.grid), coords = c("x", "y"))
  projDat <- setup.projDat(
    mesh, proj.grid,
    expert.formula,
    gating.formula
  )

  expect_equal(TRUE, projDat$doProj)
  expect_equal(rep(1, n.proj), unname(as.vector(projDat$Xd_proj)))
  expect_equal(rep(1, n.proj), unname(as.vector(projDat$Xg_proj)))
  expect_equal(Aproj, projDat$A_proj)
})
