stopifnot(
  require("testthat"),
  require("clustTMB")
)

context("test utils-mkRandom")

test_that("no random effects", {
  n.i <- 100
  n.j <- 2
  n.g <- 3
  expertformula <- ~1; expertdata <- matrix(1,n.i,1)
  gatingformula <- ~1; gatingdata <- matrix(1,n.i,1)
  spatial.list <- list(loc = NULL, mesh = NULL, init.range = list(gating.range = NULL, expert.range = NULL))
  dim.list <- list(n.i = n.i, n.j = n.j, n.t = 1,
                   n.g = n.g, n.f.sp = n.j, n.f.rand = n.j, n.v = NULL
                   )
  dim.list$nl.fix <- (dim.list$n.j^2 - dim.list$n.j) / 2
  dim.list$nl.rand <- dim.list$n.j * dim.list$n.f.rand - (dim.list$n.f.rand * (dim.list$n.f.rand - 1)) / 2
  dim.list$nl.sp <- dim.list$n.j * dim.list$n.f.sp - (dim.list$n.f.sp * (dim.list$n.f.sp - 1)) / 2
  
  true.re <- list(reStruct = matrix(0,2,3), random.names = NULL, expert.time = rep(1, n.i))
  RE <- clustTMB:::mkRandom(expertformula, gatingformula, expertdata, gatingdata, spatial.list, dim.list)
  expect_equal(true.re, RE)
})

test_that("ar1 in expert", {
  n.i <- 100
  n.j <- 2
  n.g <- 3
  expertformula <- ~ar1(0 + 1|time); expertdata <- data.frame(time = 1:n.i)
  gatingformula <- ~1; gatingdata <- matrix(1,n.i,1)
  spatial.list <- list(loc = NULL, mesh = NULL, init.range = list(gating.range = NULL, expert.range = NULL))
  dim.list <- list(n.i = n.i, n.j = n.j, n.t = n.i,
                   n.g = n.g, n.f.sp = n.j, n.f.rand = n.j, n.v = NULL
  )
  dim.list$nl.fix <- (dim.list$n.j^2 - dim.list$n.j) / 2
  dim.list$nl.rand <- dim.list$n.j * dim.list$n.f.rand - (dim.list$n.f.rand * (dim.list$n.f.rand - 1)) / 2
  dim.list$nl.sp <- dim.list$n.j * dim.list$n.f.sp - (dim.list$n.f.sp * (dim.list$n.f.sp - 1)) / 2
  
  true.re <- list(reStruct = matrix(0,2,3), random.names = "epsilon_tjg", expert.time = data.frame(time = 1:n.i))
  true.re$reStruct[2,2] <- 2
  RE <- clustTMB:::mkRandom(expertformula, gatingformula, expertdata, gatingdata, spatial.list, dim.list)
  expect_equal(true.re$reStruct, RE$reStruct)
  expect_equal(true.re$random.names, RE$random.names)
  expect_equal(true.re$expert.time$time, RE$expert.time$time)
})

test_that("gmrf in expert", {
  n.i <- 100
  n.j <- 2
  n.g <- 3
  set.seed(123)
  Loc <- data.frame(x = runif(n.i), y = runif(n.i))
  sp::coordinates(Loc) <- ~x*y
  expertformula <- ~gmrf(0 + 1|loc); expertdata <-  data.frame(x = rep(1, dim.list$n.i))
  gatingformula <- ~1; gatingdata <- matrix(1,n.i,1)
  spatial.list <- list(loc = Loc, mesh = NULL, init.range = list(gating.range = NULL, expert.range = NULL))
  dim.list <- list(n.i = n.i, n.j = n.j, n.t = n.i,
                   n.g = n.g, n.f.sp = n.j, n.f.rand = n.j, n.v = NULL
  )
  dim.list$nl.fix <- (dim.list$n.j^2 - dim.list$n.j) / 2
  dim.list$nl.rand <- dim.list$n.j * dim.list$n.f.rand - (dim.list$n.f.rand * (dim.list$n.f.rand - 1)) / 2
  dim.list$nl.sp <- dim.list$n.j * dim.list$n.f.sp - (dim.list$n.f.sp * (dim.list$n.f.sp - 1)) / 2
  
  true.re <- list(reStruct = matrix(0,2,3), random.names = "Omega_vfg", expert.time = rep(1,n.i))
  true.re$reStruct[2, 1] <- 3
  RE <- clustTMB:::mkRandom(expertformula, gatingformula, expertdata, gatingdata, spatial.list, dim.list)
  expect_equal(true.re$reStruct, RE$reStruct)
  expect_equal(true.re$random.names, RE$random.names)
  expect_equal(true.re$expert.time, RE$expert.time)
})

test_that("gmrf in gating", {
  n.i <- 100
  n.j <- 2
  n.g <- 3
  set.seed(123)
  Loc <- data.frame(x = runif(n.i), y = runif(n.i))
  sp::coordinates(Loc) <- ~x*y
  expertformula <- ~1; expertdata <- matrix(1,n.i,1)  
  gatingformula <- ~gmrf(0 + 1|loc); gatingdata <- data.frame(x = rep(1, dim.list$n.i))
  spatial.list <- list(loc = Loc, mesh = NULL, init.range = list(gating.range = NULL, expert.range = NULL))
  dim.list <- list(n.i = n.i, n.j = n.j, n.t = n.i,
                   n.g = n.g, n.f.sp = n.j, n.f.rand = n.j, n.v = NULL
  )
  dim.list$nl.fix <- (dim.list$n.j^2 - dim.list$n.j) / 2
  dim.list$nl.rand <- dim.list$n.j * dim.list$n.f.rand - (dim.list$n.f.rand * (dim.list$n.f.rand - 1)) / 2
  dim.list$nl.sp <- dim.list$n.j * dim.list$n.f.sp - (dim.list$n.f.sp * (dim.list$n.f.sp - 1)) / 2
  
  true.re <- list(reStruct = matrix(0,2,3), random.names = "Gamma_vg", expert.time = rep(1,n.i))
  true.re$reStruct[1, 1] <- 3
  RE <- clustTMB:::mkRandom(expertformula, gatingformula, expertdata, gatingdata, spatial.list, dim.list)
  expect_equal(true.re$reStruct, RE$reStruct)
  expect_equal(true.re$random.names, RE$random.names)
  expect_equal(true.re$expert.time, RE$expert.time)
})
