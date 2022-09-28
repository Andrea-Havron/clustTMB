stopifnot(
  require("testthat"),
  require("clustTMB"),
  require("TMB")
)

context("unit test utils-fit.tmb")
runExample('simple')
sdr <- sdreport(obj)
report <- obj$report(obj$env$last.par.best)

#from TMB's simple example
## Test data
set.seed(123)
y <- rep(1900:2010,each=2)
year <- factor(y)
quarter <- factor(rep(1:4,length.out=length(year)))
period <- factor((y > mean(y))+1)
## Random year+quarter effect, fixed period effect:
B <- model.matrix(~year+quarter-1)
A <- model.matrix(~period-1)
B <- as(B,"dgTMatrix")
A <- as(A,"dgTMatrix")
u <- rnorm(ncol(B)) ## logsdu=0
beta <- rnorm(ncol(A))*100
eps <- rnorm(nrow(B),sd=1) ## logsd0=0
x <- as.numeric( A %*% beta + B %*% u + eps )

args <- list(
  data=list(x=x, B=B, A=A),
  parameters=list(u=u*0, beta=beta*0, logsdu=1, logsd0=1),
  random="u",
  DLL="simple",
  silent=TRUE
)

test_that("run.model = FALSE", {
  clustTMB.mod <- fit.tmb(args, 
                          control = list(run.model = FALSE))
  
  expect_equal(names(clustTMB.mod),
               c("obj", "inits", "init.report"))
  expect_equal(obj$par, clustTMB.mod$obj$par)
  expect_equal(obj$fn(obj$par), clustTMB.mod$obj$fn())
  expect_equal(obj$gr(obj$par), clustTMB.mod$obj$gr())
  expect_equal(
    obj$report(c(u=u*0, beta=beta*0, logsdu=1, logsd0=1)), 
    clustTMB.mod$init.report
  )
  
})

test_that("run.model = TRUE, sdr.report = FALSE", {
  
  clustTMB.mod <- fit.tmb(args, 
                          control = list(
                            run.model = TRUE,
                            do.sdreport = FALSE))

  expect_equal(names(clustTMB.mod),
               c("obj", "opt", "report"))
  expect_equal(obj$par, clustTMB.mod$obj$par)
  expect_equal(obj$env$last.par.best, 
               clustTMB.mod$obj$env$last.par.best)
  expect_equal(obj$fn(opt$par)[1], 
               clustTMB.mod$obj$fn(clustTMB.mod$opt$par)[1])
  expect_equal(obj$gr(opt$par)[1], 
               clustTMB.mod$obj$gr(clustTMB.mod$opt$par)[1])
  expect_equal(opt$par, clustTMB.mod$opt$par)
  expect_equal(opt, clustTMB.mod$opt)
  expect_equal(report, clustTMB.mod$report)
  
})

test_that("run.model = TRUE, sdr.report = TRUE", {
  
  clustTMB.mod <- fit.tmb(args, 
                          control = list(
                            run.model = TRUE,
                            do.sdreport = TRUE))
  
  expect_equal(names(clustTMB.mod),
               c("obj", "opt", "report", "sdr"))
  expect_equal(obj$par, clustTMB.mod$obj$par)
  expect_equal(obj$env$last.par.best, 
               clustTMB.mod$obj$env$last.par.best)
  expect_equal(obj$fn(opt$par)[1], 
               clustTMB.mod$obj$fn(clustTMB.mod$opt$par)[1])
  expect_equal(obj$gr(opt$par)[1], 
               clustTMB.mod$obj$gr(clustTMB.mod$opt$par)[1])
  expect_equal(opt$par, clustTMB.mod$opt$par)
  expect_equal(opt, clustTMB.mod$opt)
  expect_equal(report, clustTMB.mod$report)
  expect_equal(sdr, clustTMB.mod$sdr)
  
})
