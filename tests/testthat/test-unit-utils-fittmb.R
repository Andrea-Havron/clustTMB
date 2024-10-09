stopifnot(
  require("testthat"),
  require("clustTMB"),
  require("TMB")
)

context("unit test utils-fit.tmb")


test_that("run.model = FALSE", {
  skip_on_cran()
  skip_if_not(run_model)
  clustTMB.mod <- fit.tmb(args,
    control = list(run.model = FALSE)
  )

  expect_equal(
    names(clustTMB.mod),
    c("obj", "inits", "init.report")
  )
  expect_equal(obj$par, clustTMB.mod$obj$par)
  expect_equal(obj$fn(obj$par), clustTMB.mod$obj$fn())
  expect_equal(obj$gr(obj$par), clustTMB.mod$obj$gr())
  expect_equal(
    obj$report(c(u = u * 0, beta = beta * 0, logsdu = 1, logsd0 = 1)),
    clustTMB.mod$init.report
  )
})

test_that("run.model = TRUE, sdr.report = FALSE", {
  skip_on_cran()
  skip_if_not(run_model)
  clustTMB.mod <- fit.tmb(args,
    control = list(
      run.model = TRUE,
      do.sdreport = FALSE
    )
  )

  expect_equal(
    names(clustTMB.mod),
    c("obj", "opt", "report")
  )
  expect_equal(obj$par, clustTMB.mod$obj$par)
  expect_equal(
    obj$env$last.par.best,
    clustTMB.mod$obj$env$last.par.best
  )
  expect_equal(
    obj$fn(opt$par)[1],
    clustTMB.mod$obj$fn(clustTMB.mod$opt$par)[1]
  )
  expect_equal(
    obj$gr(opt$par)[1],
    clustTMB.mod$obj$gr(clustTMB.mod$opt$par)[1]
  )
  expect_equal(opt$par, clustTMB.mod$opt$par)
  expect_equal(opt, clustTMB.mod$opt)
  expect_equal(report, clustTMB.mod$report)
})

test_that("run.model = TRUE, sdr.report = TRUE", {
  skip_on_cran()
  skip_if_not(run_model)
  clustTMB.mod <- fit.tmb(args,
    control = list(
      run.model = TRUE,
      do.sdreport = TRUE
    )
  )

  expect_equal(
    names(clustTMB.mod),
    c("obj", "opt", "report", "sdr")
  )
  expect_equal(obj$par, clustTMB.mod$obj$par)
  expect_equal(
    obj$env$last.par.best,
    clustTMB.mod$obj$env$last.par.best
  )
  expect_equal(
    obj$fn(opt$par)[1],
    clustTMB.mod$obj$fn(clustTMB.mod$opt$par)[1]
  )
  expect_equal(
    obj$gr(opt$par)[1],
    clustTMB.mod$obj$gr(clustTMB.mod$opt$par)[1]
  )
  expect_equal(opt$par, clustTMB.mod$opt$par)
  expect_equal(opt, clustTMB.mod$opt)
  expect_equal(report, clustTMB.mod$report)
  expect_equal(sdr, clustTMB.mod$sdr)
})

if (run_model) {
  dyn.unload(dynlib(testthat::test_path("ref", "simple")))
  file.remove(dynlib(testthat::test_path("ref", "simple")))
  file.remove(testthat::test_path("ref", "simple.o"))
}
