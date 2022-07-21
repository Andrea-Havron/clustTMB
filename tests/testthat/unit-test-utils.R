stopifnot(
  require("testthat"),
  require("clustTMB")
)

context("test fit.tmb")
TMB::runExample("simple")
args <- list(
  data = list(x = x, B = B, A = A),
  parameters = list(u = u * 0, beta = beta * 0, logsdu = 1, logsd0 = 1),
  random = "u",
  DLL = "simple",
  silent = TRUE
)

obj.exp <- do.call(TMB::MakeADFun, args)
test_that("test fit.tmb, run.model = FALSE", {
  obj.obs <- fit.tmb(args, control = list(run.model = FALSE))
  expect_equal(obj.exp$par, obj.obs$obj$par)
  expect_equal(obj.exp$fn(), obj.obs$obj$fn())
  expect_equal(obj.exp$gr(), obj.obs$obj$gr())
  rep.exp <- obj.exp$report()
  expect_equal(rep.exp, obj.obs$init.report)
  expect_equal(c("obj", "inits", "init.report"), names(obj.obs))
})

test_that("test fit.tmb, run.model = TRUE, sdreport = FALSE", {
  opt.obs <- fit.tmb(args, control = list(
    run.model = TRUE,
    do.sdreport = FALSE
  ))
  expect_equal(obj$par, opt.obs$obj$par)
  expect_equal(obj$fn(), opt.obs$obj$fn())
  expect_equal(obj$gr(), opt.obs$obj$gr())
  expect_equal(opt, opt.obs$opt)
  rep.exp <- obj$report()
  expect_equal(rep.exp, opt.obs$report)
  expect_equal(c("obj", "opt", "report"), names(opt.obs))
})

test_that("test fit.tmb, run_model = TRUE, sdreport = TRUE", {
  rep.exp <- obj$report()
  sdr <- sdreport(obj)
  opt.obs <- fit.tmb(args, control = list(
    run.model = TRUE,
    do.sdreport = TRUE
  ))
  expect_equal(obj$par, opt.obs$obj$par)
  expect_equal(obj$fn(), opt.obs$obj$fn())
  expect_equal(obj$gr(), opt.obs$obj$gr())
  expect_equal(opt, opt.obs$opt)
  expect_equal(rep.exp, opt.obs$report)
  expect_equal(sdr, opt.obs$sdr)
  expect_equal(c("obj", "opt", "report", "sdr"), names(opt.obs))
})

context("test fixStruct.lookup")
test_that("test fixStruct.lookup", {
  expect_equal(10, fixStruct.lookup("E"))
  expect_equal(10, fixStruct.lookup("V"))
  expect_equal(20, fixStruct.lookup("EII"))
  expect_equal(20, fixStruct.lookup("VII"))
  expect_equal(20, fixStruct.lookup("EEI"))
  expect_equal(20, fixStruct.lookup("VVI"))
  expect_equal(30, fixStruct.lookup("VVV"))
  expect_equal(30, fixStruct.lookup("EEE"))
})


context("test skewness")
test_that("test skewness", {
  set.seed(123)
  alpha.seq <- c(0.1, 1, 2, 5)
  for (a in seq_along(alpha.seq)) {
    expect_equal(2 / sqrt(alpha.seq[a]),
      mean(replicate(1000, skewness(rgamma(10000, alpha.seq[a], 1)))),
      tolerance = 0.01
    )
  }
})
