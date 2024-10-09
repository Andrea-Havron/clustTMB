require(TMB)
run_model <- FALSE
compile(testthat::test_path("ref", "simple.cpp"))
if (exists(dynlib(testthat::test_path("ref", "simple")))) {
  dyn.load(dynlib(testthat::test_path("ref", "simple")))
  run_model <- TRUE
}
## Test data
set.seed(123)
y <- rep(1900:2010, each = 2)
year <- factor(y)
quarter <- factor(rep(1:4, length.out = length(year)))
period <- factor((y > mean(y)) + 1)
## Random year+quarter effect, fixed period effect:
B <- model.matrix(~ year + quarter - 1)
A <- model.matrix(~ period - 1)
B <- as(as(as(B, "dMatrix"), "generalMatrix"), "TsparseMatrix")
A <- as(as(as(A, "dMatrix"), "generalMatrix"), "TsparseMatrix")
u <- rnorm(ncol(B)) ## logsdu=0
beta <- rnorm(ncol(A)) * 100
eps <- rnorm(nrow(B), sd = 1) ## logsd0=0
x <- as.numeric(A %*% beta + B %*% u + eps)

## Fit model
if (run_model) {
  obj <- MakeADFun(
    data = list(x = x, B = B, A = A),
    parameters = list(u = u * 0, beta = beta * 0, logsdu = 1, logsd0 = 1),
    random = "u",
    DLL = "simple",
    silent = TRUE
  )
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  sdr <- sdreport(obj)
  report <- obj$report(obj$env$last.par.best)


  args <- list(
    data = list(x = x, B = B, A = A),
    parameters = list(u = u * 0, beta = beta * 0, logsdu = 1, logsd0 = 1),
    random = "u",
    DLL = "simple",
    silent = TRUE
  )
}
