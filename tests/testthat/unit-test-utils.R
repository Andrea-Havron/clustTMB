stopifnot(
  require("testthat"),
  require("clustTMB")
)

context("test mkFac")
test_that("test mkFac", {
  d <- c(3,3)
  map.vec <- rep(factor(NA), d[1]*d[2])
  map.dim <- map.vec
  dim(map.dim) <- d
  expect_equal(map.dim, mkFac(d, map.vec))
  
  d <- 3
  map.vec <- rep(factor(NA), d)
  expect_equal(map.vec, mkFac(d, map.vec))
})

context("test families")
test_that("tweedie",{
  f <- tweedie()
  expect_equal("tweedie", f$family)
  expect_equal("log", f$link)
})

test_that("lognormal",{
  f <- lognormal()
  expect_equal("lognormal", f$family)
  expect_equal("identity", f$link)
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

context("test parm.lookup")
test_that("test parm.lookup", {
  expect_equal(23, nrow(parm.lookup()$parm))
  expect_equal(c("Fixed", "Random"), unique(parm.lookup()$parm$type))
  expect_equal(c("Matrix", "Array", "Vector"), unique(parm.lookup()$parm$str))
  expect_equal(10, nrow(parm.lookup()$key))
  expect_equal(c("Fr", "Fs", "G", "J", "Kd", "Kg", "M", "N", "T", "V"), parm.lookup()$key$dim)
})
