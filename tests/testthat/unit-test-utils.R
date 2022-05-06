stopifnot(
  require("testthat"),
  require("clustTMB")
)

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
  for(a in seq_along(alpha.seq)){
    expect_equal(2 / sqrt(alpha.seq[a]), 
                 mean(replicate(1000, skewness(rgamma(10000,alpha.seq[a],1)))), 
                 tolerance = 0.01)
  }
})