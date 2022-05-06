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
