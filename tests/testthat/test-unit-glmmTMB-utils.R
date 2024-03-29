stopifnot(
  require("testthat"),
  require("clustTMB"),
  require("glmmTMB")
)

context("test glmmTMB utils functions")

test_that("test findReTrmClasses", {
  expect_equal(
    c("na", "norm", "ar1", "gmrf", "gmrf_speedup"),
    clustTMB:::findReTrmClasses()
  )
})

t <- 1:100
formula <- ~ ar1(0 + 1 | t)
specials <- clustTMB:::findReTrmClasses()

test_that("test findbars_x, default", {
  findbars.glmm <- glmmTMB::findbars_x(formula, FALSE, specials)
  findbars.clust <- clustTMB::findbars_x(formula, FALSE, specials)
  expect_equal(findbars.glmm, findbars.clust)
  expect_equal(quote(ar1(0 + 1 | t)), findbars.clust[[1]])
})

glmmTMB:::findbars_x(~ 1 + (1 | f) + (1 | g))
clustTMB:::findbars_x(~ 1 + (1 | f) + (1 | g))
test_that("test splitForm", {
  split.glmm <- glmmTMB:::splitForm(formula)
  split.clust <- clustTMB:::splitForm(formula)
  expect_equal(split.glmm, split.clust)
})
