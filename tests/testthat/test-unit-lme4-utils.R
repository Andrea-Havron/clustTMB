stopifnot(
  require("testthat"),
  require("clustTMB")
)

context("test lme4 utils functions")

test_that("test findReTrmClasses", {
  expect_equal(
    c("na", "norm", "ar1", "gmrf", "gmrf_speedup"),
    clustTMB:::findReTrmClasses()
  )
})

t <- 1:100
formula <- ~ ar1(0 + 1 | t)
specials <- clustTMB:::findReTrmClasses()


test_that("test splitForm", {
  split.glmm <- glmmTMB:::splitForm(formula)
  split.clust <- clustTMB:::splitForm(formula)
  expect_equal(split.glmm, split.clust)
})
