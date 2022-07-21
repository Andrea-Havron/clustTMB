stopifnot(
  require("testthat"),
  require("clustTMB")
)

context("test init.method input")
init.method.names <- c(
  "hc", "quantile", "random",
  "mclust", "kmeans", "mixed"
)
test_that("test default", {
  expect_equal("hc", init.options()$init.method)
})

test_that("test correct user init.method input", {
  for (nm in seq_along(init.method.names)) {
    expect_equal(
      init.method.names[nm],
      init.options(init.method = init.method.names[nm])$init.method
    )
  }
})

test_that("test incorrect user init.method input", {
  expect_error(init.options(init.method = "a"))
  expect_error(init.options(init.method = 1))
  expect_error(init.options(init.method = c("hc", "random")))
})


context("test hc.options input")
hc.options.names <- list(
  modelName = c("VVV", "EII", "EEE", "VII", "V", "E"),
  use = c("SVD", "VARS", "STD", "SPH", "PCS", "PCR", "RND")
)
test_that("test default", {
  expect_equal("VVV", init.options()$hc.options$modelName)
  expect_equal("SVD", init.options()$hc.options$use)
  expect_equal(
    "VVV",
    init.options(
      hc.options = list(use = "SVD")
    )$hc.options$modelName
  )
  expect_equal(
    "SVD",
    init.options(
      hc.options = list(modelName = "VVV")
    )$hc.options$use
  )
})
test_that("test correct user input", {
  for (nm in seq_along(hc.options.names$modelName)) {
    expect_equal(
      hc.options.names$modelName[nm],
      init.options(
        hc.options =
          list(modelName = hc.options.names$modelName[nm])
      )$hc.options$modelName
    )
  }
  for (nm in seq_along(hc.options.names$use)) {
    expect_equal(
      hc.options.names$use[nm],
      init.options(
        hc.options =
          list(use = hc.options.names$use[nm])
      )$hc.options$use
    )
  }

  expect_equal(
    "VVV",
    init.options(
      hc.options = list("VVV", "SVD")
    )$hc.options$modelName
  )
  expect_equal(
    "VVV",
    init.options(
      hc.options = list("VVV")
    )$hc.options$modelName
  )
  expect_equal(
    "VVV",
    init.options(
      hc.options = list("SVD")
    )$hc.options$modelName
  )

  expect_equal(
    "SVD",
    init.options(
      hc.options = list("VVV", "SVD")
    )$hc.options$use
  )
  expect_equal(
    "SVD",
    init.options(
      hc.options = list("VVV")
    )$hc.options$use
  )
  expect_equal(
    "SVD",
    init.options(
      hc.options = list("SVD")
    )$hc.options$use
  )
})
test_that("test incorrect user input", {
  expect_error(init.options(hc.options = c("VVV", "SVD")))
  expect_error(init.options(hc.options = list(modelName = "a")))
  expect_error(init.options(hc.options = list(use = "a")))
  expect_error(init.options(hc.options = "VVV"))
  expect_error(init.options(hc.options = list(modelName = c("VVV", "EEE"))))
  expect_error(init.options(hc.options = list("VVV", "EEE", "SVD")))
})

# TODO: exp.init tests

context("test mix.method input")
test_that("test default input", {
  expect_equal("Gower kmeans", init.options()$mix.method)
})
test_that("test correct input", {
  mix.method.nms <- c("Gower kmeans", "Gower hclust", "kproto")
  for (nm in seq_along(mix.method.nms)) {
    expect_equal(
      mix.method.nms[nm],
      init.options(mix.method = mix.method.nms[nm])$mix.method
    )
  }
})
test_that("test incorrect input", {
  expect_error(init.options(mix.method = "a"))
  expect_error(init.options(mix.method = 1))
  expect_error(init.options(mix.method = c("Gower kmeans", "kproto")))
})

context("test user class")
test_that("test default input", {
  expect_equal(integer(0), init.options()$user.class)
})
test_that("test correct input", {
  set.seed(123)
  input.class <- rbinom(10, 3, .5)
  expect_equal(input.class + 1, init.options(
    init.method = "user",
    user.class = input.class
  )$user.class)
  expect_equal(
    as.integer(as.factor(c("a", "b", "c", "d"))),
    init.options(user.class = c("a", "b", "c", "d"))$user.class
  )
  expect_equal(
    as.integer(as.factor(c("a", "b", "c", "d"))),
    init.options(user.class = as.factor(c("a", "b", "c", "d")))$user.class
  )
  input.class <- c(2.1, 2.1, 2.1, 4.0, 4.0, 4.0, 5.3, 5.3)
  expect_equal(
    as.numeric(factor(input.class)),
    init.options(user.class = input.class)$user.class
  )
  # depends on implementation: 
  # expect_message(init.options(user.class = input.class)$user.class)
})
test_that("test incorrect input", {
  expect_error(init.options(init.method = "user"))
})

context("test default")
test_that("test default", {
  expect_equal(
    c("init.method", "hcName", "hcUse", "mix.method"),
    init.options()$defaults
  )
})
