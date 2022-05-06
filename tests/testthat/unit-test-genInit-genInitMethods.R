stopifnot(
  require("testthat"),
  require("clustTMB"),
  require("mclust")
)

# Methods based on clustTMB::init.options()
# init.methods: c("hc", "quantile", "random", 
#                 "mclust", "kmeans", "mixed", "user")

context("test init.methods")
data("faithful")
ng <- 2
test_that("test hc using default settings", {
 
  #clustTMB default use = "SVD"
  clss.true <- as.vector(
    hclass(hc(faithful, use = "SVD"), ng))
  
  clss.clustTMB <- genInitMethods(n.g = ng, n.i = nrow(faithful),
                                  n.j = ncol(faithful),
                                  control = init.options(),
                                  y = faithful)
  expect_equal(clss.true, clss.clustTMB)
})

test_that("test hc using user defined settings", {
  
  clss.true <- as.vector(
    hclass(hc(faithful, modelName = "EEE", use = "VARS"), ng)
  )
  
  clss.clustTMB <- genInitMethods(n.g = ng, n.i = nrow(faithful),
                                  n.j = ncol(faithful),
                                  control = init.options(
                                    hc.options = list(modelName = "EEE",
                                                      use = "VARS")
                                  ),
                                  y = faithful)
  expect_equal(clss.true, clss.clustTMB)
})

test_that("test quantile", {
  
  #univariate method
  clss.true <- clustTMB:::mc.qclass(faithful[,1], as.numeric(ng))
  
  clss.clustTMB <- genInitMethods(n.g = ng, n.i = nrow(faithful),
                                  n.j = ncol(faithful),
                                  control = init.options(
                                    init.method = "quantile"
                                  ),
                                  y = faithful[,1])
  expect_equal(clss.true, clss.clustTMB)
})

test_that("test random", {
  
  #random method, so set seed
  set.seed(123)
  clss.true <- sample(1:ng, nrow(faithful), replace = TRUE)
  
  set.seed(123)
  clss.clustTMB <- genInitMethods(n.g = ng, n.i = nrow(faithful),
                                  n.j = ncol(faithful),
                                  control = init.options(
                                    init.method = "random"
                                  ),
                                  y = faithful)
  expect_equal(clss.true, clss.clustTMB)
})

test_that("test kmeans", {
  
  clss.true <- cluster::pam(faithful, k = ng)$clustering
  
  clss.clustTMB <- genInitMethods(n.g = ng, n.i = nrow(faithful),
                                  n.j = ncol(faithful),
                                  control = init.options(
                                    init.method = "kmeans"
                                  ),
                                  y = faithful)
  expect_equal(clss.true, clss.clustTMB)
})

test_that("test kmeans", {
  
  clss.true <- cluster::pam(faithful, k = ng)$clustering
  
  clss.clustTMB <- genInitMethods(n.g = ng, n.i = nrow(faithful),
                                  n.j = ncol(faithful),
                                  control = init.options(
                                    init.method = "user",
                                    user.class = clss.true
                                  ),
                                  y = faithful)
  expect_equal(unname(clss.true), clss.clustTMB)
})

## TODO: test mixed 
