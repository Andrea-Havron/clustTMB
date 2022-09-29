library(clustTMB)
library(mclust)
library(sp)
library(testthat)

#load data
data("meuse")

#not working
# mclust.mod <- Mclust(meuse[,3:6], modelNames = "VVV")
# #match mclust default values
# clustTMB.mod <- clustTMB(response = as.matrix(meuse[,3:6]), 
#          G=mclust.mod$G, covariance.structure = mclust.mod$modelName)
# 
# test_that("simple integrated test", {
#   expect_equal(-mclust.mod$loglik, 
#                clustTMB.mod$opt$objective, tolerance = .01)
#   expect_equal(1,  MixSim::ClassProp(mclust.mod$classification,  
#                                      clustTMB.mod$report$classification+1))
# })

data("faithful")
mclust.mod <- Mclust(faithful, modelNames = "VVV")
clustTMB.mod <- clustTMB(response = as.matrix(faithful), 
                         G=mclust.mod$G, covariance.structure = mclust.mod$modelName)

test_that("simple integrated test", {
  expect_equal(-mclust.mod$loglik, 
               clustTMB.mod$opt$objective, tolerance = .0001)
  expect_equal(1,  MixSim::ClassProp(mclust.mod$classification,  
                                     clustTMB.mod$report$classification+1))
})
