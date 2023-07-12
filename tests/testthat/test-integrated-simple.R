library(clustTMB)
library(mclust)
library(sp)
library(testthat)
library(magrittr)
library(tidyr)
library(dplyr)
library(GGally)
#load data
library(palmerpenguins)

context("integrated simple test")

pen <- penguins_raw %>%
  dplyr::mutate(id = row_number(),
         `Culmen Length (mm)` = scale(`Culmen Length (mm)`),
         `Culmen Depth (mm)` = scale(`Culmen Depth (mm)`),
         `Flipper Length (mm)` = scale(`Flipper Length (mm)`),
         `Body Mass (g)` = scale(`Body Mass (g)`)) %>%
  dplyr::select(id, Species, Region, Island, `Culmen Length (mm)`,
         `Culmen Depth (mm)`, `Flipper Length (mm)`, 
         `Body Mass (g)`) %>%
  tidyr::drop_na(`Culmen Length (mm)`) 

dat <- pen %>% dplyr::select(`Culmen Length (mm)`, `Culmen Depth (mm)`, 
  `Flipper Length (mm)`, `Body Mass (g)`) %>%
  as.data.frame() %>% as.matrix()

mclust.mod <- Mclust(dat, G = 3, modelNames = "VII")

clustTMB.mod <- clustTMB(response = dat, G=3, 
                         covariance.structure = "VII")


test_that("simple integrated test", {
  expect_equal(-mclust.mod$loglik,
               clustTMB.mod$opt$objective, tolerance = .01)
  expect_equal(1,  MixSim::ClassProp(mclust.mod$classification,
                                     clustTMB.mod$report$classification),
               tolerance = 0.01)
  expect_equal(1,  MixSim::ClassProp(mclust.mod$classification,  
                                     pen$Species), tolerance = 0.1)
  expect_equal(1,  MixSim::ClassProp(clustTMB.mod$report$classification,  
                                     pen$Species), tolerance = 0.1)
})

data("faithful")
mclust.mod <- Mclust(faithful, G = 2, modelNames = "VVV")
clustTMB.mod <- clustTMB(response = as.matrix(faithful), 
                         G=2, covariance.structure = "VVV")

test_that("simple integrated test", {
  expect_equal(-mclust.mod$loglik, 
               clustTMB.mod$opt$objective, tolerance = .0001)
  expect_equal(1,  MixSim::ClassProp(mclust.mod$classification,  
                                     clustTMB.mod$report$classification))
})

#Add more tests for different clustTMB::fixStruct.names() using above two data sets 
#can use first column of faithful to test "E" and "V"
