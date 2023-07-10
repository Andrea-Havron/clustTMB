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

ggpairs(
  data = penguins_raw, 
  column = c(10:13), 
  mapping = aes(color = Species),
  diag = list(continuous = wrap("barDiag", color = "blue", size =4)),
  upper = list(continuous = wrap("cor", size = 4, bins = 60))
)

pen <- penguins_raw %>%
  dplyr::mutate(is = row_number(),
         `Culmen Length (mm)` = scale(`Culmen Length (mm)`),
         `Culmen Depth (mm)` = scale(`Culmen Depth (mm)`),
         `Flipper Length (mm)` = scale(`Flipper Length (mm)`),
         `Body Mass (g)` = scale(`Body Mass (g)`)) %>%
  dplyr::select(id, Species, Region, Island, `Culmen Length (mm)`,
         `Culmen Depth (mm)`, `Flipper Length (mm)`, 
         `Body Mass (g)`) %>%
  tidyr::drop_na(`Culmen Length (mm)`) 
ggpairs(
    data = pen, 
    columns = c(5:8), 
    mapping = aes(color = mclust.mod$classification),
    diag = list(continuous = wrap("barDiag", color = "blue", size =4)),
    upper = list(continuous = wrap("cor", size = 4, bins = 60))
  )

mclust.mod <- Mclust(pen[,5:8], modelNames = "VVV")


clustTMB.mod <- clustTMB(response = as.matrix(pen[,5:8]),
                         G=mclust.mod$G, covariance.structure = mclust.mod$modelName)


test_that("simple integrated test", {
  expect_equal(-mclust.mod$loglik,
               clustTMB.mod$opt$objective, tolerance = .01)
  expect_equal(1,  MixSim::ClassProp(mclust.mod$classification,
                                     clustTMB.mod$report$classification+1))
  expect_equal(1,  MixSim::ClassProp(mclust.mod$classification,  
                                     pen$Species), tolerance = 0.1)
  expect_equal(1,  MixSim::ClassProp(mclust.mod$classification,  
                                     pen$Species), tolerance = 0.1)
})

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
