library(clustTMB)
library(MoEClust)
library(testthat)
library(MixSim)

# load data
data("CO2data")

CO2 <- CO2data$CO2
GNP <- CO2data$GNP

# Set-up MoEClust models
m1 <- MoE_clust(CO2, G = 1:3)
m2 <- MoE_clust(CO2, G = 2:3, gating = ~GNP, network.data = CO2data)
m3 <- MoE_clust(CO2, G = 1:3, expert = ~GNP, network.data = CO2data)
m4 <- MoE_clust(CO2, G = 2:3, gating = ~GNP, expert = ~GNP, network.data = CO2data)

# match MoEclust alues
cT.m1 <- clustTMB(
  response = CO2,
  G = m1$G, covariance.structure = m1$modelName
)
cT.m2 <- clustTMB(
  response = CO2,
  gatingformula = ~GNP,
  gatingdata = CO2data,
  G = m2$G, covariance.structure = m2$modelName
)
cT.m3 <- clustTMB(
  response = CO2,
  expertformula = ~GNP,
  expertdata = CO2data,
  G = m3$G, covariance.structure = m3$modelName
)
cT.m4 <- clustTMB(
  response = CO2,
  gatingformula = ~GNP, expertformula = ~GNP,
  expertdata = CO2data,
  gatingdata = CO2data,
  G = m4$G, covariance.structure = m4$modelName
)



test_that("covariate test, m1", {
  expect_equal(-summary(m1)$loglik,
    cT.m1$opt$objective,
    tolerance = .01
  )
  expect_equal(1, ClassProp(
    as.vector(m1$classification),
    cT.m1$report$classification + 1
  ))
})

test_that("covariate test, m2", {
  expect_equal(-summary(m2)$loglik,
    cT.m2$opt$objective,
    tolerance = .01
  )
  expect_equal(1, ClassProp(
    as.vector(m2$classification),
    cT.m2$report$classification + 1
  ))
})

test_that("covariate test, m3", {
  expect_equal(-summary(m3)$loglik,
    cT.m3$opt$objective,
    tolerance = .01
  )
  expect_equal(1, ClassProp(
    as.vector(m3$classification),
    cT.m3$report$classification + 1
  ))
})

# tolerance higher
test_that("covariate test, m4", {
  expect_equal(-summary(m4)$loglik,
    cT.m4$opt$objective,
    tolerance = .1
  )
  expect_equal(1, ClassProp(
    as.vector(m4$classification),
    cT.m4$report$classification + 1
  ))
})
