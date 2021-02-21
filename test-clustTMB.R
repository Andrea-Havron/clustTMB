#test clustTMB covariate model
stopifnot(require("testthat"),
          require("clustTMB"),
          require("MoEClust"))

data(CO2data)
CO2 <- CO2data$CO2
GNP <- CO2data$GNP


m5    <- MoE_clust(CO2, G=2:3, equalPro=TRUE, verbose=FALSE) ##! need to set up option to map out betag for equalPro model
m6    <- MoE_clust(CO2, G=2:3, expert= ~ GNP, equalPro=TRUE, verbose=FALSE)
comp  <- MoE_compare(m1, m2, m3, m4, m5, m6, optimal.only=TRUE)
summary(comp)

Context('Compare clustTMB to MoEClust')
test_that('no covariates', {
  m1    <- MoE_clust(CO2, G=1:3, verbose=FALSE)
  m1.cT <- clustTMB(CO2, G = m1$G, covariance.structure = 'E')
  expect_equal(as.vector(unlist(summary(m1)[7])), -m1.cT$opt$objective, tolerance = 0.01)
})
test_that('gating covariate', {
  m2    <- MoE_clust(CO2, G=2:3, gating= ~ GNP, verbose=FALSE)
  m2.cT <- clustTMB(CO2, gatingformula = ~GNP, gatingdata = CO2data,  G = m2$G, covariance.structure = 'E')
  expect_equal(as.vector(unlist(summary(m2)[7])), -m2.cT$opt$objective, tolerance = 0.1)
})
test_that('expert covariate', {
  m3    <- MoE_clust(CO2, G=1:3, expert= ~ GNP, verbose=FALSE)
  m3.cT <- clustTMB(CO2, expertformula = ~GNP, expertdata = CO2data,  G = m3$G, covariance.structure = summary(m3)$modelName)
  expect_equal(as.vector(unlist(summary(m3)[7])), -m3.cT$opt$objective, tolerance = 0.1)
})
testthat('expert and gating covariate', {
  m4    <- MoE_clust(CO2, G=2:3, gating= ~ GNP, expert= ~ GNP, verbose=FALSE)
  m4.cT <- clustTMB(CO2, gatingformula = ~GNP, expertformula = ~GNP, gatingdata = CO2data,
                    expertdata = CO2data,  G = m4$G,
                    covariance.structure = summary(m4)$modelName)
  expect_equal(as.vector(unlist(summary(m4)[7])), -m4.cT$opt$objective, tolerance = 0.1)
})


