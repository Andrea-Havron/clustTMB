stopifnot(require("testthat"),
          require("clustTMB"),
          require("mclust"),
          require('mvnfast'))

# multivariate normal data
data('faithful')
clss <- unmap(hclass(hc(faithful), 2))
init.parm.mclust <- list(mu = list(g1 = apply(faithful[clss[,1]==1,],2, mean),
                               g2 = apply(faithful[clss[,2]==1,],2,mean)),
                     Sigma = list(g1 = apply(faithful[clss[,1]==1,],2,var),
                                  g2 = apply(faithful[clss[,2]==1,],2,var)),
                     Corr = list(g1 = cor(faithful[clss[,1]==1,])[2,1],
                                 g2 = cor(faithful[clss[,2]==1,])[2,1]),
                     Pi = list(g1 = sum(clss[,1])/nrow(faithful),
                               g2 = sum(clss[,2])/nrow(faithful)))
dim.list <- list(n.i = nrow(faithful), n.j = ncol(faithful), n.t = 1,
                 n.g = 2, n.f.sp = ncol(faithful), n.f.rand = ncol(faithful), n.v = NULL)
Dat <- mkDat(response = as.matrix(faithful), time.vector = rep(1, nrow(faithful)),
             expert.dat = as.matrix(rep(1, nrow(faithful))),
             gating.dat = as.matrix(rep(1, nrow(faithful))),
             ll.method = 0,
             fixStruct = 'VVV',
             rrStruct = rep(0,2),
             reStruct = matrix(0,2,3),
             dim.list =  dim.list
)
dim.list$n.v <- Dat$spde$n_s
init.parm.clustTMB <- genInit('hc', Dat, family = gaussian(link = "identity"), dim.list, data.trans = NA,
                              hc.options = c("VVV", "SVD"), #default matches mclust
                              mix.method = NULL, r.tol = 1e-10, true.class = NULL)

context("equal mvn init parameters")
test_that("mu", {
  expect_equal(as.vector(init.parm.mclust$mu$g1), init.parm.clustTMB$parms$betad[1,,1])
  expect_equal(as.vector(init.parm.mclust$mu$g2), init.parm.clustTMB$parms$betad[1,,2])
})
test_that("Sigma", {
  expect_equal(as.vector(init.parm.mclust$Sigma$g1), exp(init.parm.clustTMB$parms$theta[,1]))
  expect_equal(as.vector(init.parm.mclust$Sigma$g2), exp(init.parm.clustTMB$parms$theta[,2]))
})
test_that("Corr", {
  expect_equal(as.vector(unlist(init.parm.mclust$Corr)),
               as.vector(1/(1+exp(-init.parm.clustTMB$parms$logit_corr_fix))*2-1))
})
test_that("Pi", {
  expect_equal(as.vector(unlist(init.parm.mclust$Pi)),
               c(1/(1+exp(-init.parm.clustTMB$parms$betag)),
                 1-1/(1+exp(-init.parm.clustTMB$parms$betag))))
})
test_that("Class", {
  expect_equal(map(clss), init.parm.clustTMB$class)
})
test_that("Random Effects dim", {
  expect_equal(c(dim.list$n.t, dim.list$n.g-1), dim(init.parm.clustTMB$parms$upsilon_tg))
  expect_equal(c(dim.list$n.t, dim.list$n.j, dim.list$n.g), dim(init.parm.clustTMB$parms$epsilon_tjg))
  expect_equal(c(dim.list$n.i, dim.list$n.g-1), dim(init.parm.clustTMB$parms$u_ig))
  expect_equal(c(dim.list$n.i, dim.list$n.f.rand, dim.list$n.g), dim(init.parm.clustTMB$parms$v_ifg))
  expect_equal(c(dim.list$n.v, dim.list$n.g-1), dim(init.parm.clustTMB$parms$Gamma_vg))
  expect_equal(c(dim.list$n.v, dim.list$n.f.sp, dim.list$n.g), dim(init.parm.clustTMB$parms$Omega_vfg))
})


clss <- unmap(mclust:::qclass(faithful$waiting, 2))
init.parm.mclust <- list(mu = list(g1 = mean(faithful$waiting[clss[,1]==1]),
                                   g2 = mean(faithful$waiting[clss[,2]==1])),
                         Sigma = list(g1 = var(faithful$waiting[clss[,1]==1]),
                                      g2 = var(faithful$waiting[clss[,2]==1])),
                         Corr = list(g1 = 0, g2 = 0),
                         Pi = list(g1 = sum(clss[,1])/nrow(faithful),
                                   g2 = sum(clss[,2])/nrow(faithful)))
dim.list <- list(n.i = nrow(faithful), n.j = 1, n.t = 1,
                 n.g = 2, n.f.sp = 1, n.f.rand = 1, n.v = NULL)
Dat <- mkDat(response = as.matrix(faithful$waiting), time.vector = rep(1, nrow(faithful)),
             expert.dat = as.matrix(rep(1, nrow(faithful))),
             gating.dat = as.matrix(rep(1, nrow(faithful))),
             ll.method = 0,
             fixStruct = 'E',
             rrStruct = rep(0,2),
             reStruct = matrix(0,2,3),
             dim.list =  dim.list
)
dim.list$n.v <- Dat$spde$n_s
init.parm.clustTMB <- genInit(init.method = 'mc.qclass', Dat, family = gaussian(link = "identity"),
                              dim.list, data.trans = NA, hc.options = c("VVV", "SVD"), #default matches mclust
                              mix.method = NULL, r.tol = 1e-10, true.class = NULL)

context("equal uniNormal init parameters")
test_that("mu", {
  expect_equal(as.vector(init.parm.mclust$mu$g1), init.parm.clustTMB$parms$betad[1,,1])
  expect_equal(as.vector(init.parm.mclust$mu$g2), init.parm.clustTMB$parms$betad[1,,2])
})
test_that("Sigma", {
  expect_equal(as.vector(init.parm.mclust$Sigma$g1), exp(init.parm.clustTMB$parms$theta[,1]))
  expect_equal(as.vector(init.parm.mclust$Sigma$g2), exp(init.parm.clustTMB$parms$theta[,2]))
})
test_that("Corr", {
  expect_equal(as.vector(unlist(init.parm.mclust$Corr)),
               as.vector(1/(1+exp(-init.parm.clustTMB$parms$logit_corr_fix))*2-1))
})
test_that("Pi", {
  expect_equal(as.vector(unlist(init.parm.mclust$Pi)),
               c(1/(1+exp(-init.parm.clustTMB$parms$betag)),
                 1-1/(1+exp(-init.parm.clustTMB$parms$betag))))
})
test_that("Class", {
  expect_equal(map(clss), init.parm.clustTMB$class)
})


#test MakeADFun
test_that("Random Effects dim", {
  expect_equal(c(dim.list$n.t, dim.list$n.g-1), dim(init.parm.clustTMB$parms$upsilon_tg))
  expect_equal(c(dim.list$n.t, dim.list$n.j, dim.list$n.g), dim(init.parm.clustTMB$parms$epsilon_tjg))
  expect_equal(c(dim.list$n.i, dim.list$n.g-1), dim(init.parm.clustTMB$parms$u_ig))
  expect_equal(c(dim.list$n.i, dim.list$n.f.rand, dim.list$n.g), dim(init.parm.clustTMB$parms$v_ifg))
  expect_equal(c(dim.list$n.v, dim.list$n.g-1), dim(init.parm.clustTMB$parms$Gamma_vg))
  expect_equal(c(dim.list$n.v, dim.list$n.f.sp, dim.list$n.g), dim(init.parm.clustTMB$parms$Omega_vfg))
})

# context("univariate normal data with covaraites")
# data(CO2data)
# CO2 <- CO2data$CO2
# GNP <- CO2data$GNP
#
# clss <- unmap(mclust:::qclass(CO2, 2))
#
# dim.list <- list(n.i = length(CO2), n.j = 1, n.t = 1,
#                  n.g = 2, n.f.sp = 1, n.f.rand = 1, n.v = NULL)
# Dat <- mkDat(response = as.matrix(CO2), time.vector = rep(1, length(CO2)),
#              expert.dat = as.matrix(cbind(rep(1,length(CO2)), GNP)),
#              gating.dat = as.matrix(cbind(rep(1,length(CO2)), GNP)),
#              ll.method = 1,
#              fixStruct = 'E',
#              rrStruct = rep(0,2),
#              reStruct = matrix(0,2,3),
#              dim.list =  dim.list
# )
# dim.list$n.v <- Dat$spde$n_s
# init.parm.clustTMB <- genInit(init.method = NULL, Dat, dim.list, data.trans = NA, hc.options = c("VVV", "SVD"), #default matches mclust
#                               mix.method = NULL, r.tol = 1e-10, true.class = NULL)
# test_that("mu", {
#   expect_equal(as.vector(init.parm.mclust$mu$g1), init.parm.clustTMB$parms$betad[1,,1])
#   expect_equal(as.vector(init.parm.mclust$mu$g2), init.parm.clustTMB$parms$betad[1,,2])
# })
# test_that("Sigma", {
#   expect_equal(as.vector(init.parm.mclust$Sigma$g1), exp(init.parm.clustTMB$parms$theta[,1]))
#   expect_equal(as.vector(init.parm.mclust$Sigma$g2), exp(init.parm.clustTMB$parms$theta[,2]))
# })
# test_that("Corr", {
#   expect_equal(as.vector(unlist(init.parm.mclust$Corr)),
#                as.vector(1/(1+exp(-init.parm.clustTMB$parms$logit_corr_fix))*2-1))
# })
# test_that("Pi", {
#   expect_equal(as.vector(unlist(init.parm.mclust$Pi)),
#                c(1/(1+exp(-init.parm.clustTMB$parms$betag)),
#                  1-1/(1+exp(-init.parm.clustTMB$parms$betag))))
# })
# test_that("Class", {
#   expect_equal(map(clss), init.parm.clustTMB$class)
# })
