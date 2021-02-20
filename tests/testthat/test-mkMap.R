stopifnot(require("testthat"),
          require("clustTMB"),
          require('mvnfast'))

context('match mkMap with init parm')
n.i <- 100
n.j <- 4
n.g <- 3
y <- rmvn(n.i, rep(0, n.j), diag(n.j))

dim.list <- list(n.i = n.i, n.j = n.j, n.t = 1,
                 n.g = n.g, n.f.sp = n.j, n.f.rand = n.j, n.v = NULL)


test_that('mvn with no random effect, no rank reduction', {
  Dat <- mkDat(response = as.matrix(y), time.vector = rep(1, dim.list$n.i),
               expert.dat = as.matrix(rep(1, dim.list$n.i)),
               gating.dat = as.matrix(rep(1, dim.list$n.i)),
               ll.method = 0,
               fixStruct = 'VVV',
               rrStruct = rep(0,2),
               reStruct = matrix(0,2,3),
               dim.list =  dim.list
  )
  dim.list$n.v <- Dat$spde$n_s
  init.parm <- genInit('hc', Dat, family = gaussian(link = "identity"),
                       dim.list, data.trans = NA, hc.options = c("VVV", "SVD"), #default matches mclust
                       mix.method = NULL, r.tol = 1e-10, true.class = NULL)
  map.list <- mkMap(Dat$family, Dat$fixStruct, Dat$rrStruct, Dat$reStruct, dim.list)
  map.names <- names(map.list)
  exp.map.names <- c("thetaf", "ld_rand", "ld_sp", "Hg_input", "Hd_input",
                     "ln_kappag", "ln_kappad", "ln_taud", "logit_rhog", "logit_rhod",
                     "ln_sigmaup",  "ln_sigmaep", "ln_sigmau", "ln_sigmav", "upsilon_tg",
                     "epsilon_tjg", "u_ig", "v_ifg", "Gamma_vg", "Omega_vfg"   )
  for(m in 1:length(map.names)){
    expect_equal(sapply(map.list[map.names[m]], dim), sapply(init.parm$parms[map.names[m]], dim))
  }
  expect_equal(sort(map.names), sort(exp.map.names))
})

test_that('mvn with gating random effects, no rank reduction', {
  Dat <- mkDat(response = as.matrix(y), time.vector = rep(1, dim.list$n.i),
               expert.dat = as.matrix(rep(1, dim.list$n.i)),
               gating.dat = as.matrix(rep(1, dim.list$n.i)),
               ll.method = 0,
               fixStruct = 'VVV',
               rrStruct = rep(0,2),
               reStruct = matrix(0,2,3),
               dim.list =  dim.list
  )
  dim.list$n.v <- Dat$spde$n_s
  exp.map.names <- c("thetaf", "ld_rand", "ld_sp", "Hg_input", "Hd_input",
                     "ln_kappag", "ln_kappad", "ln_taud", "logit_rhog", "logit_rhod",
                     "ln_sigmaup",  "ln_sigmaep", "ln_sigmau", "ln_sigmav", "upsilon_tg",
                     "epsilon_tjg", "u_ig", "v_ifg", "Gamma_vg", "Omega_vfg"   )
  reNum <- c(3, 2, 1)
  parmName <- list("Gamma_vg", c("logit_rhog", "ln_sigmaup", "upsilon_tg"), c("ln_sigmau","u_ig"))
  for(j in 1:3){
    Dat$reStruct <- matrix(0,2,3)
    Dat$reStruct[1,j] <- reNum[j]
    init.parm <- genInit('hc', Dat, family = gaussian(link = "identity"),
                         dim.list, data.trans = NA, hc.options = c("VVV", "SVD"), #default matches mclust
                         mix.method = NULL, r.tol = 1e-10, true.class = NULL)
    map.list <- mkMap(Dat$family, Dat$fixStruct, Dat$rrStruct, Dat$reStruct, dim.list)
    map.names <- names(map.list)
    exp.map.names <- c("thetaf", "ld_rand", "ld_sp", "Hg_input", "Hd_input",
                       "ln_kappag", "ln_kappad", "ln_taud", "logit_rhog", "logit_rhod",
                       "ln_sigmaup",  "ln_sigmaep", "ln_sigmau", "ln_sigmav", "upsilon_tg",
                       "epsilon_tjg", "u_ig", "v_ifg", "Gamma_vg", "Omega_vfg"   )
    for(i in 1:length(parmName[[j]])){
      exp.map.names <- exp.map.names[exp.map.names != parmName[[j]][i]]
    }

    for(m in 1:length(map.names)){
      expect_equal(sapply(map.list[map.names[m]], dim), sapply(init.parm$parms[map.names[m]], dim))
    }
    expect_equal(sort(map.names), sort(exp.map.names))
    if(j==1) expect_equal(sum(as.numeric(map.list$ln_kappag)), dim.list$n.g-1)
  }
})

test_that('mvn with expert random effects, no rank reduction', {
  Dat <- mkDat(response = as.matrix(y), time.vector = rep(1, dim.list$n.i),
               expert.dat = as.matrix(rep(1, dim.list$n.i)),
               gating.dat = as.matrix(rep(1, dim.list$n.i)),
               ll.method = 0,
               fixStruct = 'VVV',
               rrStruct = rep(0,2),
               reStruct = matrix(0,2,3),
               dim.list =  dim.list
  )
  dim.list$n.v <- Dat$spde$n_s
  reNum <- c(3, 2, 1)
  parmName <- list("Omega_vfg", c("logit_rhod", "ln_sigmaep", "epsilon_tjg"), c("ln_sigmav","v_ifg"))
  for(j in 1:3){
    Dat$reStruct <- matrix(0,2,3)
    Dat$reStruct[2,j] <- reNum[j]
    init.parm <- genInit('hc', Dat, family = gaussian(link = "identity"),
                         dim.list, data.trans = NA, hc.options = c("VVV", "SVD"), #default matches mclust
                         mix.method = NULL, r.tol = 1e-10, true.class = NULL)
    map.list <- mkMap(Dat$family, Dat$fixStruct, Dat$rrStruct, Dat$reStruct, dim.list)
    map.names <- names(map.list)
    exp.map.names <- c("thetaf", "ld_rand", "ld_sp", "Hg_input", "Hd_input",
                       "ln_kappag", "ln_kappad", "ln_taud", "logit_rhog", "logit_rhod",
                       "ln_sigmaup",  "ln_sigmaep", "ln_sigmau", "ln_sigmav", "upsilon_tg",
                       "epsilon_tjg", "u_ig", "v_ifg", "Gamma_vg", "Omega_vfg"   )
    for(i in 1:length(parmName[[j]])){
      exp.map.names <- exp.map.names[exp.map.names != parmName[[j]][i]]
    }

    for(m in 1:length(map.names)){
      expect_equal(sapply(map.list[map.names[m]], dim), sapply(init.parm$parms[map.names[m]], dim))
    }
    expect_equal(sort(map.names), sort(exp.map.names))
    if(j==1){
      expect_equal(sum(as.numeric(map.list$ln_kappad)), dim.list$n.g*dim.list$n.j)
      expect_equal(sum(as.numeric(map.list$ln_taud)), dim.list$n.g*dim.list$n.j)
    }
  }
})

test_that('mvn with expert random effects and rank reduction', {

  ## random reduction
  dim.list <- list(n.i = n.i, n.j = n.j, n.t = 1,
                   n.g = n.g, n.f.sp = n.j, n.f.rand = n.j-1, n.v = NULL)
  Dat <- mkDat(response = as.matrix(y), time.vector = rep(1, dim.list$n.i),
               expert.dat = as.matrix(rep(1, dim.list$n.i)),
               gating.dat = as.matrix(rep(1, dim.list$n.i)),
               ll.method = 0,
               fixStruct = 'VVV',
               rrStruct = rep(0,2),
               reStruct = matrix(0,2,3),
               dim.list =  dim.list
  )
  dim.list$n.v <- Dat$spde$n_s

  Dat$reStruct <- matrix(0,2,3)
  Dat$reStruct[2,3] <- 1
  Dat$rrStruct[1] <- 1
  init.parm <- genInit('hc', Dat, family = gaussian(link = "identity"),
                       dim.list, data.trans = NA, hc.options = c("VVV", "SVD"), #default matches mclust
                       mix.method = NULL, r.tol = 1e-10, true.class = NULL)
  map.list <- mkMap(Dat$family, Dat$fixStruct, Dat$rrStruct, Dat$reStruct, dim.list)
  map.names <- names(map.list)
  exp.map.names <- c("thetaf", "ld_sp", "Hg_input", "Hd_input",
                     "ln_kappag", "ln_kappad", "ln_taud", "logit_rhog", "logit_rhod",
                     "ln_sigmaup",  "ln_sigmaep", "ln_sigmau", "ln_sigmav", "upsilon_tg",
                     "epsilon_tjg", "u_ig",  "Gamma_vg", "Omega_vfg" )

  for(m in 1:length(map.names)){
    expect_equal(sapply(map.list[map.names[m]], dim), sapply(init.parm$parms[map.names[m]], dim))
  }
  expect_equal(sort(map.names), sort(exp.map.names))

  expect_equal(sum(is.na(map.list$ln_sigmaep)), dim.list$n.g*dim.list$n.j)

  ## spatial reduction
  dim.list <- list(n.i = n.i, n.j = n.j, n.t = 1,
                   n.g = n.g, n.f.sp = n.j-1, n.f.rand = n.j, n.v = NULL)
  Dat <- mkDat(response = as.matrix(y), time.vector = rep(1, dim.list$n.i),
               expert.dat = as.matrix(rep(1, dim.list$n.i)),
               gating.dat = as.matrix(rep(1, dim.list$n.i)),
               ll.method = 0,
               fixStruct = 'VVV',
               rrStruct = rep(0,2),
               reStruct = matrix(0,2,3),
               dim.list =  dim.list
  )
  dim.list$n.v <- Dat$spde$n_s
  Dat$reStruct <- matrix(0,2,3)
  Dat$reStruct[2,1] <- 3
  Dat$rrStruct[2] <- 1
  init.parm <- genInit('hc', Dat, family = gaussian(link = "identity"),
                       dim.list, data.trans = NA, hc.options = c("VVV", "SVD"), #default matches mclust
                       mix.method = NULL, r.tol = 1e-10, true.class = NULL)
  map.list <- mkMap(Dat$family, Dat$fixStruct, Dat$rrStruct, Dat$reStruct, dim.list)
  map.names <- names(map.list)
  exp.map.names <- c("thetaf", "ld_rand", "Hg_input", "Hd_input",
                     "ln_kappag", "ln_kappad", "ln_taud", "logit_rhog", "logit_rhod",
                     "ln_sigmaup",  "ln_sigmaep", "ln_sigmau", "ln_sigmav", "upsilon_tg",
                     "epsilon_tjg", "u_ig", "v_ifg", "Gamma_vg"  )

  for(m in 1:length(map.names)){
    expect_equal(sapply(map.list[map.names[m]], dim), sapply(init.parm$parms[map.names[m]], dim))
  }
  expect_equal(sort(map.names), sort(exp.map.names))

  expect_equal(sum(as.numeric(map.list$ln_kappad)), dim.list$n.g*dim.list$n.j)
  expect_equal(sum(is.na(map.list$ln_taud)), dim.list$n.g*dim.list$n.j)

})

