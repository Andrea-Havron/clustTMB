stopifnot(
  require("testthat"),
  require("clustTMB"),
  require("mvnfast")
)

context("match mkMap with init parm")
n.i <- 100
n.j <- 4
n.g <- 3
y <- rmvn(n.i, rep(0, n.j), diag(n.j))

dim.list <- list(
  n.i = n.i, n.j = n.j, n.t = 1,
  n.g = n.g, n.f.sp = n.j, n.f.rand = n.j, n.v = NULL
)
dim.list$nl.fix <- (dim.list$n.j^2 - dim.list$n.j) / 2
dim.list$nl.rand <- dim.list$n.j * dim.list$n.f.rand -
  (dim.list$n.f.rand * (dim.list$n.f.rand - 1)) / 2
dim.list$nl.sp <- dim.list$n.j * dim.list$n.f.sp -
  (dim.list$n.f.sp * (dim.list$n.f.sp - 1)) / 2


test_that("mvn with no random effect, no rank reduction", {
  Dat <- mkDat(
    response = as.matrix(y), time.vector = rep(1, dim.list$n.i),
    expert.dat = as.data.frame(rep(1, dim.list$n.i)),
    gating.dat = as.data.frame(rep(1, dim.list$n.i)),
    expert.formula = ~1,
    gating.formula = ~1,
    ll.method = 0,
    fixStruct = "VVV",
    rrStruct = rep(0, 2),
    reStruct = matrix(0, 2, 3),
    dim.list = dim.list
  )
  dim.list$n.v <- Dat$spde$n_s
  init.parm <- genInit(Dat, family = gaussian(link = "identity"), dim.list)
  map.list <- mkMap(
    Dat$family,
    Dat$fixStruct,
    Dat$rrStruct,
    Dat$reStruct,
    dim.list
  )
  map.names <- names(map.list)
  exp.map.names <- c(
    "thetaf", "ld_rand", "ld_sp", "Hg_input", "Hd_input",
    "ln_kappag", "ln_kappad", "ln_taud", "logit_rhog", "logit_rhod",
    "ln_sigmaup", "ln_sigmaep", "ln_sigmav", "upsilon_tg",
    "epsilon_tjg", "v_ifg", "Gamma_vg", "Omega_vfg"
  ) # rm ln_sigmau, u_ig
  for (m in seq_along(map.names)) {
    mapparm <- unname(unlist(map.list[map.names[m]]))
    initparm <- unname(unlist(init.parm$parms[map.names[m]]))
    expect_equal(
      vapply(mapparm, length, rep(0, 1)),
      vapply(initparm, length, rep(0, 1))
    )
  }
  expect_equal(sort(map.names), sort(exp.map.names))
})

test_that("mvn with gating random effects, no rank reduction", {
  Dat <- mkDat(
    response = as.matrix(y), time.vector = rep(1, dim.list$n.i),
    expert.dat = as.data.frame(rep(1, dim.list$n.i)),
    gating.dat = as.data.frame(rep(1, dim.list$n.i)),
    expert.formula = ~1,
    gating.formula = ~1,
    ll.method = 0,
    fixStruct = "VVV",
    rrStruct = rep(0, 2),
    reStruct = matrix(0, 2, 3),
    dim.list = dim.list
  )
  dim.list$n.v <- Dat$spde$n_s
  exp.map.names <- c(
    "thetaf", "ld_rand", "ld_sp", "Hg_input", "Hd_input",
    "ln_kappag", "ln_kappad", "ln_taud", "logit_rhog", "logit_rhod",
    "ln_sigmaup", "ln_sigmaep", "ln_sigmav", "upsilon_tg",
    "epsilon_tjg", "v_ifg", "Gamma_vg", "Omega_vfg"
  ) # rm ln_sigmau, u_ig
  reNum <- c(3, 2, 1)
  parmName <- list(
    "Gamma_vg",
    c("logit_rhog", "ln_sigmaup", "upsilon_tg"),
    c("ln_sigmau", "u_ig")
  )
  for (j in 1:3) {
    Dat$reStruct <- matrix(0, 2, 3)
    Dat$reStruct[1, j] <- reNum[j]
    init.parm <- genInit(Dat, family = gaussian(link = "identity"), dim.list)
    map.list <- mkMap(
      Dat$family,
      Dat$fixStruct,
      Dat$rrStruct,
      Dat$reStruct,
      dim.list
    )
    map.names <- names(map.list)
    exp.map.names <- c(
      "thetaf", "ld_rand", "ld_sp", "Hg_input", "Hd_input",
      "ln_kappag", "ln_kappad", "ln_taud", "logit_rhog", "logit_rhod",
      "ln_sigmaup", "ln_sigmaep", "ln_sigmav", "upsilon_tg",
      "epsilon_tjg", "v_ifg", "Gamma_vg", "Omega_vfg"
    ) # rm ln_sigmau, u_ig
    for (i in seq_along(parmName[[j]])) {
      exp.map.names <- exp.map.names[exp.map.names != parmName[[j]][i]]
    }

    for (m in seq_along(map.names)) {
      mapparm <- unname(unlist(map.list[map.names[m]]))
      initparm <- unname(unlist(init.parm$parms[map.names[m]]))
      expect_equal(
        vapply(mapparm, length, rep(0, 1)),
        vapply(initparm, length, rep(0, 1))
      )
    }
    expect_equal(sort(map.names), sort(exp.map.names))
    if (j == 1) {
      expect_equal(sum(as.numeric(map.list$ln_kappag)), dim.list$n.g - 1)
    }
  }
})

test_that("mvn with expert random effects, no rank reduction", {
  Dat <- mkDat(
    response = as.matrix(y), time.vector = rep(1, dim.list$n.i),
    expert.dat = as.data.frame(rep(1, dim.list$n.i)),
    gating.dat = as.data.frame(rep(1, dim.list$n.i)),
    expert.formula = ~1,
    gating.formula = ~1,
    ll.method = 0,
    fixStruct = "VVV",
    rrStruct = rep(0, 2),
    reStruct = matrix(0, 2, 3),
    dim.list = dim.list
  )
  dim.list$n.v <- Dat$spde$n_s
  reNum <- c(3, 2, 1)
  parmName <- list(
    "Omega_vfg",
    c("logit_rhod", "ln_sigmaep", "epsilon_tjg"),
    c("ln_sigmav", "v_ifg")
  )
  for (j in 1:3) {
    Dat$reStruct <- matrix(0, 2, 3)
    Dat$reStruct[2, j] <- reNum[j]
    init.parm <- genInit(Dat, family = gaussian(link = "identity"), dim.list)
    map.list <- mkMap(
      Dat$family,
      Dat$fixStruct,
      Dat$rrStruct,
      Dat$reStruct,
      dim.list
    )
    map.names <- names(map.list)
    exp.map.names <- c(
      "thetaf", "ld_rand", "ld_sp", "Hg_input", "Hd_input",
      "ln_kappag", "ln_kappad", "ln_taud", "logit_rhog", "logit_rhod",
      "ln_sigmaup", "ln_sigmaep", "ln_sigmav", "upsilon_tg",
      "epsilon_tjg", "v_ifg", "Gamma_vg", "Omega_vfg"
    ) # rm ln_sigmau, u_ig
    for (i in seq_along(parmName[[j]])) {
      exp.map.names <- exp.map.names[exp.map.names != parmName[[j]][i]]
    }

    for (m in seq_along(map.names)) {
      mapparm <- unname(unlist(map.list[map.names[m]]))
      initparm <- unname(unlist(init.parm$parms[map.names[m]]))
      expect_equal(
        vapply(mapparm, length, rep(0, 1)),
        vapply(initparm, length, rep(0, 1))
      )
    }
    expect_equal(sort(map.names), sort(exp.map.names))
    if (j == 1) {
      expect_equal(
        sum(as.numeric(map.list$ln_kappad)),
        dim.list$n.g * dim.list$n.j
      )
      expect_equal(
        sum(as.numeric(map.list$ln_taud)),
        dim.list$n.g * dim.list$n.j
      )
    }
  }
})

test_that("mvn with expert random effects and rank reduction", {
  ## random reduction
  dim.list <- list(
    n.i = n.i, n.j = n.j, n.t = 1,
    n.g = n.g, n.f.sp = n.j, n.f.rand = n.j - 1, n.v = NULL
  )

  dim.list$nl.fix <- (dim.list$n.j^2 - dim.list$n.j) / 2
  dim.list$nl.rand <- dim.list$n.j * dim.list$n.f.rand -
    (dim.list$n.f.rand * (dim.list$n.f.rand - 1)) / 2
  dim.list$nl.sp <- dim.list$n.j * dim.list$n.f.sp -
    (dim.list$n.f.sp * (dim.list$n.f.sp - 1)) / 2


  Dat <- mkDat(
    response = as.matrix(y), time.vector = rep(1, dim.list$n.i),
    expert.dat = as.data.frame(rep(1, dim.list$n.i)),
    gating.dat = as.data.frame(rep(1, dim.list$n.i)),
    expert.formula = ~1,
    gating.formula = ~1,
    ll.method = 0,
    fixStruct = "VVI",
    rrStruct = rep(0, 2),
    reStruct = matrix(0, 2, 3),
    dim.list = dim.list
  )
  dim.list$n.v <- Dat$spde$n_s

  Dat$reStruct <- matrix(0, 2, 3)
  Dat$reStruct[2, 3] <- 1
  Dat$rrStruct[1] <- 1
  init.parm <- genInit(Dat, family = gaussian(link = "identity"), dim.list)
  map.list <- mkMap(
    Dat$family,
    Dat$fixStruct,
    Dat$rrStruct,
    Dat$reStruct,
    dim.list
  )
  map.names <- names(map.list)
  exp.map.names <- c(
    "thetaf", "ld_sp", "Hg_input", "Hd_input", "logit_corr_fix",
    "ln_kappag", "ln_kappad", "ln_taud", "logit_rhog", "logit_rhod",
    "ln_sigmaup", "ln_sigmaep", "ln_sigmav", "upsilon_tg",
    "epsilon_tjg", "Gamma_vg", "Omega_vfg"
  ) # rm ln_sigmau, u_ig

  for (m in seq_along(map.names)) {
    mapparm <- unname(unlist(map.list[map.names[m]]))
    initparm <- unname(unlist(init.parm$parms[map.names[m]]))
    expect_equal(
      vapply(mapparm, length, rep(0, 1)),
      vapply(initparm, length, rep(0, 1))
    )
  }
  expect_equal(sort(map.names), sort(exp.map.names))

  expect_equal(
    sum(is.na(map.list$ln_sigmav)),
    dim.list$n.g * dim.list$n.f.rand
  )

  expect_equal(
    dim(init.parm$parms$v_ifg),
    c(dim.list$n.i, dim.list$n.f.rand, dim.list$n.g)
  )

  expect_equal(
    dim(init.parm$parms$ld_rand)[1],
    c(dim.list$n.j * dim.list$n.f.rand -
      (dim.list$n.f.rand * (dim.list$n.f.rand - 1)) / 2)
  )

  ## spatial reduction
  dim.list <- list(
    n.i = n.i, n.j = n.j, n.t = 1,
    n.g = n.g, n.f.sp = n.j - 1, n.f.rand = n.j, n.v = NULL
  )

  dim.list$nl.fix <- (dim.list$n.j^2 - dim.list$n.j) / 2
  dim.list$nl.rand <- dim.list$n.j * dim.list$n.f.rand -
    (dim.list$n.f.rand * (dim.list$n.f.rand - 1)) / 2
  dim.list$nl.sp <- dim.list$n.j * dim.list$n.f.sp -
    (dim.list$n.f.sp * (dim.list$n.f.sp - 1)) / 2

  Dat <- mkDat(
    response = as.matrix(y), time.vector = rep(1, dim.list$n.i),
    expert.dat = as.data.frame(rep(1, dim.list$n.i)),
    gating.dat = as.data.frame(rep(1, dim.list$n.i)),
    expert.formula = ~1,
    gating.formula = ~1,
    ll.method = 0,
    fixStruct = "VVV",
    rrStruct = rep(0, 2),
    reStruct = matrix(0, 2, 3),
    dim.list = dim.list
  )
  dim.list$n.v <- Dat$spde$n_s
  Dat$reStruct <- matrix(0, 2, 3)
  Dat$reStruct[2, 1] <- 3
  Dat$rrStruct[2] <- 1
  init.parm <- genInit(Dat, family = gaussian(link = "identity"), dim.list)
  map.list <- mkMap(
    Dat$family,
    Dat$fixStruct,
    Dat$rrStruct,
    Dat$reStruct,
    dim.list
  )
  map.names <- names(map.list)
  exp.map.names <- c(
    "thetaf", "ld_rand", "Hg_input", "Hd_input",
    "ln_kappag", "ln_kappad", "ln_taud", "logit_rhog", "logit_rhod",
    "ln_sigmaup", "ln_sigmaep", "ln_sigmav", "upsilon_tg",
    "epsilon_tjg", "v_ifg", "Gamma_vg"
  ) # rm ln_sigmau, u_ig

  for (m in seq_along(map.names)) {
    mapparm <- unname(unlist(map.list[map.names[m]]))
    initparm <- unname(unlist(init.parm$parms[map.names[m]]))
    expect_equal(
      vapply(mapparm, length, rep(0, 1)),
      vapply(initparm, length, rep(0, 1))
    )
  }
  expect_equal(sort(map.names), sort(exp.map.names))

  expect_equal(
    sum(as.numeric(map.list$ln_kappad)),
    dim.list$n.g * dim.list$n.f.sp
  )
  expect_equal(
    sum(is.na(map.list$ln_taud)),
    dim.list$n.g * dim.list$n.f.sp
  )

  expect_equal(
    dim(init.parm$parms$Omega_vfg),
    c(dim.list$n.v, dim.list$n.f.sp, dim.list$n.g)
  )

  expect_equal(
    dim(init.parm$parms$ld_sp)[1],
    c(dim.list$n.j * dim.list$n.f.sp -
      (dim.list$n.f.sp * (dim.list$n.f.sp - 1)) / 2)
  )
})

context("test user specified start/map")
test_that("map tests", {
  expect_error(clustTMB(y,
    covariance.structure = "VVV",
    G = n.g, Map = list(Sigma = runif(10, 0, 1))
  ))
  expect_error(clustTMB(y,
    covariance.structure = "VVV",
    G = n.g, Map = list(thetaf = runif(10, 0, 1))
  ))
  map.thetaf <- factor(matrix(1.6, n.j, n.g))
  expect_error(clustTMB(log(y - min(y) + 1),
    family = tweedie(link = "log"),
    covariance.structure = "VII", G = n.g, Map = list(thetaf = map.thetaf)
  ))
  dim(map.thetaf) <- c(n.j, n.g)
  mod <- clustTMB(log(y - min(y) + 1),
    family = tweedie(link = "log"),
    covariance.structure = "VII", G = n.g, Map = list(thetaf = map.thetaf),
    control = run.options(check.input = TRUE)
  )
  expect_equal(map.thetaf, mod$map$thetaf)
  mod <- clustTMB(log(y - min(y) + 1),
    family = tweedie(link = "log"),
    covariance.structure = "EII", G = n.g, Map = list(thetaf = map.thetaf)
  )
  expect_equal(n.j * n.g + 1 + n.g - 1 + 1, length(mod$opt$par))
  mod <- clustTMB(log(y - min(y) + 1),
    family = tweedie(link = "log"),
    covariance.structure = "VII", G = n.g, Map = list(thetaf = map.thetaf)
  )
  expect_equal(n.j * n.g + n.g + n.g - 1 + 1, length(mod$opt$par))
  mod <- clustTMB(log(y - min(y) + 1),
    family = tweedie(link = "log"),
    covariance.structure = "EEI", G = n.g, Map = list(thetaf = map.thetaf)
  )
  expect_equal(n.j * n.g + n.j + n.g - 1 + 1, length(mod$opt$par))
  mod <- clustTMB(log(y - min(y) + 1),
    family = tweedie(link = "log"),
    covariance.structure = "VVI", G = n.g, Map = list(thetaf = map.thetaf)
  )
  expect_equal(n.j * n.g * 2 + n.g - 1 + 1, length(mod$opt$par))
  init.thetaf <- matrix(1.6, n.j, n.g)
  map.thetaf <- factor(matrix(NA, n.j, n.g))
  dim(map.thetaf) <- c(n.j, n.g)
  mod <- clustTMB(log(y - min(y) + 1),
    family = tweedie(link = "log"),
    covariance.structure = "VII", G = n.g,
    Start = list(thetaf = init.thetaf), Map = list(thetaf = map.thetaf)
  )
  expect_equal(init.thetaf, attributes(mod$obj$env$parameters$thetaf)$shape)
})
test_that("start tests", {
  expect_error(clustTMB(y,
    covariance.structure = "VVV",
    G = n.g, Start = list(theta = matrix(1.6, n.j, n.g))
  ))
  # FIXME: expect_condition(
  #   try(
  #     clustTMB(
  #       y,
  #       covariance.structure = 'VVV',
  #       G = n.g,start = list(u_ig = matrix(rnorm(n.i*(n.g-1)),n.i,n.g-1))
  #       )
  #     )
  #   )
  expect_error(clustTMB(log(y - min(y) + 1),
    family = tweedie(link = "log"), G = n.g,
    Start = list(thetaf = rep(1.6, 3)),
    covariance.structure = "VII"
  ))
  init.thetaf <- matrix(1.6, n.j, n.g)
  mod <- clustTMB(log(y - min(y) + 1),
    family = tweedie(link = "log"), G = n.g,
    Start = list(thetaf = init.thetaf), covariance.structure = "VII",
    control = run.options(check.input = TRUE)
  )
  expect_equal(init.thetaf, mod$inits$parms$thetaf)
})
# test_that('incorrect parm dimension', {
#   expect_error()
# })
context("tweedie tests")
test_that("tweedie covariance structure", {
  expect_error(
    clustTMB(log(y - min(y) + 1),
      family = tweedie(link = "log"),
      covariance.structure = "VVV"
    )
  )
})
