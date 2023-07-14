stopifnot(
  require("testthat"),
  require("clustTMB"),
  require("mclust"),
  require("MoEClust"),
  require("mvnfast")
)

# multivariate normal data
data("faithful")
clss <- unmap(hclass(hc(faithful), 2))
init.parm.mclust <- list(
  mu = list(
    g1 = apply(faithful[clss[, 1] == 1, ], 2, mean),
    g2 = apply(faithful[clss[, 2] == 1, ], 2, mean)
  ),
  Sigma = list(
    g1 = apply(faithful[clss[, 1] == 1, ], 2, var),
    g2 = apply(faithful[clss[, 2] == 1, ], 2, var)
  ),
  Corr = list(
    g1 = cor(faithful[clss[, 1] == 1, ])[2, 1],
    g2 = cor(faithful[clss[, 2] == 1, ])[2, 1]
  ),
  Pi = list(g1 = 0.5, g2 = 0.5)
)
# Pi = list(g1 = sum(clss[,1])/nrow(faithful),
#           g2 = sum(clss[,2])/nrow(faithful)))
dim.list <- list(
  n.i = nrow(faithful), n.j = ncol(faithful), n.t = 1,
  n.g = 2, n.f.sp = ncol(faithful), n.f.rand = ncol(faithful), n.v = NULL
)
dim.list$nl.fix <- (dim.list$n.j^2 - dim.list$n.j) / 2
dim.list$n.f.rand <- dim.list$n.j
dim.list$nl.rand <- dim.list$n.j * dim.list$n.f.rand -
  (dim.list$n.f.rand * (dim.list$n.f.rand - 1)) / 2
dim.list$n.f.sp <- dim.list$n.j
dim.list$nl.sp <- dim.list$n.j * dim.list$n.f.sp -
  (dim.list$n.f.sp * (dim.list$n.f.sp - 1)) / 2
Dat <- mkDat(
  response = as.matrix(faithful), time.vector = rep(1, nrow(faithful)),
  expert.dat = as.data.frame(rep(1, nrow(faithful))),
  gating.dat = as.data.frame(rep(1, nrow(faithful))),
  expert.formula = ~1,
  gating.formula = ~1,
  ll.method = 0,
  fixStruct = "VVV",
  rrStruct = rep(0, 2),
  reStruct = matrix(0, 2, 3),
  dim.list = dim.list
)
dim.list$n.v <- Dat$spde$n_s
init.parm.clustTMB <- genInit(Dat,
  family = gaussian(link = "identity"), dim.list,
  control = init.options(hc.options = list(
    modelName = "VVV", use = "VARS"
  ))
)

context("equal mvn init parameters")
test_that("mu", {
  expect_equal(
    as.vector(init.parm.mclust$mu$g1),
    init.parm.clustTMB$parms$betad[1, , 1]
  )
  expect_equal(
    as.vector(init.parm.mclust$mu$g2),
    init.parm.clustTMB$parms$betad[1, , 2]
  )
})
test_that("Sigma", {
  expect_equal(
    as.vector(init.parm.mclust$Sigma$g1),
    exp(init.parm.clustTMB$parms$theta[, 1])
  )
  expect_equal(
    as.vector(init.parm.mclust$Sigma$g2),
    exp(init.parm.clustTMB$parms$theta[, 2])
  )
})
# test_that("Corr", { Not same now input is based on cholesky decomp
#   expect_equal(as.vector(unlist(init.parm.mclust$Corr)),
#                as.vector(1/(1+exp(-init.parm.clustTMB$parms$logit_corr_fix))*2-1))
# })
test_that("Pi", {
  expect_equal(
    as.vector(unlist(init.parm.mclust$Pi)),
    c(
      1 / (1 + exp(-init.parm.clustTMB$parms$betag)),
      1 - 1 / (1 + exp(-init.parm.clustTMB$parms$betag))
    )
  )
})
test_that("Class", {
  expect_equal(map(clss), init.parm.clustTMB$class)
})
test_that("Random Effects dim", {
  expect_equal(
    c(dim.list$n.t, dim.list$n.g - 1),
    dim(init.parm.clustTMB$parms$upsilon_tg)
  )
  expect_equal(
    c(dim.list$n.t, dim.list$n.j, dim.list$n.g),
    dim(init.parm.clustTMB$parms$epsilon_tjg)
  )
  # expect_equal(
  # c(dim.list$n.i, dim.list$n.g-1),
  # dim(init.parm.clustTMB$parms$u_ig)
  # )
  expect_equal(
    c(dim.list$n.i, dim.list$n.f.rand, dim.list$n.g),
    dim(init.parm.clustTMB$parms$v_ifg)
  )
  expect_equal(
    c(dim.list$n.v, dim.list$n.g - 1),
    dim(init.parm.clustTMB$parms$Gamma_vg)
  )
  expect_equal(
    c(dim.list$n.v, dim.list$n.f.sp, dim.list$n.g),
    dim(init.parm.clustTMB$parms$Omega_vfg)
  )
})


clss <- unmap(mclust:::qclass(faithful$waiting, 2))
init.parm.mclust <- list(
  mu = list(
    g1 = mean(faithful$waiting[clss[, 1] == 1]),
    g2 = mean(faithful$waiting[clss[, 2] == 1])
  ),
  Sigma = list(
    g1 = var(faithful$waiting[clss[, 1] == 1]),
    g2 = var(faithful$waiting[clss[, 2] == 1])
  ),
  Corr = list(g1 = 0, g2 = 0),
  Pi = c(0.5, 0.5)
)
# Pi = list(g1 = sum(clss[,1])/nrow(faithful),
#           g2 = sum(clss[,2])/nrow(faithful)))
dim.list <- list(
  n.i = nrow(faithful), n.j = 1, n.t = 1,
  n.g = 2, n.f.sp = 1, n.f.rand = 1, n.v = NULL
)

dim.list$nl.fix <- 1
dim.list$nl.rand <- dim.list$n.j * dim.list$n.f.rand -
  (dim.list$n.f.rand * (dim.list$n.f.rand - 1)) / 2
dim.list$nl.sp <- dim.list$n.j * dim.list$n.f.sp -
  (dim.list$n.f.sp * (dim.list$n.f.sp - 1)) / 2
Dat <- mkDat(
  response = as.matrix(faithful$waiting), time.vector = rep(1, nrow(faithful)),
  expert.dat = as.data.frame(rep(1, nrow(faithful))),
  gating.dat = as.data.frame(rep(1, nrow(faithful))),
  expert.formula = ~1,
  gating.formula = ~1,
  ll.method = 0,
  fixStruct = "E",
  rrStruct = rep(0, 2),
  reStruct = matrix(0, 2, 3),
  dim.list = dim.list
)
dim.list$n.v <- Dat$spde$n_s
init.parm.clustTMB <- genInit(Dat,
  family = gaussian(link = "identity"), dim.list,
  control = init.options(init.method = "quantile")
)

context("equal uniNormal init parameters")
test_that("mu", {
  expect_equal(
    as.vector(init.parm.mclust$mu$g1),
    init.parm.clustTMB$parms$betad[1, , 1]
  )
  expect_equal(
    as.vector(init.parm.mclust$mu$g2),
    init.parm.clustTMB$parms$betad[1, , 2]
  )
})
test_that("Sigma", {
  expect_equal(
    as.vector(init.parm.mclust$Sigma$g1),
    exp(init.parm.clustTMB$parms$theta[, 1])
  )
  expect_equal(
    as.vector(init.parm.mclust$Sigma$g2),
    exp(init.parm.clustTMB$parms$theta[, 2])
  )
})
test_that("Corr", {
  expect_equal(
    as.vector(unlist(init.parm.mclust$Corr)),
    as.vector(1 / (1 + exp(-init.parm.clustTMB$parms$logit_corr_fix)) * 2 - 1)
  )
})
test_that("Pi", {
  expect_equal(
    as.vector(unlist(init.parm.mclust$Pi)),
    c(
      1 / (1 + exp(-init.parm.clustTMB$parms$betag)),
      1 - 1 / (1 + exp(-init.parm.clustTMB$parms$betag))
    )
  )
})
test_that("Class", {
  expect_equal(map(clss), init.parm.clustTMB$class)
})


# test MakeADFun
test_that("Random Effects dim", {
  expect_equal(
    c(dim.list$n.t, dim.list$n.g - 1),
    dim(init.parm.clustTMB$parms$upsilon_tg)
  )
  expect_equal(
    c(dim.list$n.t, dim.list$n.j, dim.list$n.g),
    dim(init.parm.clustTMB$parms$epsilon_tjg)
  )
  # expect_equal(
  # c(dim.list$n.i, dim.list$n.g-1),
  # dim(init.parm.clustTMB$parms$u_ig)
  # )
  expect_equal(
    c(dim.list$n.i, dim.list$n.f.rand, dim.list$n.g),
    dim(init.parm.clustTMB$parms$v_ifg)
  )
  expect_equal(
    c(dim.list$n.v, dim.list$n.g - 1),
    dim(init.parm.clustTMB$parms$Gamma_vg)
  )
  expect_equal(
    c(dim.list$n.v, dim.list$n.f.sp, dim.list$n.g),
    dim(init.parm.clustTMB$parms$Omega_vfg)
  )
})


context("univariate normal data with covaraites")
data(CO2data)
CO2 <- CO2data$CO2
GNP <- CO2data$GNP
G <- 2

z.co2 <- unmap(hclass(hc(cbind(CO2, GNP), use = "VARS"), G))
y.co2 <- as.matrix(CO2)
x.co2 <- as.matrix(GNP)
new.z <- clustTMB:::run.mahala(z.co2, y.co2, x.co2)
# MoEclust Mahalanobis distance
## MoEClust code
init.z <- function(y., dat, g, max.init = 1000) {
  expN <- CO2 ~ GNP
  n <- nrow(dat)
  z.tmp <- unmap(hclass(hc(dat, use = "VARS"), g))
  n <- nrow(dat)
  z.mat <- z.alloc <- matrix(0L, nrow = n * g, ncol = g)
  muX <- vector("numeric", g)

  tmp.z <- matrix(NA, nrow = n, ncol = g)
  mahala <- res.G <- Efit <- list()
  xN <- as.matrix(y.)
  # xN        <- X[!noise,, drop=FALSE]
  # expnoise  <- expx.covs[!noise,, drop=FALSE]
  #  expN      <- stats::update.formula(expert, xN ~ .)
  ix <- 0L
  ne <- ncol(dat)
  oldcrit <- Inf
  newcrit <- .Machine$double.xmax
  while (!identical(tmp.z, z.tmp) &&
    newcrit <= oldcrit && ix <= max.init) {
    old.z <- tmp.z
    tmp.z <- z.tmp
    oldcrit <- newcrit
    ix <- ix + 1L
    for (k in 1:g) {
      sub <- z.tmp[, k] == 1
      exp <- tryCatch(stats::lm(expN, data = dat, subset = sub))
      if (inherits(exp, "try-error")) {
        init.exp <- FALSE
        break
      } else {
        Efit[[k]] <- exp
      }
      pred <- tryCatch(suppressWarnings(stats::predict(exp, newdata = dat)))
      if (inherits(pred, "try-error")) {
        init.exp <- FALSE
      } else {
        pred <- as.matrix(pred)

        res <-
          res.G[[k]] <- xN - pred
        mahala[[k]] <- MoEClust::MoE_mahala(exp, res, squared = TRUE, identity = TRUE)
      }
    }

    maha <- do.call(cbind, mahala)
    if (anyNA(maha)) {
      init.exp <- FALSE
      break
    } else {
      mahamin <- apply(maha, 1, min)
      newcrit <- pmin(sum(mahamin), oldcrit)
      z.tmp <- t(apply(maha, 1, function(x) ifelse(x == min(x), 1, 0)))
      if (identical(z.tmp, old.z)) {
        break
      }
    }
  }
  return(z.tmp)
}
adj.z <- init.z(y. = data.frame(CO2 = CO2), dat = data.frame(CO2, GNP), 2)
test_that("class", {
  expect_equal(new.z, adj.z)
})

context("multivariate normal data with covaraites")
data(ais)
hema <- ais[, 3:7]
G <- 3

z.ais <- unmap(hclass(hc(cbind(hema, ais$sex, ais$BMI), use = "VARS"), G))
y.ais <- as.matrix(hema)
x.ais <- cbind(ais$sex, ais$BMI)
new.z <- run.mahala(z.ais, y.ais, x.ais)
# MoEclust Mahalanobis distance
## MoEClust code
init.z <- function(y., dat, g, max.init = 1000) {
  expN <- y.ais ~ ais$sex + ais$BMI
  n <- nrow(dat)
  z.tmp <- unmap(hclass(hc(dat, use = "VARS"), g))
  n <- nrow(dat)
  z.mat <- z.alloc <- matrix(0L, nrow = n * g, ncol = g)
  muX <- vector("numeric", g)

  tmp.z <- matrix(NA, nrow = n, ncol = g)
  mahala <- res.G <- Efit <- list()
  xN <- as.matrix(y.)
  # xN        <- X[!noise,, drop=FALSE]
  # expnoise  <- expx.covs[!noise,, drop=FALSE]
  #  expN      <- stats::update.formula(expert, xN ~ .)
  ix <- 0L
  ne <- ncol(dat)
  oldcrit <- Inf
  newcrit <- .Machine$double.xmax
  while (!identical(tmp.z, z.tmp) &&
    newcrit <= oldcrit && ix <= max.init) {
    old.z <- tmp.z
    tmp.z <- z.tmp
    oldcrit <- newcrit
    ix <- ix + 1L
    for (k in 1:g) {
      sub <- z.tmp[, k] == 1
      exp <- tryCatch(stats::lm(expN, data = dat, subset = sub))
      if (inherits(exp, "try-error")) {
        init.exp <- FALSE
        break
      } else {
        Efit[[k]] <- exp
      }
      pred <- tryCatch(suppressWarnings(stats::predict(exp, newdata = dat)))
      if (inherits(pred, "try-error")) {
        init.exp <- FALSE
      } else {
        pred <- as.matrix(pred)

        res <-
          res.G[[k]] <- xN - pred
        mahala[[k]] <- MoE_mahala(exp, res, squared = TRUE, identity = TRUE)
      }
    }

    maha <- do.call(cbind, mahala)
    if (anyNA(maha)) {
      init.exp <- FALSE
      break
    } else {
      mahamin <- apply(maha, 1, min)
      newcrit <- pmin(sum(mahamin), oldcrit)
      z.tmp <- t(apply(maha, 1, function(x) ifelse(x == min(x), 1, 0)))
      if (identical(z.tmp, old.z)) {
        break
      }
    }
  }
  return(z.tmp)
}
adj.z <- init.z(
  y. = as.matrix(hema),
  dat = data.frame(y.ais, ais$sex, ais$BMI), 3
)
test_that("class", {
  expect_equal(new.z, adj.z)
})


# test_that("mu", {
#   expect_equal(
# as.vector(init.parm.mclust$mu$g1),
# init.parm.clustTMB$parms$betad[1,,1]
# )
#   expect_equal(
# as.vector(init.parm.mclust$mu$g2),
# init.parm.clustTMB$parms$betad[1,,2]
# )
# })
# test_that("Sigma", {
#   expect_equal(
# as.vector(init.parm.mclust$Sigma$g1),
# exp(init.parm.clustTMB$parms$theta[,1])
# )
#   expect_equal(
# as.vector(init.parm.mclust$Sigma$g2),
# exp(init.parm.clustTMB$parms$theta[,2])
# )
# })
# test_that("Corr", {
#   expect_equal(
# as.vector(unlist(init.parm.mclust$Corr)),
# as.vector(1/(1+exp(-init.parm.clustTMB$parms$logit_corr_fix))*2-1)
# )
# })
# test_that("Pi", {
#   expect_equal(as.vector(unlist(init.parm.mclust$Pi)),
#                c(1/(1+exp(-init.parm.clustTMB$parms$betag)),
#                  1-1/(1+exp(-init.parm.clustTMB$parms$betag))))
# })
# test_that("Class", {
#   expect_equal(map(clss), init.parm.clustTMB$class)
# })
