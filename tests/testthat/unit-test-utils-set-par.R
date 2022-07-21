stopifnot(
  require("testthat"),
  require("clustTMB"),
  require("mclust"),
  require("MoEClust")
)

context("test set parameters, no covariates")
data("faithful")
ng <- 2
y <- faithful
test_that("test expmod = FALSE", {
  mod <- Mclust(faithful, G = 2)
  classify <- mod$classification
  Class <- unmap(classify)

  for (g in 1:ng) {
    for (j in 1:ncol(y)) {
      # Subset data by cluster and column
      y.sub <- y[Class[, g] == 1, j]
      Pars <- set.MuVarPow(
        Class. = Class[, g],
        ysub = y.sub,
        expmod = FALSE,
        Xd = matrix(1, length(y.sub), 1),
        family. = gaussian()
      )
      expect_equal(unname(mod$parameters$mean[j, g]),
        Pars$mu_init,
        tolerance = .001
      )
      expect_equal(unname(diag(mod$parameters$variance$sigma[, , g])[j]),
        Pars$var_init,
        tolerance = 0.01
      )

      dat <- list(family = 100)
      Inits <- set.BetaTheta(dat, Pars)
      expect_equal(Pars$mu_init, Inits$betad)
      expect_equal(log(Pars$var_init), Inits$theta)
      expect_equal(Pars$power_init, Inits$thetaf)

      dat <- list(family = 700)
      Inits <- set.BetaTheta(dat, Pars)
      expect_equal(log(Pars$mu_init), Inits$betad)
      power.init <- Pars$power_init
      if (power.init >= 2) power.init <- 1.95
      if (power.init <= 1) power.init <- 1.05
      var.init <- log(Pars$var_init / Pars$mu_init^power.init)
      expect_equal(var.init, Inits$theta)
      expect_equal(power.init, 1 / (1 + exp(-Inits$thetaf)) + 1)

      dat <- list(family = 300)
      Inits <- set.BetaTheta(dat, Pars)
      theta.init <- log(Pars$mu_init^2 / Pars$var_init)
      expect_equal(theta.init, Inits$theta)
    }
  }
})



test_that("test set.MuVarPow, expmod = TRUE", {
  data(ais)
  y <- ais[, 3:7]
  ng <- 3

  mod <- MoE_clust(y,
    G = ng, gating = ~1, expert = ~sex,
    modelNames = c("VEE"), network.data = ais
  )

  classify <- mod$classification
  Class <- unmap(classify)
  for (g in 1:ng) {
    for (j in 1:ncol(y)) {
      # Subset data by cluster and column
      y.sub <- y[Class[, g] == 1, j]
      x.sub <- as.matrix(ais$sex)[Class[, g] == 1, ]
      mod.sub <- glm(y.sub ~ x.sub)
      Pars <- set.MuVarPow(
        Class. = Class[, g],
        ysub = y.sub,
        expmod = TRUE,
        Xd = as.matrix(ais$sex),
        family. = gaussian()
      )
      expect_equal(
        unname(mod.sub$coefficients),
        Pars$mu_init
      )
      expect_equal(
        unname(var(mod.sub$residuals)),
        Pars$var_init
      )
    }
  }
})
