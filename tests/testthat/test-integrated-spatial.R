## Spatial integrated test
## - Uses the sdmTMB package to simulate spatial clusters and generate cluster IDs
## - Uses the MixSim and mvtnorm package to simulate observations based on cluster IDs

test_that("test integrated spatial", {

  library(sdmTMB)
  library(clustTMB)
  library(MixSim)
  library(mvtnorm)
  library(ggplot2)
  library(fmesher)

  ## Simulate spatial clusters
  # generate location data
  set.seed(123)
  Loc <- data.frame(x = runif(500), y = runif(500))

  Mesh.sim <- fmesher::fm_mesh_2d(Loc, max.edge = c(.03, .05))
  Loc.sim <- data.frame(x = Mesh.sim$loc[, 1], y = Mesh.sim$loc[, 2])
  mesh.sim <- make_mesh(
    data = Loc.sim, xy_cols = c("x", "y"),
    mesh = Mesh.sim
  )

  # Generate three spatial fields:
  sim_dat_1 <- sdmTMB_simulate(
    formula = ~1,
    data = Loc.sim,
    mesh = mesh.sim,
    family = gaussian(),
    range = 1 / 3,
    phi = 0.1,
    sigma_O = 1,
    seed = 1,
    B = 1
  )

  sim_dat_2 <- sdmTMB_simulate(
    formula = ~1,
    data = Loc.sim,
    mesh = mesh.sim,
    family = gaussian(),
    range = 1 / 3,
    phi = 0.1,
    sigma_O = 1,
    seed = 2,
    B = 1
  )

  sim_dat_3 <- sdmTMB_simulate(
    formula = ~1,
    data = Loc.sim,
    mesh = mesh.sim,
    family = gaussian(),
    range = 1 / 3,
    phi = 0.1,
    sigma_O = 1,
    seed = 3,
    B = 1
  )

  # Spatial omega values
  omega_clusters <- cbind(sim_dat_1$omega_s, sim_dat_2$omega_s, sim_dat_3$omega_s)
  # Probability vectors
  cluster_pi <- t(apply(omega_clusters, 1, function(x) exp(x) / (sum(exp(x)))))
  # Cluster ID
  cluster_id <- apply(cluster_pi, 1, function(x) which(x == max(x)))

  sim1 <- sim_dat_1[which(cluster_id == 1), ]
  sim1$cluster_id <- 1
  sim2 <- sim_dat_2[which(cluster_id == 2), ]
  sim2$cluster_id <- 2
  sim3 <- sim_dat_3[which(cluster_id == 3), ]
  sim3$cluster_id <- 3

  cluster_dat <- rbind(
    sim1, sim2, sim3
  ) %>% dplyr::select(x, y, omega_s, cluster_id)

  cluster_samp <- cluster_dat[cluster_dat$x > 0 & cluster_dat$x < 1 & cluster_dat$y > 0 & cluster_dat$y < 1, ]
  cluster_samp_ <- cluster_samp[sample(1:nrow(cluster_samp), 500), ]
  Loc <- cluster_samp_[, 1:2]
  Mesh.fit <- fmesher::fm_mesh_2d(Loc, max.edge = c(.1, .5))

  ggplot(cluster_dat, aes(x = x, y = y, color = factor(cluster_id))) +
    geom_point()
  table(cluster_id) / nrow(Loc.sim)
  table(cluster_samp_$cluster_id) / nrow(cluster_samp_)
  # convert locations to spatial coordinates
  sp::coordinates(Loc) <- ~ x * y

  ## Simulate clustered observations
  ## Generate multivariate mean and covariances for
  ## 3 groups and 4 responses using 1% overlap in clusters
  set.seed(123)
  Q <- MixSim::MixSim(BarOmega = .01, K = 3, p = 4)

  # cluster probability from simulated spatial data (sample):
  cluster_counts <- table(cluster_samp_$cluster_id)

  simdat <- data.frame()
  id <- c()
  for (g in 1:3) {
    set.seed(g * 100)
    simdat <- rbind(
      simdat,
      mvtnorm::rmvnorm(
        n = cluster_counts[g],
        mean = Q$Mu[g, ],
        sigma = as.matrix(Q$S[, , g])
      )
    )
    id <- c(id, rep(g, cluster_counts[g]))
  }

  simdat$id <- id

  ## Fit model
  mod <- suppressWarnings(
    clustTMB(
      response = as.matrix(simdat[, 1:4]),
      family = gaussian(),
      gatingformula = ~ gmrf(0 + 1 | loc),
      G = 3, covariance.structure = "VVV",
      spatial.list = list(loc = Loc, mesh = Mesh.fit)
    )
  )

  expect_equal(mod$opt$convergence, 0)
  expect_equal(mod$sdr$pdHess, TRUE)
  expect_equal(mod$obj$gr(mod$opt$par) < .01, rep(TRUE, length(mod$opt$par)))
  # expect classification to be within 1%
  expect_equal(1, MixSim::ClassProp(mod$report$classification, simdat$id),
    tolerance = .01
  )


  # re-order by estimate due to non-identifiability of cluster ID
  beta.df <- as.data.frame(summary(mod$sdr, "fixed")[3:14, ])
  row.names(beta.df) <- NULL
  beta.ord <- beta.df[order(beta.df$Estimate), ]

  beta.ord$true <- sort(Q$Mu)

  expect_equal(
    abs(beta.ord$true - beta.ord$Estimate) <
      qnorm(.975) * beta.ord$`Std. Error`,
    rep(TRUE, nrow(beta.ord))
  )
})
