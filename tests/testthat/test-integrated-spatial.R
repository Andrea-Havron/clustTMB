## Spatial integrated test 
## - Uses the sdmTMB package to simulate spatial clusters and generate cluster IDs
## - Uses the MixSim and mvtnorm package to simulate observations based on cluster IDs


library(sdmTMB)
library(clustTMB)
library(MixSim)
library(mvtnorm)
library(ggplot2)

## Simulate spatial clusters
# generate location data
set.seed(123)
Loc <- data.frame(x = runif(300), y = runif(300))

Mesh <- INLA::inla.mesh.2d(Loc, max.edge = c(.1,.5))
mesh <- make_mesh(data = Loc, xy_cols = c("x", "y"), mesh = Mesh)

# Generate three spatial fields:
sim_dat_1 <- sdmTMB_simulate(
  formula = ~ 1,
  data = Loc,
  mesh = mesh,
  family = gaussian(),
  range = 2/3,
  phi = 0.1,
  sigma_O = 1,
  seed = 11,
  B = 1
)

sim_dat_2 <- sdmTMB_simulate(
  formula = ~ 1,
  data = Loc,
  mesh = mesh,
  family = gaussian(),
  range = 2/3,
  phi = 0.1,
  sigma_O = 1,
  seed = 22,
  B = 1
)

sim_dat_3 <- sdmTMB_simulate(
  formula = ~ 1,
  data = Loc,
  mesh = mesh,
  family = gaussian(),
  range = 2/3,
  phi = 0.1,
  sigma_O = 1,
  seed = 33,
  B = 1
)

#Spatial omega values
omega_clusters <- cbind(sim_dat_1$omega_s, sim_dat_2$omega_s, sim_dat_3$omega_s)
#Probability vectors
cluster_pi <- t(apply(omega_clusters, 1, function(x) exp(x)/(sum(exp(x)))))
#Cluster ID
cluster_id <- apply(cluster_pi, 1, function(x) which(x == max(x)))

sim1 <- sim_dat_1[which(cluster_id == 1),]
sim1$cluster_id <- 1
sim2 <- sim_dat_2[which(cluster_id == 2),]
sim2$cluster_id <- 2
sim3 <- sim_dat_3[which(cluster_id == 3),]
sim3$cluster_id <- 3

cluster_dat <-  rbind(
  sim1, sim2, sim3
) %>% dplyr::select(x, y, omega_s, cluster_id)

ggplot(cluster_dat, aes(x=x, y=y, color = factor(cluster_id))) + geom_point()

#convert locations to spatial coordinates
sp::coordinates(Loc) <- ~ x*y

## Simulate clustered observations
## Generate multivariate mean and covariances for 
## 3 groups and 4 responses using 1% overlap in clusters
set.seed(42)
Q <- MixSim::MixSim(BarOmega = .01, K = 3, p = 4)

#cluster probability from simulated spatial data:
cluster_counts <- table(cluster_id)

set.seed(123)
simdat <- data.frame()
id <- c()
for(g in 1:3){
  simdat <- rbind(simdat,
                  mvtnorm::rmvnorm(n = cluster_counts[g], 
                                   mean = Q$Mu[g,], 
                                   sigma = as.matrix(Q$S[,,g]))
                  )
  id <- c(id, rep(g, cluster_counts[g]))
}

simdat$id <- id

## Fit model
mod <- clustTMB(response = as.matrix(simdat[,1:4]),
                family = gaussian(), 
                gatingformula = ~ gmrf(0 + 1|loc),
                G = 3, covariance.structure = "VVV",
                spatial.list = list(loc = Loc, mesh = Mesh)
                )
test_that("test convergence", {
  expect_equal(mod$opt$convergence, 0)
  expect_equal(mod$sdr$pdHess, TRUE)
  expect_equal(mod$obj$gr(mod$opt$par)< .01, rep(TRUE, length(mod$opt$par)))
})

test_that("test classification",{
  #expect classification to be within 1%
  expect_equal(1, MixSim::ClassProp(mod$report$classification, simdat$id),
               tolerance = .01)
})

test_that("test mu",{
  beta1 <- summary(mod$sdr, "fixed")[3:6,]
  beta2 <- summary(mod$sdr, "fixed")[11:14,]
  beta3 <- summary(mod$sdr, "fixed")[7:10,]
  
  expect_equal(unname(abs(Q$Mu[1,] - beta1[,1]) < qnorm(.975)*beta1[,2]), rep(TRUE, 4))
  expect_equal(unname(abs(Q$Mu[2,] - beta2[,1]) < qnorm(.975)*beta2[,2]), rep(TRUE, 4))
  expect_equal(unname(abs(Q$Mu[3,] - beta3[,1]) < qnorm(.975)*beta3[,2]), rep(TRUE, 4))
})

#Add test on covariance



