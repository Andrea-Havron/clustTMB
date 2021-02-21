#' Title
#'
#' @param n Simulation size of non-spatial data
#' @param J Dimension of response column
#' @param G Number of clusters
#' @param omega Controls the degree of distribution overlap. See ?MixSim::MixSim variable BarOmega.
#' @param family Distribution family
#' @param space True/False boolean indicating whether function will simulate data with spatial random effects in the gating model component
#' @param spatial.list  List of spatial settings needed to estimate spatial model. Range: distance at which spatial correlation is negligible (~12%), default is the inverse of the no. of clusters; x.range and y.range: min and max coordinates of simulated area; x.incr and y.incr: units of increment in the simulated spatial grid; sp.var.part: variance partition of data
#'
#' @return
#' @export
#'
#' @examples
sim.dat.fn <- function(n = NA, J, G, omega, family, space,
                       spatial.list = list(Range = 1/G,
                                           x.range = c(0,1),
                                           y.range = c(0,1),
                                           x.incr = 0.01,
                                           y.incr = 0.01,
                                           sp.var.part = 0.85))
  {

  Q <- MixSim(BarOmega = omega, K = G, p = J)
  if(family == "Gaussian"){
    if(space == FALSE){
      A <- simdataset(n = n, Pi = Q$Pi, Mu = Q$Mu, S = Q$S)
      dat <- data.frame(A$X, id = A$id)
    }
    if(space == TRUE){
      mu.y <- t(Q$Mu)
      Kappa = sqrt(8)/spatial.list$Range
      var.y <- apply(Q$S,3,diag)
      #Partition variances
      sp.var <- sp.var.part*max(var.y)
      adj.var.y <- abs(var.y - sp.var)
      adj.mu.y <- abs(mu.y - sp.var*2)

      #recalculate covariance matrix for obs error after adjusting for spatial variance
      corr.mat <- cov.mat <- array(0, dim = c(J,J,G))
      for(g in 1:G){
        for(j1 in 1:J){
          for(j2 in 1:J){
            corr.mat[j1,j2,g] <- Q$S[j1,j2,g]/(sqrt(var.y[j1,g])*sqrt(var.y[j2,g]))
            cov.mat[j1,j2,g] <- corr.mat[j1,j2,g] *(sqrt(adj.var.y[j1,g])*sqrt(adj.var.y[j2,g]))
          }
        }
      }

      ## Simulate data outside TMB model ================================
      nu <- 1
      grid.xy = expand.grid(x=seq(spatial.list$x.range[1],spatial.list$x.range[2],spatial.list$x.incr),
                            y = seq(spatial.list$y.range[1],spatial.list$y.range[2],spatial.list$y.incr))
      n.grid <- nrow(grid.xy)
      Sigma <- matrix(NA,n.grid,n.grid)
      d <- as.matrix(dist(grid.xy, upper = TRUE))
      diag(Sigma) <- sp.var
      for(j1 in 1:(n.grid-1)){
        for(j2 in (j1+1):n.grid){
          Sigma[j1,j2] <- sp.var * 2^(1-nu)/gamma(nu) * (Kappa * d[j1,j2]) * besselK(Kappa * d[j1,j2], nu)
          Sigma[j2,j1] <- Sigma[j1,j2]
        }
      }
      omega.g <- matrix(NA, G, n.grid)
      class(omega.g) <- "numeric"
      rmvn(G, rep(0,n.grid), Sigma, ncores = 8,  A = omega.g)
      omega.g <- t(omega.g)
      rm(Sigma);gc()

      pi.g <- matrix(0,n.grid,G)
      X.true <- matrix(0,n.grid,J)
      id.true <- rep(0,n.grid)
      for(j1 in 1:n.grid){
        for(j2 in 1:K){
          pi.g[j1,j2] <- exp(0 + omega.g[j1,j2])/sum(exp(0 + omega.g[j1,]))
        }
        id.true[j1] <- which(pi.g[j1,]==max(pi.g[j1,]))
        X.true[j1,] <- mvtnorm::rmvnorm(1,adj.mu.y[,id.true[j1]], cov.mat[,,id.true[j1]])
      }
      dat <- data.frame(X.true, idx = 1:nrow(grid.xy), x = grid.xy[,1], y = grid.xy[,2],  id = id.true)
    }
  }

  if(family == "Tweedie"){
    #Simulate using copula approach to convert output from MixSim into Compound Poisson Gamma
    if(space == FALSE){
      while(max(apply(Q$S, 3, diag))>10){
        Q <- MixSim(BarOmega = omega, K = G, p = J, int = c(1,5))
      }
      var.mat <- apply(Q$S, 3, diag)
      corr.mat <- Q$S
      for(g in 1:G){
        for(j1 in 1:J){
          for(j2 in 1:J){
            corr.mat[j1,j2,g] <- corr.mat[j1,j2,g]/(sqrt(var.mat[j1,g]*var.mat[j2,g]))
          }
        }
      }

      mu.y <- t(exp(Q$Mu))
      power.y <- matrix(runif(J*G,1.05,1.95),J,G)
      phi.y <- var.mat

      Ng <- rep(1, G) + drop(rmultinom(1, n - G, Q$Pi))
      id <- NULL
      X <- NULL
      for (g in 1:G) {
        id <- c(id, rep(g, Ng[g]))
        x.j <- matrix(NA, Ng[g], J)
        z <- rmvn(Ng[g], rep(0,J), corr.mat[,,g])
        cdf = pnorm(z)
        cdf[which(cdf>0.999999999999999, arr.ind = TRUE)] <- 0.999999999999999 #qtweedie will throw error if any cdf close to 1
        for(j in 1:J){
          x.j[,j] <- qtweedie(cdf[,j], mu = mu.y[j,g],phi = phi.y[j,g], power = power.y[j,g] )
        }
        X <- rbind(X, x.j)
      }

      #apply two transformations: boxcox and log(y+1)
      x.boxcox <- X
      x.log <- X
      for(j in 1:J){
        x.log[,j] <- log(x.log[,j]+1)
      }
      for(j in 1:J){
        boxlambda <- boxcoxfit(x.boxcox[,j], add.to.data = 1)$lambda
        x.boxcox[,j] <- ((x.boxcox[,j] + 1)^boxlambda - 1)/boxlambda
      }

      dat <- data.frame(X, x.boxcox = x.boxcox, x.log = x.log, id = id)
    }
    if(space == TRUE){
      while(max(apply(Q$S, 3, diag))>10){
        Q <- MixSim(BarOmega = omega, K = G, p = J, int = c(1,5))
      }
      mu.y <- t(Q$Mu)
      Range = 1/G
      Kappa = sqrt(8)/Range
      var.y <- apply(Q$S,3,diag)
      #Partition variances
      sp.var <- .85*max(var.y)
      adj.var.y <- abs(var.y - sp.var)
      adj.mu.y <- abs(mu.y - sp.var*2)

      mu.y <- exp(adj.mu.y)*10
      power.y <- matrix(runif(J*G,1.05,1.95),J,G)
      phi.y <- matrix(runif(J*G,1.05,4.95),J,G)

      nu <- 1
      grid.xy = expand.grid(x=seq(0,1,.01), y = seq(0,1,.01))
      n.grid <- nrow(grid.xy)
      Sigma <- matrix(NA,n.grid,n.grid)
      d <- as.matrix(dist(grid.xy, upper = TRUE))
      diag(Sigma) <- sp.var
      for(j1 in 1:(n.grid-1)){
        for(j2 in (j1+1):n.grid){
          Sigma[j1,j2] <- sp.var * 2^(1-nu)/gamma(nu) * (Kappa * d[j1,j2]) * besselK(Kappa * d[j1,j2], nu)
          Sigma[j2,j1] <- Sigma[j1,j2]
        }
      }
      omega.g <- matrix(NA, G, n.grid)
      class(omega.g) <- "numeric"
      rmvn(G, rep(0,n.grid), Sigma, ncores = 8,  A = omega.g)
      omega.g <- t(omega.g)
      rm(Sigma);gc()

      pi.g <- matrix(0,n.grid,G)
      X.true <- matrix(0,n.grid,J)
      id.true <- rep(0,n.grid)
      for(j1 in 1:n.grid){
        for(j2 in 1:G){
          pi.g[j1,j2] <- exp(0 + omega.g[j1,j2])/sum(exp(0 + omega.g[j1,]))
        }
        id.true[j1] <- which(pi.g[j1,]==max(pi.g[j1,]))
      }
      for(g in 1:G){
        idx <- which(id.true == g)
        n.g <- length(idx)
        z = rmvn(n.g, rep(0,J), cov.mat[,,g])
        cdf = pnorm(z)
        cdf[which(cdf>0.999999, arr.ind = TRUE)] <- 0.999999 #qtweedie will throw error if any cdf close to 1
        for(j in 1:J){
          X.true[idx,j] <- qtweedie(cdf[,j], mu = mu.y[j,g], phi = phi.y[j,g], power = power.y[j,g] )
        }
      }
      #apply two transformations: boxcox and log(y+1)
      x.boxcox <- X.true
      x.log <- X.true
      for(j in 1:J){
        x.log[,j] <- log(x.log[,j]+1)
      }
      for(j in 1:J){
        boxlambda <- boxcoxfit(x.boxcox[,j], add.to.data = 1)$lambda
        x.boxcox[,j] <- ((x.boxcox[,j] + 1)^boxlambda - 1)/boxlambda
      }

      dat <- data.frame(X.true, x.boxcox = x.boxcox, x.log = x.log, idx = 1:nrow(grid.xy),
                         x = grid.xy[,1], y = grid.xy[,2],  id = id.true)

    }
  }

 return(dat)

}

