stopifnot(
  require("testthat"),
  require("clustTMB"),
  require("mclust"),
  require("MoEClust"),
  require("mvnfast"),
  require("FMsmsnReg")
)

# multivariate normal data
data("faithful")


CO2 <- MoEClust::CO2data$CO2
GNP <- MoEClust::CO2data$GNP

## MoEClust code
init.z <- function(y., dat, g, expN, max.init = 1000) {
  n <- nrow(dat)
  z.tmp <- unmap(hclass(hc(dat), G=g))
  n <- nrow(dat)
  z.mat <- z.alloc <- matrix(0L, nrow = n * g, ncol = g)
  muX <- vector("numeric", g)
  
  tmp.z <- matrix(NA, nrow = n, ncol = g)
  mahala <- res.G <- Efit <- list()
  xN <- as.matrix(y.)
  # xN        <- X[!noise,, drop=FALSE]
  # expnoise  <- expx.covs[!noise,, drop=FALSE]
 # expN  <- stats::update.formula(expN, xN ~ .)
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
      dat$sub <- z.tmp[, k] == 1
      exp <- tryCatch(stats::lm(expN, data = dat, subset = sub))
      if (inherits(exp, "try-error")) {
        init.exp <- FALSE
        break
      } else {
        Efit[[k]] <- exp
      }
      #pred.dat <- list(xN, dat)
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

