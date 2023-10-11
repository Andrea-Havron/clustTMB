#' genInit: generates initial values for cluster model
#'
#' @param Data List containing TMB data objects
#' @param family Distribution family
#' @param dim.list List of model dimensions
#' @param control Calls init.options() to generate settings for initial values. Arguments of init.options() can be specified by the user.
#' 1. init.method - Single character string indicating initial clustering method. Methods include: hc, quantile, random, mclust, kmeans, mixed, user. Defaults to 'hc'. In the case where data are univariate and there are no covariates in the gating/expert formula, this defaults to 'quantile'
#' 2. hc.options - Named list of two character strings specifying hc modelName and hcUse when init.method = 'hc'. The default modelName is 'VVV' and the default use is 'SVD' unless gating/expert covariates specified, in which case default in VARS. See ?mclust::mclust.options for complete list of options.
#' 3. mix.method - String stating initialization method for mixed-type data (init.method = 'mixed'). Current default when Tweedie family specified. Options include: Gower kmeans (default), Gower hclust, and kproto.
#' 4. user - Numeric or character vector defining user specified initial classification. init.method must be set to 'user' when using this option.
#'
#' @importFrom stats cutree gaussian hclust runif rmultinom nlminb cor glm var
#' @importFrom mclust unmap hclass hc hcVVV hcE hcEII hcEEE hcVII hcV
#' @importFrom cluster daisy pam
#' @importFrom clustMixType kproto
#'
#' @return list
#' @keywords internal
#' @noRd
genInit <- function(Data, family = NULL, dim.list, control = init.options()) {
  # list2env(dim.list, environment(genStart))
  n.i <- dim.list$n.i
  n.j <- dim.list$n.j
  n.t <- dim.list$n.t
  n.k.g <- ncol(Data$Xg)
  n.k.e <- ncol(Data$Xd)
  n.g <- dim.list$n.g
  n.f.rand <- dim.list$n.f.rand
  n.f.sp <- dim.list$n.f.sp
  n.f <- dim.list$n.f
  n.v <- dim.list$n.v
  nl.fix <- dim.list$nl.fix
  nl.rand <- dim.list$nl.rand
  nl.sp <- dim.list$nl.sp

  # Apply any data transformations
  y <- Data$Y
  if (Data$family == 300 | Data$family == 600) {
    y <- log(y)
  }
  X.g <- subset(Data$Xg, select = colnames(Data$Xg) != ("(Intercept)"))
  X.d <- subset(Data$Xd, select = colnames(Data$Xd) != ("(Intercept)"))
  gate.mod <- exp.mod <- FALSE
  if (ncol(X.g) > 0) {
    gate.mod <- TRUE
    y <- cbind(y, X.g)
  }
  if (ncol(X.d) > 0) {
    exp.mod <- TRUE
    y <- cbind(y, X.d)
  }

  ## reset defaults based on family, dimension, and expert/gating models
  control <- reset.defaults(Data$family, control, gate.mod, exp.mod, n.j)

  y <- y[, !duplicated(t(y))]

  # Apply classification method
  classify <- genInitMethods(n.g, n.i, n.j, control, y)
  Class <- unmap(classify)
  pi.init <- apply(Class, 2, sum) / n.i

  # setup ParList
  ParList <- list(
    betag = matrix(0, n.k.g, (n.g - 1)),
    betad = array(0, dim = c(n.k.e, n.j, n.g)),
    betapz = array(numeric(0)),
    theta = matrix(0, n.j, n.g),
    thetaf = matrix(0, n.j, n.g),
    logit_corr_fix = matrix(0, nl.fix, n.g),
    ld_rand = matrix(0, nl.rand, n.g),
    ld_sp = matrix(0, nl.sp, n.g),
    Hg_input = matrix(0, 2, (n.g - 1)),
    Hd_input = array(0, dim = c(2, n.f.sp, n.g)),
    ln_kappag = rep(0, (n.g - 1)),
    ln_kappad = matrix(0, n.f.sp, n.g),
    ln_taud = matrix(0, n.f.sp, n.g),
    logit_rhog = rep(0, (n.g - 1)),
    logit_rhod = matrix(0, n.j, n.g),
    ln_sigmaup = rep(0, (n.g - 1)),
    ln_sigmaep = matrix(0, n.j, n.g),
    # ln_sigmau = rep(0,(n.g-1)),
    ln_sigmav = matrix(0, n.f.rand, n.g),
    upsilon_tg = array(0, dim = c(n.t, (n.g - 1))),
    epsilon_tjg = array(0, dim = c(n.t, n.j, n.g)),
    # u_ig = array(0, dim = c(n.i,(n.g-1))),
    v_ifg = array(0, dim = c(n.i, n.f.rand, n.g)),
    Gamma_vg = array(0, dim = c(n.v, (n.g - 1))),
    Omega_vfg = array(0, dim = c(n.v, n.f.sp, n.g))
  )
  # if(sum(Data$reStruct[1,])==0 & !gate.mod){
  #   ParList$betag[1,] <- matrix(log(pi.init[1:(n.g-1)]/(1 - sum(pi.init[1:(n.g-1)]))), nrow = 1)
  # }
  # if(sum(Data$reStruct[1,])==0 & gate.mod){
  #   mod <- multinom(Class~X.g)
  #   for(g in 1:(n.g-1)){
  #     #nnet flips labels, make coefficients negative to correct
  #     ParList$betag[,g] <- -as.vector(t(summary(mod)$coefficients)[,g])
  #   }
  # }
  #
  # if(sum(Data$reStruct[1,])>0 & gate.mod){
  #   mod <- multinom(Class~X.g-1)
  #   for(g in 1:(n.g-1)){
  #     #nnet flips labels, make coefficients negative to correct
  #     ParList$betag[,g] <- -c(0,as.vector(t(summary(mod)$coefficients)[,g]))
  #   }
  # }

  # Set initial values for kappa and tau if n.r provided
  if (Data$reStruct[1, 1] > 2 & !is.null(dim.list$n.r.g)) {
    ParList$ln_kappag <- rep(log(sqrt(8) / (dim.list$n.r.g / 2)), (n.g - 1))
  }
  if (Data$reStruct[2, 1] > 2 & !is.null(dim.list$n.r.e)) {
    ParList$ln_kappad <- matrix(
      log(sqrt(8) / (dim.list$n.r.e / 2)),
      n.j,
      n.g
    )
    ParList$ln_taud <- matrix(
      1 / (2 * sqrt(pi) * sqrt(8) / (dim.list$n.r.e / 2)),
      n.j,
      n.g
    )
  }

  # re-intitiate data for initial value estimation
  y <- Data$Y
  if (Data$family == 300 | Data$family == 600) {
    y <- log(y)
  }

  # Update initial values based on classification
  res.mat <- matrix(0, n.i, n.j)
  if (exp.mod & control$exp.init$mahala) {
    Class <- run.mahala(Class, as.matrix(y), as.matrix(X.d))
  }
  if (exp.mod & any(apply(Class, 2, sum) == 0)) {
    stop("initalization method results in an empty or unit cluster
         which is not suitable when intializing the expert model")
  }

  # Set initial values in ParList
  for (g in 1:n.g) {
    for (j in 1:n.j) {
      # Subset data by cluster and column
      y.sub <- y[Class[, g] == 1, j]
      Stats <- set.MuVarPow(Class[, g], y.sub, exp.mod, X.d, family)
      mu.init <- Stats$mu_init
      var.init <- Stats$var_init
      power.init <- Stats$power_init
      if (exp.mod) {
        res.mat[Class[, g] == 1, j] <- Stats$residuals
      }

      Inits <- set.BetaTheta(Data, Stats)
      ParList$betad[, j, g] <- Inits$betad
      ParList$theta[j, g] <- Inits$theta
      ParList$thetaf[j, g] <- Inits$thetaf
    }

    y.mat <- y[Class[, g] == 1, ]
    res <- res.mat[Class[, g] == 1, ]

    if (Data$fixStruct != 10) {
      cor.mat <- cor(y.mat)
      if (exp.mod) {
        cor.mat <- cor(res)
      }
      # Apply correction if NA in cor.mat
      cor.mat <- cormat.correction(cor.mat, y.mat, n.j)

      corvec <- cor.mat[lower.tri(cor.mat)]

      Loadings <- set.Loadings(Data, cor.mat, corvec, dim.list)
      ParList$logit_corr_fix[, g] <- Loadings$logit_corr_fix
      ParList$ld_rand[, g] <- Loadings$ld_rand
      ParList$ld_sp[, g] <- Loadings$ld_sp
    }
  } # end g loop

  gen.init <- list(parms = ParList, class = classify)
  return(gen.init)
}

# genInit helper functions

#' Apply classification method dependent on init.method
#'
#' @param n.g Number of clusters
#' @param n.i Number of observations
#' @param n.j Number of columns
#' @param control Classification settings from init.options()
#' @param y Observations
#'
#' @return classification vector
#' @noRd
genInitMethods <- function(n.g, n.i, n.j,
                           control, y) {
  # Apply classification method
  if (control$init.method == "random") {
    classify <- sample(1:n.g, n.i, replace = TRUE)
  }

  if (control$init.method == "quantile") {
    classify <- mc.qclass(y, as.numeric(n.g))
  }

  if (control$init.method == "hc") {
    classify <- as.vector(hclass(
      hc(y,
        modelName = control$hc.options$modelName,
        use = control$hc.options$use
      ), n.g
    ))
  }

  if (control$init.method == "kmeans") {
    classify <- cluster::pam(y, k = n.g)$clustering
  }

  if (control$init.method == "mixed") {
    tmp.pa <- matrix(NA, n.i, n.j)
    y1 <- as.data.frame(matrix(NA, n.i, n.j))
    for (j in 1:n.j) {
      tmp.pa[, j] <- ifelse(y[, j] == 0, "A", "P")
      y1[, j] <- as.factor(tmp.pa[, j])
    }
    y <- data.frame(cbind(y1, y))


    if (control$mix.method != "kproto") {
      diss <- cluster::daisy(y, metric = "gower")
      if (control$mix.method == "Gower kmeans") {
        classify <- cluster::pam(diss, k = n.g)$clustering
      }
      if (control$mix.method == "Gower hclust") {
        classify <- cutree(hclust(diss), n.g)
      }
    } else {
      classify <- kproto(y, n.g,
        iter.max = 1000,
        nstart = 100, verbose = FALSE
      )$cluster
    }
  }

  if (control$init.method == "user") {
    classify <- control$user.class
    if (length(unique(classify)) != n.g) {
      stop("Number of unique classes does not equal
           number of clusters specified in model")
    }
  }

  return(classify)
}

#' Reset defaults based on family, dimension, and expert/gating models
#'
#' @param fam distribution family
#' @param control list of options for setting up initial classification values
#' @param gate.mod true if covariates in the gating model
#' @param exp.mod true if covariates in the expert model
#' @param n.j number of response columns
#'
#' @return list of options for setting up initial classification values
#' @noRd
reset.defaults <- function(fam, control, gate.mod, exp.mod, n.j) {
  if (fam == 700) {
    if (control$init.method != "mixed") {
      if (is.element("init.method", control$defaults)) {
        control$init.method <- "mixed"
        control$mix.method <- "Gower kmeans"
      } else {
        warning("mixed init.method recommended when Tweedie family specified")
      }
    }
  }
  if (fam != 700) {
    if (n.j == 1 & is.element("init.method", control$defaults)) {
      control$init.method <- "quantile"
    }
    # default when covariate in expert or gating model and not Tweedie
    if (gate.mod | exp.mod) {
      # TODO: check and fix comment/code default when data
      # are univariate and not Tweedie and when no covariates in expert/gating
      if (n.j == 1 & is.element("init.method", control$defaults)) {
        control$init.method <- "hc"
      }
      control$hc.options$use <- "VARS"
    }
  }
  return(control)
}


#' Set initial values for mu, var, and power
#'
#' @param Class. Inital classificaiton
#' @param ysub Subset of observations based on column and cluster
#' @param expmod True if covariates in the expert model
#' @param Xd Expert covariates
#' @param family. Distribution family and link function
#'
#' @return List of initial values for mu, var, and power
#' @noRd
set.MuVarPow <- function(Class., ysub, expmod, Xd, family.) {
  out <- list(
    mu_init = 0.01,
    var_init = 0.01,
    power_init = 1.05,
    residuals = NA
  )

  if (expmod) {
    xsub <- Xd[Class. == 1, ]
    mod <- glm(ysub ~ xsub, family = family.)
    coeff <- as.vector(mod$coefficients)
    coeff[is.na(coeff)] <- 0
    out$mu_init <- coeff
    out$residuals <- mod$residuals
    out$var_init <- var(mod$residuals)
  }
  if (!expmod) {
    if (sum(ysub) != 0) {
      out$mu_init <- mean(ysub)
      out$var_init <- var(ysub)
      out$power_init <- skewness(ysub) *
        out$mu_init / sqrt(out$var_init) # Clark and Thayer, 2004
    }
  }
  return(out)
}

#' Set initial values for betad, theta, and thetaf
#'
#' @param Data. Initial Data list
#' @param inits Initial mu, var, and power
#'
#' @return List of initial values for betad, theta, and thetaf
#' @noRd
set.BetaTheta <- function(Data., inits) {
  out <- list(
    betad = inits$mu_init,
    theta = log(inits$var_init),
    thetaf = inits$power_init
  )

  if (Data.$family == 700) { # Tweedie
    ## ! ideally this will be based on link function not family
    out$betad <- log(inits$mu_init)
    if (inits$power_init >= 2) inits$power_init <- 1.95
    if (inits$power_init <= 1) inits$power_init <- 1.05
    out$thetaf <- log((1 - inits$power_init) / (inits$power_init - 2))
    ## ! adjust varaince when random effects -
    ## ParList$theta[j,g] <- log(var/mu^power.est/exp(1)) / 10
    ## exp(1) accounts for var=1 attributed to spatial
    out$theta <- log(inits$var_init / inits$mu_init^inits$power_init)
  }
  if (Data.$family == 300) {
    out$theta <- log(inits$mu_init^2 / inits$var_init)
  }

  return(out)
}

#' Set initial values for loadings parameters
#'
#' @param Data. Initial data list
#' @param cormat. Correlation matrix for observation subset
#' @param corvec. Vector of off-diaganal correlation parameters for observation subset
#' @param dimlist. List of parameter dimensions
#'
#' @return List of initial loadings parameters
#' @noRd
set.Loadings <- function(Data., cormat., corvec., dimlist.) {
  out <- list(
    logit_corr_fix = rep(0, dimlist.$nl.fix),
    ld_rand = rep(0, dimlist.$nl.rand),
    ld_sp = rep(0, dimlist.$nl.sp)
  )

  # L.mat <- t(chol(cormat))
  if (Data.$fixStruct == 30) {
    #   off.diag <- c()
    #   cnt <- 1
    #   Norm <- 1/diag(L.mat)
    #   L.mat[upper.tri(L.mat)] <- 0
    #   for(i in 2:nrow(L.mat)){
    #     for(j in 1:(i-1)){
    #       off.diag[cnt] <- L.mat[i,j]*Norm[i]
    #       cnt <- cnt+1
    #     }
    #   }
    out$logit_corr_fix <- log((corvec. + 1) / (1 - corvec.)) # off.diag
  }

  if (sum(Data.$rrStruct) > 0) {
    L.mat <- t(chol(cormat.))
    keep <- lower.tri(L.mat, diag = TRUE)
    if (Data.$rrStruct[1] == 1) {
      keep.rand <- keep
      keep.rand[, (dimlist.$n.f.rand + 1):dimlist.$n.j] <- FALSE
      out$ld_rand <- L.mat[keep.rand]
    }
    if (Data.$rrStruct[2] == 1) {
      keep.sp <- keep
      keep.sp[, (dimlist.$n.f.sp + 1):dimlist.$n.j] <- FALSE
      out$ld_sp <- L.mat[keep.sp]
    }
    if (sum(Data.$rrStruct == 2)) {
      # equal probability correlation results from
      # spatial or random rank reduction
      out$ld_rand <- out$ld_rand / 2
      out$ld_sp <- out$ld_sp / 2
    }
  }

  return(out)
}

#' Apply correction if NA in correlation matrix of observation subset
#'
#' @param cormat. Correlation matrix of observation subset
#' @param ymat. Observation subset based on cluster
#' @param nj. Number of response columns in observation matrix
#'
#' @return Corrected correlation matrix without NA values
#' @noRd
cormat.correction <- function(cormat., ymat., nj.) {
  if (sum(is.na(cormat.)) > 0) {
    idx.na <- which(is.na(cormat.), arr.ind = TRUE)
    tmp.pa <- matrix(0, nrow(ymat.), nj.)
    for (j in 1:nj.) {
      tmp.pa[ymat.[, j] > 0, j] <- 1
    }

    for (n in seq_along(idx.na[, 1])) {
      # Set up confusion matrix between 2 columns with NA correlation
      tmp.confusion <- cbind(
        c(
          nrow(
            tmp.pa[tmp.pa[, idx.na[n, 1]] == 1 &
              tmp.pa[, idx.na[n, 2]] == 1, ]
          ),
          nrow(
            tmp.pa[tmp.pa[, idx.na[n, 1]] == 0 &
              tmp.pa[, idx.na[n, 2]] == 1, ]
          )
        ),
        c(
          nrow(
            tmp.pa[tmp.pa[, idx.na[n, 1]] == 1 &
              tmp.pa[, idx.na[n, 2]] == 0, ]
          ),
          nrow(
            tmp.pa[tmp.pa[, idx.na[n, 1]] == 0 &
              tmp.pa[, idx.na[n, 2]] == 0, ]
          )
        )
      )
      # tmp.ctab[tmp.ctab==0] <- 1
      # Calculate Matthews correlation coefficient
      denom <- sum(tmp.confusion[1, ]) * sum(tmp.confusion[2, ]) *
        sum(tmp.confusion[, 1]) * sum(tmp.confusion[, 2])
      # if 0 occurs in any of the sums,
      # the denominator can arbitrarily be set to 1
      if (denom == 0) denom <- 1
      cormat.[idx.na[n, 2], idx.na[n, 1]] <-
        (tmp.confusion[1, 1] * tmp.confusion[2, 2] -
          tmp.confusion[1, 2] * tmp.confusion[2, 1]) /
          sqrt(denom)
    }
  }

  return(cormat.)
}

#' mc.qclass: quantile function from mclust. Defaults used to initiate 'E' or 'V' models when no covariates in expert/gating model
#'
#' @param x A numeric vector of observations for classification
#' @param k Integer specifying the number of mixtures
#'
#' @importFrom stats sd quantile
#' @return classification vector
#'
#' @keywords internal
#' @noRd
mc.qclass <- function(x, k) {
  x <- as.vector(x)
  eps <- sd(x) * sqrt(.Machine$double.eps)
  q <- NA
  n <- k
  while (length(q) < (k + 1)) {
    n <- n + 1
    q <- unique(quantile(x, seq(from = 0, to = 1, length = n)))
  }
  if (length(q) > (k + 1)) {
    dq <- diff(q)
    nr <- length(q) - k - 1
    q <- q[-order(dq)[1:nr]]
  }
  q[1] <- min(x) - eps
  q[length(q)] <- max(x) + eps
  cl <- rep(0, length(x))
  for (i in 1:k) {
    cl[x >= q[i] & x < q[i + 1]] <- i
  }
  return(cl)
}

#' Initialization options with S3 classes
#'
#' @param init.method Name of method used to set initial values. If init.method = 'user', must define 'user.class' with a classification vector.
#' @param hc.options Model names and use when init.method is 'hc' following conventions of mclust::mclust.options()
#' @param exp.init Turn on mahala initialization when expert network
#' @param mix.method Initialization methods when data are mixed. Default method when data are Tweedie distributed.
#' @param user.class Vector of classification vector set by user and required when init.method = 'user'
#'
#' @return list of initialization specifications
#' @export
#'
#' @examples
#' init.options()
#' init.options(init.method = "hc")
#' init.options(init.method = "mixed")
#' init.options(init.method = "user", user.class = c(1, 1, 2, 1, 3, 3, 1, 2))
init.options <- function(init.method = "hc",
                         hc.options = list(
                           modelName = "VVV",
                           use = "SVD"
                         ),
                         exp.init = list(mahala = TRUE),
                         mix.method = "Gower kmeans",
                         user.class = integer()) {
  # track defaults
  defaults <- character()
  if (missing(init.method)) defaults <- c(defaults, "init.method")
  if (missing(hc.options)) defaults <- c(defaults, "hcName", "hcUse")
  if (missing(mix.method)) defaults <- c(defaults, "mix.method")

  # match up unnamed list to correct names
  hc.options <- name.hc.options(hc.options)
  # assign hc.option defaults if one missing
  if (!("modelName" %in% names(hc.options))) {
    defaults <- c(defaults, "hcName")
  }
  if (!("use" %in% names(hc.options))) {
    defaults <- c(defaults, "hcUse")
  }

  # convert user.class to integer
  user.class <- as.integer(as.factor(user.class))

  # Set up S3 classes

  new_init <- function(method, names) {
    stopifnot(is.character(method))
    stopifnot(is.character(names))
    method <- match.arg(method, names)

    structure(method,
      class = "Init"
    )
  }

  hc_init <- function(methods = list()) {
    stopifnot(is.list(methods))
    stopifnot(length(methods) <= 2)

    methods <- list(
      modelName = match.arg(
        methods$modelName,
        c(
          "VVV", "EII", "EEE",
          "VII", "V", "E"
        )
      ),
      use = match.arg(
        methods$use,
        c(
          "SVD", "VARS", "STD", "SPH",
          "PCS", "PCR", "RND"
        )
      )
    )

    structure(methods,
      class = "HcInit"
    )
  }

  validate_user_class <- function(user,
                                  method) {
    if (length(user) == 0 & method == "user") {
      stop("user.class must be a vector of classification characters
           or integers when 'init.method = user'")
    }
    stopifnot(is.integer(user))
  }

  user_class <- function(user, method = character()) {
    validate_user_class(user, method)

    structure(user,
      class = "user"
    )
  }

  # Create S3 objects
  new_init_method <- new_init(
    method = init.method,
    names = c(
      "hc", "quantile", "random", "mclust",
      "kmeans", "mixed", "user"
    )
  )

  new_mix_method <- new_init(
    method = mix.method,
    names = c("Gower kmeans", "Gower hclust", "kproto")
  )

  new_hc_init <- hc_init(methods = hc.options)

  new_user_class <- user_class(user = user.class, method = init.method)


  out <- list(
    init.method = new_init_method[1],
    hc.options = new_hc_init[1:2],
    exp.init = exp.init,
    mix.method = new_mix_method[1],
    user.class = unclass(new_user_class),
    defaults = defaults
  )
  class(out) <- "InitOptions"

  return(out)
}

#' init.options() helper function: assigns names to hc.options list if missing
#'
#' @param methods list of hc.option methods
#'
#' @return named list of hc.option methods
#' @noRd
name.hc.options <- function(methods) {
  modelName <- c("VVV", "EII", "EEE", "VII", "V", "E")
  use <- c("SVD", "VARS", "STD", "SPH", "PCS", "PCR", "RND")

  # match up unnamed list to correct names
  nms <- names(methods)
  if (is.null(nms) & length(methods) > 0) {
    for (i in seq_along(methods)) {
      if (methods[[i]] %in% modelName) {
        names(methods)[i] <- "modelName"
      }
      if (methods[[i]] %in% use) {
        names(methods)[i] <- "use"
      }
    }
  }
  return(methods)
}

#' Updates initial class when covariates in expert formula using the Mahalanobis distance criteria
#'
#' @param z. initial classification matrix
#' @param y. response matrix
#' @param x. covariates
#' @param max.it maximum interation
#'
#' @importFrom mclust map unmap
#' @importFrom stats lm predict
#' @importFrom MoEClust MoE_mahala
#'
#' @return updated classification matrix
#' @keywords internal
#' @noRd
run.mahala <- function(z., y., x., family, max.it = 1000) {
  # modified from MoEClust
  init.exp <- TRUE
  stop.crit <- FALSE
  cnt <- 1
  M <- matrix(NA, nrow(z.), ncol(z.))
  pred <- mahala <- list()
  # df <- data.frame(y., x.)

  while (!stop.crit) {
    for (g in seq_along(z.[1, ])) {
      sub <- which(z.[, g] == 1)
      mod <- tryCatch(lm(y. ~ x., subset = sub))
      if (inherits(mod, "try-error")) {
        init.exp <- FALSE
        break
      } else {
        pred[[g]] <- cbind(rep(1, nrow(z.)), x.) %*% mod$coefficients
        res <- y. - pred[[g]]
        ## ! write my own function
        mahala[[g]] <- MoE_mahala(mod, res, squared = TRUE, identity = TRUE)
        M[, g] <- mahala[[g]]
      }
    }
    if (anyNA(M)) {
      init.exp <- FALSE
      break
    } else {
      new.z <- rep(0, nrow(z.))
      for (i in seq_along(z.[, 1])) {
        new.z[i] <- which(M[i, ] == min(M[i, ]))
      }
      if (identical(map(z.), new.z) | cnt == max.it) {
        stop.crit <- TRUE
      }
      z. <- unmap(new.z)
      cnt <- cnt + 1
    }
  }
  return(z.)
}
