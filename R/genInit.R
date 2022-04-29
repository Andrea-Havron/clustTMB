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
#' @importFrom R7 `@`
#'
#' @return list
#' @keywords internal
#' @noRd
genInit <- function(Data, family = NULL, dim.list, control = init.options()) {
  
 
  # list2env(dim.list, environment(genStart))
  n.i <- dim.list@n.i
  n.j <- dim.list@n.j
  n.t <- dim.list@n.t
  n.k.g <- ncol(Data$Xg)
  n.k.e <- ncol(Data$Xd)
  n.g <- dim.list@n.g
  n.f.rand <- dim.list@n.f.rand
  n.f.sp <- dim.list@n.f.sp
  n.v <- dim.list@n.v
  nl.fix <- dim.list@nl.fix
  nl.rand <- dim.list@nl.rand
  nl.sp <-  dim.list@nl.sp

  # Apply any data transformations
  y <- Data$Y
  if (Data$family == 300 | Data$family == 600) {
    y <- log(y)
  }
  X.g <- subset(Data$Xg, select = colnames(Data$Xg) != ("(Intercept)"))
  X.d <- subset(Data$Xd, select = colnames(Data$Xd) != ("(Intercept)"))
  gate.mod <- exp.mod <- FALSE
  if (ncol(X.g) > 0) gate.mod <- TRUE
  if (ncol(X.d) > 0) exp.mod <- TRUE

  control <- reset.defaults(Data$family, control, gate.mod, exp.mod, n.j)


  if (gate.mod) {
    y <- cbind(y, X.g)
  }
  if (exp.mod) {
    y <- cbind(y, X.d)
  }
  y <- y[, !duplicated(t(y))]

  classify <- genInitMethods(n.g, n.i, n.j,
                             control, y)
  Class <- unmap(classify)
  pi.init <- apply(Class, 2, function(x) sum(x) / n.i)
  

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
  if (Data$reStruct[1, 1] > 2) {
    if (!is.null(dim.list@n.r.g)) {
      ParList$ln_kappag <- rep(log(sqrt(8) / (dim.list@n.r.g / 2)), (n.g - 1))
    }
  }
  if (Data$reStruct[2, 1] > 2) {
    if (!is.null(dim.list@n.r.e)) {
      ParList$ln_kappad <- matrix(log(sqrt(8) / (dim.list@n.r.e / 2)), n.j, n.g)
      ParList$ln_taud <- matrix(1 / (2 * sqrt(pi) * sqrt(8) / (dim.list@n.r.e / 2)), n.j, n.g)
    }
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
  for (g in 1:n.g) {
    for (j in 1:n.j) {
      y.sub <- y[Class[, g] == 1, j]
      if (exp.mod) {
        class.sum <- apply(Class, 2, sum)
        if (any(class.sum <= 1)) stop("initalization method results in an empty or unit cluster which is not suitable when intializing the expert model")
        x.sub <- X.d[Class[, g] == 1, ]
        mod <- glm(y.sub ~ x.sub, family = family)
        coeff <- as.vector(mod$coefficients)
        coeff[is.na(coeff)] <- 0
        mu.init <- coeff
        res.mat[Class[, g] == 1, j] <- mod$residuals
        var.init <- var(mod$residuals)
      } else {
        if (sum(y.sub) == 0) {
          mu.init <- 0.01
          var.init <- 0.01
          power.est <- 1.05
        } else {
          mu.init <- mean(y.sub)
          var.init <- var(y.sub)
          power.est <- skewness(y.sub) * mu.init / sqrt(var.init) # Clark and Thayer, 2004
        }
      }

      if (Data$family == 700) { # Tweedie
        ParList$betad[, j, g] <- log(mu.init) ## ! ideally this will be based on link function not family
        if (power.est >= 2) power.est <- 1.95
        if (power.est <= 1) power.est <- 1.05
        ParList$thetaf[j, g] <- log((1 - power.est) / (power.est - 2))
        ## ! adjust varaince when random effects -  ParList$theta[j,g] <- log(var/mu^power.est/exp(1)) / 10 #exp(1) accounts for var=1 attributed to spatial
        ParList$theta[j, g] <- log(var.init / mu.init^power.est)
      } else {
        ParList$betad[, j, g] <- mu.init
        ParList$theta[j, g] <- log(var.init)
      }
      if (Data$family == 300) {
        ParList$theta[j, g] <- log(mu.init^2 / var.init)
      }
    }

    y.mat <- y[Class[, g] == 1, ]
    if (Data$fixStruct != 10) {
      if (exp.mod) {
        res <- res.mat[Class[, g] == 1, ]
        cor.mat <- cor(res)
      } else {
        cor.mat <- cor(y.mat)
      }
      # Apply correction if NA in cor.mat
      if (sum(is.na(cor.mat)) > 0) {
        idx.na <- which(is.na(cor.mat), arr.ind = TRUE)
        tmp.pa <- matrix(0, nrow(y.mat), n.j)
        for (j in 1:n.j) {
          tmp.pa[y.mat[, j] > 0, j] <- 1
        }

        for (n in 1:nrow(idx.na)) {
          # Set up confusion matrix between 2 columns with NA correlation
          tmp.confusion <- cbind(
            c(
              nrow(tmp.pa[tmp.pa[, idx.na[n, 1]] == 1 & tmp.pa[, idx.na[n, 2]] == 1, ]),
              nrow(tmp.pa[tmp.pa[, idx.na[n, 1]] == 0 & tmp.pa[, idx.na[n, 2]] == 1, ])
            ),
            c(
              nrow(tmp.pa[tmp.pa[, idx.na[n, 1]] == 1 & tmp.pa[, idx.na[n, 2]] == 0, ]),
              nrow(tmp.pa[tmp.pa[, idx.na[n, 1]] == 0 & tmp.pa[, idx.na[n, 2]] == 0, ])
            )
          )
          # tmp.ctab[tmp.ctab==0] <- 1
          # Calculate Matthews correlation coefficient
          denom <- sum(tmp.confusion[1, ]) * sum(tmp.confusion[2, ]) * sum(tmp.confusion[, 1]) * sum(tmp.confusion[, 2])
          # if 0 occurs in any of the sums, the denominator can arbitrarily be set to 1
          if (denom == 0) denom <- 1
          cor.mat[idx.na[n, 2], idx.na[n, 1]] <-
            (tmp.confusion[1, 1] * tmp.confusion[2, 2] - tmp.confusion[1, 2] * tmp.confusion[2, 1]) /
              sqrt(denom)
        }
      } # end correction
      corvec <- cor.mat[lower.tri(cor.mat)]
      # L.mat <- t(chol(cor.mat))
      if (Data$fixStruct == 30) {
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
        ParList$logit_corr_fix[, g] <- log((corvec + 1) / (1 - corvec)) # off.diag
      }

      if (sum(Data$rrStruct) > 0) {
        L.mat <- t(chol(cor.mat))
        keep <- lower.tri(L.mat, diag = TRUE)
        if (Data$rrStruct[1] == 1) {
          keep.rand <- keep
          keep.rand[, (n.f.rand + 1):n.j] <- FALSE
          ParList$ld_rand[, g] <- L.mat[keep.rand]
        }
        if (Data$rrStruct[2] == 1) {
          keep.sp <- keep
          keep.sp[, (n.f.sp + 1):n.j] <- FALSE
          ParList$ld_sp[, g] <- L.mat[keep.sp]
        }
        if (sum(Data$rrStruct == 2)) {
          # equal probability correlation results from spatial or random rank reduction
          ParList$ld_rand[, g] <- ParList$ld_rand[, g] / 2
          ParList$ld_sp[, g] <- ParList$ld_sp[, g] / 2
        }
      }
    }
  } # end g loop

  gen.init <- list(parms = ParList, class = classify)
  return(gen.init)
}



#' An R7 class to generate the model dimension list
#'
#' @importFrom R7 new_class class_double class_any prop
#' @slot n.i Number of observations
#' @slot n.j Number of response columns
#' @slot n.t Number of time steps
#' @slot n.k.g Number of covariates in gating model
#' @slot n.k.e Number of covariates in expert model
#' @slot n.g Number of clusters
#' @slot n.r.g Spatial range in gating model
#' @slot n.r.e Spatial range in expert model
#' @slot n.f.rand Number of factors in random rank reduction
#' @slot n.f.sp Number of factors in spatial rank reduction
#' @slot n.v Number of vertices in spatial INLA mesh
#' @slot nl.fix Number of loadings in covariance matrix
#' @slot nl.rand Number of factor loading in random rank reduction
#' @slot nl.sp Number of factor loading in spatial rank reduction
DimList <- new_class("DimList",
                     properties = list(
                       n.i = class_integer,
                       n.j = class_integer,
                       n.t = class_integer,
                       n.k.g = class_integer,
                       n.k.e = class_integer,
                       n.g = class_numeric,
                       n.r.g = class_any,
                       n.r.e = class_any,
                       n.f.rand = class_numeric,
                       n.f.sp = class_numeric,
                       n.v = class_numeric,
                       nl.fix = class_numeric,
                       nl.rand = class_numeric,
                       nl.sp = class_numeric
                     ),
                     validator = function(self){
                       if(round(prop(self, "n.g")) != prop(self,"n.g")){
                         "number of clusters, G, must be an integer"
                       }
                     }
)


#genInit helper functions

#' Apply classification method dependent on init.method
#'
#' @param n.g Number of clusters
#' @param n.i Number of observations
#' @param n.j Number of columns
#' @param control Classification settings from init.options()
#' @param y Observations
#'
#' @return classification vector
genInitMethods <- function(n.g, n.i, n.j,
                           control, y){
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
      classify <- kproto(y, n.g, iter.max = 1000, nstart = 100, verbose = FALSE)$cluster
    }
  }
  
  if (control$init.method == "user") {
    classify <- control$user.class
    if (length(unique(classify)) != n.g) {
      stop("Number of unique classes does not equal number of clusters specified in model")
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
#'
#' @return list of options for setting up initial classification values
reset.defaults <- function(fam, control, gate.mod, exp.mod, n.j){
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
      # TODO: check and fix comment/code default when data are univariate and not Tweedie and when no covariates in expert/gating
      if (n.j == 1 & is.element("init.method", control$defaults)) {
        control$init.method <- "hc"
      }
      control$hc.options$use <- "VARS"
    }
  }
  return(control)
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
                           use = "SVD"),
                         exp.init = list(mahala = TRUE),
                         mix.method = "Gower kmeans",
                         user.class = integer()) {
  
  #track defaults
  defaults <- character()
  if (missing(init.method)) defaults <- c(defaults, "init.method")
  if (missing(hc.options)) defaults <- c(defaults, "hcName", "hcUse")
  if (missing(mix.method)) defaults <- c(defaults, "mix.method")
  
  #match up unnamed list to correct names
  hc.options <- name.hc.options(hc.options)
  #assign hc.option defaults if one missing
  if ( !("modelName" %in% names(hc.options)) ){
    defaults <- c(defaults, "hcName")
  }
  if ( !("use" %in% names(hc.options)) ){
    defaults <- c(defaults, "hcUse")
  }
  
  #convert user.class to integer
  user.class <- as.integer(as.factor(user.class))
 
  #Set up S3 classes
  
  new_init <- function(method,names){
    stopifnot(is.character(method))
    stopifnot(is.character(names))
    method <- match.arg(method, names)
    
    structure(method,
              class = "Init")
  }
  
  hc_init <- function(methods = list()){

    stopifnot(is.list(methods))
    stopifnot(length(methods) <= 2)
    
    methods <- list(
      modelName = match.arg(methods$modelName,  
                            c( "VVV", "EII", "EEE", 
                               "VII", "V", "E")),
      use = match.arg(methods$use, 
                      c("SVD", "VARS", "STD", "SPH",
                        "PCS", "PCR", "RND"))
      )
    
    structure(methods,
              class = "HcInit")
    
  }
  
  validate_user_class <- function(user,
                                  method){
    if( length(user)==0 & method == "user" ){
      stop("user.class must be a vector of classification characters or integers when 'init.method = user'")
    }
    stopifnot(is.integer(user))
  }
  
  user_class <- function(user, method = character()){
    
    validate_user_class(user, method)
    
    structure(user,
              class = "user")
  }
  
  #Create S3 objects
  new_init_method <- new_init(method = init.method, 
                          names = c("hc", "quantile", "random", "mclust", 
                                    "kmeans", "mixed", "user"))

  new_mix_method <- new_init(method = mix.method,
                         names = c("Gower kmeans", "Gower hclust", "kproto"))
  
  new_hc_init = hc_init(methods = hc.options)

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


#' Initialization options with S4 classes
#'
#' @param init.method Name of method used to set initial values. If init.method = 'user', must define 'user.class' with a classification vector.
#' @param hc.options Model names and use when init.method is 'hc' following conventions of mclust::mclust.options()
#' @param exp.init Turn on mahala initialization when expert network
#' @param mix.method Initialization methods when data are mixed. Default method when data are Tweedie distributed.
#' @param user.class Vector of classification vector set by user and required when init.method = 'user'
#' 
#'
#' @return list of initialization specifications
#' @export
#'
#' @examples
#' init.options()
#' init.options(init.method = "hc")
#' init.options(init.method = "mixed")
#' init.options(init.method = "user", user.class = c(1, 1, 2, 1, 3, 3, 1, 2))
init.options.s4 <- function(init.method = "hc",
                         hc.options = list(
                           modelName = "VVV",
                           use = "SVD"),
                         exp.init = list(mahala = TRUE),
                         mix.method = "Gower kmeans",
                         user.class = integer()) {
  
  
  #track defaults
  defaults <- character()
  if (missing(init.method)) defaults <- c(defaults, "init.method")
  if (missing(hc.options)) defaults <- c(defaults, "hcName", "hcUse")
  if (missing(mix.method)) defaults <- c(defaults, "mix.method")
  
  #match up unnamed list to correct names
  hc.options <- name.hc.options(hc.options)
  #assign hc.option defaults if one missing
  if ( !("modelName" %in% names(hc.options)) ){
    defaults <- c(defaults, "hcName")
  }
  if ( !("use" %in% names(hc.options)) ){
    defaults <- c(defaults, "hcUse")
  }
  
  #convert user.class to integer
  user.class <- as.integer(as.factor(user.class))
  
  #set Class objects
  Init <- setClass(
    "Init",
    
    slots = c(method = "character", 
              names = "character"),
    
  )

  HcInit <- setClass(
    "HcInit",
    
    slots = c(methods = "list"),

  )
  setValidity("HcInit", function(object){
    if(length(object@methods) > 2){
      stop("hc.options must contain a single argument each for modelNames and use")
    }
  })
  
  UserInit <- setClass(
    "UserInit",
    
    slots = c(user = "integer",
              method = "character")
  )
  
  setValidity("UserInit", function(object){
    if( length(object@user)==0 & object@method == "user" ){
      stop("user.class must be a vector of classification characters or integers when 'init.method = user'")
    } else {
      TRUE
    }
  })
  
  # set generic getInit Class and methods
  setGeneric( name = "getInit",
              def = function(object){
                standardGeneric("getInit")
              }
              ) 
  
  setMethod(f = "getInit",
            signature = "Init",
            definition = function(object){
              return(match.arg(object@method, object@names))
            }
            )

  setMethod(f = "getInit",
            signature = "HcInit",
            definition = function(object){
            
              methods <- list(
                modelName = match.arg(object@methods$modelName, 
                                      c( "VVV", "EII", "EEE",  
                                         "VII", "V", "E")),
                use = match.arg(object@methods$use,
                                c("SVD", "VARS", "STD", "SPH", 
                                  "PCS", "PCR", "RND"))
              )
              
              out <- methods
              return(out)
            })
  
  
  #Create new objects and apply method
  init_method <- Init(method = init.method, 
                      names = c("hc", "quantile", "random", "mclust", 
                                        "kmeans", "mixed", "user"))
  init.method <- getInit(init_method)


  hc_options <- HcInit(methods = hc.options)
  hc.options <-  getInit(hc_options)
  
  
  mix_method <- Init(method = mix.method, 
                        names =  c("Gower kmeans", 
                                "Gower hclust", "kproto"))
  mix.method <- getInit(mix_method)
  

  user_class <- UserInit(user = user.class, 
                          method = init.method)

  out <- list(
    init.method = init.method,
    hc.options = hc.options,
    exp.init = exp.init,
    mix.method = mix.method,
    user.class = user_class@user,
    defaults = defaults
  )
  class(out) <- "InitOptions"
  
  return(out)
}

#' Initialization options with R7 classes
#'
#' @param init.method Name of method used to set initial values. If init.method = 'user', must define 'user.class' with a classification vector.
#' @param hc.options Model names and use when init.method is 'hc' following conventions of mclust::mclust.options()
#' @param exp.init Turn on mahala initialization when expert network
#' @param mix.method Initialization methods when data are mixed. Default method when data are Tweedie distributed.
#' @param user.class Vector of classification vector set by user and required when init.method = 'user'
#' @importFrom R7 new_class new_generic R7_dispatch method class_any class_character class_integer class_numeric
#'
#' @return list of initialization specifications
#' @export
#'
#' @examples
#' init.options()
#' init.options(init.method = "hc")
#' init.options(init.method = "mixed")
#' init.options(init.method = "user", user.class = c(1, 1, 2, 1, 3, 3, 1, 2))
init.options.r7 <-  function(init.method = "hc",
                             hc.options = list(
                               modelName = "VVV",
                               use = "SVD"),
                             exp.init = list(mahala = TRUE),
                             mix.method = "Gower kmeans",
                             user.class = integer()) {
  
  #track defaults
  defaults <- character()
  if (missing(init.method)) defaults <- c(defaults, "init.method")
  if (missing(hc.options)) defaults <- c(defaults, "hcName", "hcUse")
  if (missing(mix.method)) defaults <- c(defaults, "mix.method")
  
  #match up unnamed list to correct names
  hc.options <- name.hc.options(hc.options)
  #assign hc.option defaults if one missing
  if ( !("modelName" %in% names(hc.options)) ){
    defaults <- c(defaults, "hcName")
  }
  if ( !("use" %in% names(hc.options)) ){
    defaults <- c(defaults, "hcUse")
  }
  
  #convert user.class to integer
  user.class <- as.integer(as.factor(user.class))
  
  #Set up classes
  Init <- new_class("Init",
                    properties = list(
                      method = class_any, 
                      namevector = class_character
                    ),
                    validator = function(self){
                      if(!is.character(self@method)){
                        "method from init.options() must be a character"
                      }
                    })
  HcInit <- new_class("HcInit",
                     properties = list(
                      methods = class_any
                     ),
                     validator = function(self){
                       if(!is.list(self@methods)){
                         "hc.options from init.options() must be a list"
                       }
                       if(length(self@methods) > 2){
                         "hc.options from init.options() must contain a single argument each for modelNames and use"
                       }
                     })
  UserInit <- new_class("UserInit",
                        properties = list(
                          user = class_integer,
                          method = class_character
                        ),
                        validator = function(self){
                          if( length(self@user)==0 & self@method == "user" ){
                            "user.class must be a vector of classification factors, characters or integers when 'init.method = user'"
                          }
                        })
  
  #Set up generic function
  getInit <- new_generic("getInit", "object", function(object){
    R7_dispatch()
  })
  
  #Set methods for generic function
  method(getInit, Init) <- function(object){
    match.arg(object@method, object@namevector)
  }
  
  method(getInit, HcInit) <- function(object){
    list(
      modelName = match.arg(object@methods$modelName, 
                            c( "VVV", "EII", "EEE",  
                               "VII", "V", "E")),
      use = match.arg(object@methods$use,
                      c("SVD", "VARS", "STD", "SPH", 
                        "PCS", "PCR", "RND"))
    )
  }
  
  #Create new objects and apply method
  init.method <- getInit(
    Init(
      method = init.method,
      namevector = c("hc", "quantile", "random", "mclust", 
        "kmeans", "mixed", "user")
      )
  )
  
  hc.options <- getInit(
    HcInit(
      methods = hc.options
    )
  )
  
  mix.method <- getInit(
    Init(
      method = mix.method,
      namevector = c("Gower kmeans", 
                     "Gower hclust", "kproto")
    )
  )
  
  user.class <- UserInit(user = user.class, method = init.method)
  
  out <- list(
    init.method = init.method,
    hc.options = hc.options,
    exp.init = exp.init,
    mix.method = mix.method,
    user.class = user.class@user,
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
name.hc.options <- function(methods){
  modelName <- c( "VVV", "EII", "EEE",  "VII", "V", "E")
  use <- c("SVD", "VARS", "STD", "SPH",  "PCS", "PCR", "RND")
  
  #match up unnamed list to correct names
  nms <- names(methods)
  if(is.null(nms) & length(methods) > 0){
    for(i in 1:length(methods)){
      if(methods[[i]] %in% modelName){
        names(methods)[i] <- "modelName"
      }
      if(methods[[i]] %in% use){
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
#'
run.mahala <- function(z., y., x., family, max.it = 1000) {
  # modified from MoEClust
  init.exp <- TRUE
  stop.crit <- FALSE
  cnt <- 1
  M <- matrix(NA, nrow(z.), ncol(z.))
  pred <- mahala <- list()
  # df <- data.frame(y., x.)

  while (!stop.crit) {
    for (g in 1:ncol(z.)) {
      sub <- which(z.[, g] == 1)
      mod <- tryCatch(lm(y. ~ x., subset = sub))
      if (inherits(mod, "try-error")) {
        init.exp <- FALSE
        break
      } else {
        pred[[g]] <- cbind(rep(1, nrow(z.)), x.) %*% mod$coefficients
        res <- y. - pred[[g]]
        mahala[[g]] <- MoE_mahala(mod, res, squared = TRUE, identity = TRUE) ## ! write my own function
        M[, g] <- mahala[[g]]
      }
    }
    if (anyNA(M)) {
      init.exp <- FALSE
      break
    } else {
      new.z <- rep(0, nrow(z.))
      for (i in 1:nrow(z.)) {
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
