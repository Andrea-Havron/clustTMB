#' genInit: generates initial values for cluster model
#'
#' @param init.method Method as string used to generate initial values. Options include: "random", "hc"
#' @param Data List containing TMB data objects
#' @param family Distribution family
#' @param dim.list List of model dimensions
#' @param control Calls init.options() to generate settings for initial values. Arguments of init.options() can be specified by the user.
#' 1. init.method - Single character string indicating intial clustering method. Mehtods include: hc, quantile, random, mclust, kmeans, mixed, user. Defaults to 'hc'. In the case where data are univariate and there are no covariates in the gating/expert formula, this defaults to 'quantile'
#' 2. hc.options - Named list of two character strings specifying hc modelName and hcUse when init.method = 'hc'. The default modelName is 'VVV' and the default use is 'SVD' unless gating/expert covariates specified, in which case default in VARS. See ?mclust::mclust.options for complete list of options.
#' 3. mix.method - String stating initialization method for mixed-type data (init.method = 'mixed'). Current default when Tweedie family specified. Options include: Gower kmeans (default), Gower hclust, and kproto.
#' 4. user - Numeric or character vector defining user specified intial classification. init.method must be set to 'user' when using this option.
#'
#' @importFrom stats cutree gaussian hclust runif rmultinom nlminb cor glm
#' @importFrom mclust unmap hclass hc hcVVV hcE hcEII hcEEE hcVII hcV
#' @importFrom e1071 skewness
#' @importFrom cluster daisy pam
#' @importFrom geoR boxcoxfit
#' @importFrom clustMixType kproto
#' @importFrom nnet multinom
#'
#' @return list
#' @keywords internal
#' @noRd
#'
#' @examples
genInit <- function(Data, family = NULL, dim.list, control = init.options()){

  #list2env(dim.list, environment(genStart))
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
  n.r <- dim.list$n.r
  nl.fix <- ifelse(n.j>1, (n.j^2-n.j)/2, 1)
  nl.rand <- n.j*n.f.rand - (n.f.rand*(n.f.rand-1))/2
  nl.sp <- n.j*n.f.sp - (n.f.sp*(n.f.sp-1))/2


#  if(n.j == 1) hc.options[1] <- 'E'

  # ## Apply any data transformations
  y <- Data$Y
  # if(!is.na(data.trans)){
  #   if(data.trans == "log+1"){
  #     for(j in 1:n.j){
  #       u <- runif(n.i, 0.1,1)
  #       y[,j] <- log(Data$Y[,j]+u)
  #     }
  #   }
  #   if(data.trans == "log"){
  #     y <- log(Data$Y)
  #   }
  #   if(data.trans == "Yeojin Boxcox"){
  #     for(j in 1:n.j){
  #       u <- runif(n.i, 0.1,1)
  #       boxlambda <- geoR::boxcoxfit(Data$Y[,j] + u, add.to.data = 0)$lambda
  #       y[,j] <- ((Data$Y[,j] + u)^boxlambda - 1)/boxlambda
  #     }
  #   }
  # }

  X.g <- subset(Data$Xg, select = colnames(Data$Xg) != ('(Intercept)'))
  X.d <- subset(Data$Xd, select = colnames(Data$Xd) != ('(Intercept)'))
  gate.mod <- exp.mod <- FALSE
  if(ncol(X.g)>0) gate.mod <- TRUE
  if(ncol(X.d)>0) exp.mod <- TRUE

  ## reset defaults based on family, dimension, and expert/gating models
  if(Data$family == 700){
    if(control$init.method != 'mixed'){
      if(is.element('init.method', control$defaults)){
        control$init.method <- 'mixed'
        control$mix.method <- 'Gower kmeans'
      } else{
        warning('mixed init.method recommended when Tweedie family specified')
      }
    }
  } else {
    if(dim.list$n.j == 1 & is.element('init.method', control$defaults)){
      control$init.method <- 'quantile'
    }
    # default when covariate in expert or gating model and not Tweedie
    if( gate.mod | exp.mod  ){
      # default when data are univariate and not Tweedie and when no covariates in expert/gating
      if(dim.list$n.j == 1 & is.element('init.method', control$defaults)){
        control$init.method <- 'hc'
      }
      control$hc.options$use <- 'VARS'
    }
  }

  if(gate.mod){
    y <- cbind(y, X.g)
  }
  if(exp.mod){
    y <- cbind(y, X.d)
  }
  y <- y[,!duplicated(t(y))]

  #Apply classification method
  if(control$init.method =="random"){
    pi.init <- rep(1/n.g,n.g)
    Class <- t(rmultinom(n.i, 1, pi.init)) ## FIXME: write my own rmultinom function to avoid dependency here
  }

  if(control$init.method == 'quantile'){
    classify <- mc.qclass(y, as.numeric(n.g))
    Class <- matrix(0, n.i, n.g)
    for(i in 1:n.i){
      Class[i,classify[i]] <- 1
    }
    pi.init <- apply(Class,2,sum)/n.i
  }

  if(control$init.method == "hc"){
    classify <- as.vector(hclass(
      hc(y, modelName = control$hc.options$modelName,
         use = control$hc.options$use), n.g
      ))
    Class <- unmap(classify)
    pi.init <-  apply(Class,2,function(x) sum(x)/nrow(Data$Y))

  }

  if(control$init.method == "kmeans"){
    classify <- cluster::pam(diss, k = n.g)$clustering
    Class <- unmap(classify)
    pi.init <- apply(Class,2,function(x) sum(x)/nrow(Data$Y))
  }

  if(control$init.method == "mixed"){
    tmp.pa <- matrix(NA, n.i, n.j)
    y1 <- as.data.frame(matrix(NA, n.i, n.j) )
    for(j in 1:n.j){
      tmp.pa[,j] <- ifelse(y[,j]==0, "A", "P")
      y1[,j] <- as.factor(tmp.pa[,j])
    }
    y <- data.frame(cbind(y1, y))


    if(control$mix.method != "kproto"){
      diss <- cluster::daisy(y, metric = "gower")
      if(control$mix.method == "Gower kmeans"){
        classify <- cluster::pam(diss, k = n.g)$clustering
      }
      if(control$mix.method == "Gower hclust"){
        classify <- cutree(hclust(diss),n.g)
      }
    } else {
      classify <- kproto(y, n.g, iter.max = 1000, nstart = 100, verbose = FALSE)$cluster
    }
    Class <- unmap(classify)
    pi.init <- apply(Class,2,function(x) sum(x)/nrow(Data$Y))
  }

  if(control$init.method == "user"){
    classify <- control$user.class
    if(length(unique(classify)) != n.g){
      stop("Number of unique classes does not equal number of clusters specified in model")
    }

    Class <- unmap(classify)
    pi.init <- apply(Class,2,function(x) sum(x)/nrow(Data$Y))
  }


  #setup ParList
  ParList <- list(
    betag = matrix(0,n.k.g,(n.g-1)),
    betad = array(0, dim = c(n.k.e,n.j,n.g)),
    betapz = array(numeric(0)),
    theta = matrix(0,n.j,n.g),
    thetaf = matrix(0,n.j,n.g),
    logit_corr_fix = matrix(0,nl.fix,n.g),
    ld_rand = matrix(0,nl.rand,n.g),
    ld_sp = matrix(0,nl.sp,n.g),
    Hg_input = matrix(0,2,(n.g-1)),
    Hd_input = array(0, dim = c(2,n.j,n.g)), #fix in cpp - should be n.f, not n.j
    ln_kappag = rep(0, (n.g-1)),
    ln_kappad = matrix(0,n.j,n.g), #fix in cpp - should be n.f, not n.j
    ln_taud = matrix(0,n.j,n.g), #fix in cpp - should be n.f, not n.j
    logit_rhog = rep(0,(n.g-1)),
    logit_rhod = matrix(0,n.j,n.g),
    ln_sigmaup = rep(0,(n.g-1)),
    ln_sigmaep = matrix(0,n.j,n.g),
    ln_sigmau = rep(0,(n.g-1)),
    ln_sigmav = matrix(0,n.j,n.g),
    upsilon_tg = array(0, dim = c(n.t,(n.g-1))),
    epsilon_tjg = array(0, dim = c(n.t,n.j,n.g)),
    u_ig = array(0, dim = c(n.i,(n.g-1))),
    v_ifg = array(0, dim = c(n.i,n.j,n.g)),
    Gamma_vg = array(0, dim = c(n.v,(n.g-1))),
    Omega_vfg = array(0, dim = c(n.v,n.f.sp,n.g))
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

  #Set initial values for kappa and tau if n.r provided
  if(Data$reStruct[1,1] > 2){
    if(!is.null(dim.list$n.r)){
      ParList$ln_kappag = rep(log(sqrt(8)/(n.r/2)), (n.g-1))
    }
  }
  if(Data$reStruct[2,1] > 2){
    if(!is.null(dim.list$n.r)){
      ParList$ln_kappad = matrix(log(sqrt(8)/(n.r/2)), n.j,n.g)
      ParList$ln_taud = matrix( 1/(2*sqrt(pi)*sqrt(8)/(n.r/2)) ,n.j,n.g)
    }
  }


  #Update initial values based on classification
  res.mat <- matrix(0, n.i, n.j)
  if(exp.mod & control$exp.init$mahala){
    Class <- run.mahala(Class, as.matrix(Data$Y), as.matrix(X.d))
  }
  for(g in 1:n.g){
    for(j in 1:n.j){
      y.sub <- Data$Y[Class[,g] == 1,j]
      if(exp.mod){
        class.sum <- apply(Class,2,sum)
        if(any(class.sum <= 1)) stop("initalization method results in an empty or unit cluster which is not suitable when intializing the expert model")
        x.sub <- X.d[Class[,g] == 1,]
        mod <- glm(y.sub ~ x.sub, family = family)
        coeff <- as.vector(mod$coefficients)
        coeff[is.na(coeff)] <- 0
        mu.init <- coeff
        res.mat[Class[,g] == 1,j] <- mod$residuals
        var.init <- var(mod$residuals)
      } else {
        if(sum(y.sub) == 0){
          mu.init <- 0.01
          var.init <- 0.01
          power.est <- 1.05
        } else {
          mu.init <- mean(y.sub)
          var.init <- var(y.sub)
          power.est <-  e1071::skewness(y.sub)*mu.init/sqrt(var.init)  # Clark and Thayer, 2004
        }
      }

      if(Data$family == 700){ #Tweedie
        ParList$betad[,j,g] <- log(mu.init) ##! ideally this will be based on link function not family
        if(power.est >= 2) power.est <- 1.95
        if(power.est <= 1) power.est <- 1.05
        ParList$thetaf[j,g] <- log((1-power.est)/(power.est - 2))
        ##! adjust varaince when random effects -  ParList$theta[j,g] <- log(var/mu^power.est/exp(1)) / 10 #exp(1) accounts for var=1 attributed to spatial
        ParList$theta[j,g] <- log(var.init/mu.init^power.est)
      } else {
        ParList$betad[,j,g] <- mu.init
        ParList$theta[j,g] <- log(var.init)
      }
    }

    y.mat <- Data$Y[Class[,g] == 1,]
    if(Data$fixStruct != 10){
      if(exp.mod){
        res <- res.mat[Class[,g] == 1,]
        cor.mat <- cor(res)
      } else {
        cor.mat <- cor(y.mat)
      }
      if(sum(is.na(cor.mat))>0){
        idx.na <- which(is.na(cor.mat), arr.ind = TRUE)
        tmp.pa <- matrix(0, nrow(y.mat), n.j)
        for(j in 1:n.j){
          tmp.pa[y.mat[,j]>0,j] <- 1
        }

        for(n in 1:nrow(idx.na)){
          #Set up confusion matrix between 2 columns with NA correlation
          tmp.confusion <- cbind(
            c(nrow(tmp.pa[tmp.pa[,idx.na[n,1]]==1 & tmp.pa[,idx.na[n,2]]==1,]),
              nrow(tmp.pa[tmp.pa[,idx.na[n,1]]==0 & tmp.pa[,idx.na[n,2]]==1,])),
            c(nrow(tmp.pa[tmp.pa[,idx.na[n,1]]==1 & tmp.pa[,idx.na[n,2]]==0,]),
              nrow(tmp.pa[tmp.pa[,idx.na[n,1]]==0 & tmp.pa[,idx.na[n,2]]==0,]))
          )
          #tmp.ctab[tmp.ctab==0] <- 1
          #Calculate Matthews correlation coefficient
          denom = sum(tmp.confusion[1,])*sum(tmp.confusion[2,])*sum(tmp.confusion[,1])*sum(tmp.confusion[,2])
          #if 0 occurs in any of the sums, the denominator can arbitrarily be set to 1
          if(denom == 0) denom <- 1
          cor.mat[idx.na[n,2], idx.na[n,1]] <-
            (tmp.confusion[1,1]*tmp.confusion[2,2] - tmp.confusion[1,2]*tmp.confusion[2,1]) /
            sqrt(denom)
        }
      }
      corvec <- cor.mat[lower.tri(cor.mat)]
      if(Data$fixStruct == 30){
        ParList$logit_corr_fix[,g] <- log((corvec+1)/(1-corvec))
      }

      if(sum(Data$rrStruct)>0){
        L.mat <- t(chol(cor.mat))
        keep <- lower.tri(L.mat, diag = TRUE)
        if(Data$rrStruct[1] == 1){
          keep.rand <- keep
          keep.rand[,(n.f.rand+1):n.j] <- FALSE
          ParList$ld_rand[,g] = L.mat[keep.rand]
        }
        if(Data$rrStruct[2] == 1){
          keep.sp <- keep
          keep.sp[,(n.f.sp+1):n.j] <- FALSE
          ParList$ld_sp[,g] = L.mat[keep.sp]
        }
        if(sum(Data$rrStruct == 2)){
          #equal probability correlation results from spatial or random rank reduction
          ParList$ld_rand[,g] <- ParList$ld_rand[,g]/2
          ParList$ld_sp[,g] <- ParList$ld_sp[,g]/2
        }
      }
    }
  } # end g loop

  gen.init <- list(parms = ParList, class = classify)
  return(gen.init)

}

#' mc.qclass: quantile function from mclust. Defaults used to initate 'E' or 'V' models when no covariates in expert/gating model
#'
#' @param x A numeric vector of observations for classification
#' @param k Integer specifying the number of mixtures
#'
#' @return classification vector
#'
#' @keywords internal
#' @noRd
mc.qclass <- function (x, k)
{
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

#' Initialization options
#'
#' @param init.method
#' @param hc.options
#' @param mix.method
#' @param user.class
#' @param defaults
#'
#' @return
#' @export
#'
#' @examples
init.options <- function(init.method = c("hc", "quantile", "random", "mclust", "kmeans", "mixed", "user"),
                         hc.options = list(
                           modelName = c("VVV", "EII", "EEE", "VII", "V", "E"),
                           use = c("SVD", "VARS", "STD", "SPH", "PCS", "PCR",  "RND")),
                         exp.init = list(mahala = TRUE),
                         mix.method = c("Gower kmeans", "Gower hclust", "kproto"),
                         user.class = NULL, defaults = c()){
  if(missing(init.method)) defaults <- c(defaults, "init.method")
  if(!missing(init.method) & length(init.method) > 1 | !is.character(init.method) ){
    stop(" 'init.method' must be a single character string, see ?init.options() for valid init.method")
  }
  if(!missing(init.method)){
    if(!is.element(init.method, c("quantile", "hc", "random", "mclust", "kmeans", "mixed", "user"))){
      stop(" specified init.method not supported, see ?init.options() for valid init.method")
    }
  }
  init.method <- match.arg(init.method)

  if(init.method == 'hc'){
    if(length(hc.options) > 2 | !is.list(hc.options)){
      stop("hc.options must be list of two items: modelName and use, see ?init.options() for details")
    }
    if(missing(hc.options)){
      hc.options <- list(modelName = "VVV", use = "SVD")
      defaults <-c(defaults, 'hcName', 'hcUse')
    }
    nm <- names(hc.options)
    if(!is.null(nm) & (!(all(is.element(nm, c("modelName", "use")))))){
      stop("hc.options must be a named list of two items: modelName and use, see ?init.options() for details")
    }
    if(is.null(nm) & length(hc.options) == 2){
      names(hc.options) <- c("modelName", "use")
    }
    if(is.null(nm) & length(hc.options) == 1){
      stop("ambiguous specification of hc.options, see ?init.options() for details")
    }
    if(length(nm) == 1){
      if(nm == 'modelName'){
        hc.options$use <- "SVD"
        defaults <- c(defaults, 'hcUse')
      }
      if(nm == "use"){
        hc.options$modelName <- "VVV"
        defaults <- c(defaults, 'hcName')
      }
    }
    if(!is.element(hc.options$modelName, c("VVV", "EII", "EEE", "VII"))){
      stop(" hc.options modelName not supported, see ?init.options()")
    }
    if(!is.element(hc.options$use, c("VARS", "STD", "SPH", "PCS", "PCR", "SVD", "RND"))){
      stop(" hc.options use not supported, see ?init.options()")
    }

  }
  if(init.method == "mixed"){
    if(!missing(mix.method) & length(mix.method) > 1 | !is.character(mix.method)){
      stop("mix.method needs to be a single character string, see ?init.options() for details")
    }
    if(missing(mix.method)){
      mix.method <- match.arg(mix.method)
      defaults <- c(defaults, 'mix.method')
    }
    if(!is.element(mix.method, c("Gower kmeans", "Gower hclust", "kproto"))){
      stop("mix.method not supported, see ?init.options() for details")
    }
  }
  if(init.method == "user"){
    if(is.null(user.class)){
      stop("user.class must be a vector of classification characters or integers when 'init.method = user'")
    }
    if(!is.vector(user.class)){
      stop("user.class must be a vector of classification characters or integers when 'init.method = user'")
    }
    if(is.character(user.class)){
      user.class <- as.numeric(as.factor(user.class))
    }
  }

  return(list(init.method = init.method,
              hc.options = hc.options,
              exp.init = exp.init,
              mix.method = mix.method,
              user.class = user.class,
              defaults = defaults))


}

#' Updates initial class when covariates in expert formula using the Mahalanobis distance criteria
#'
#' @param z. initial classification matrix
#' @param y. response matrix
#' @param x. covariates
#' @param max.it maximum interation
#'
#' @importFrom mclust map unmap
#' @importFrom stats predict
#' @importFrom MoEClust MoE_mahala
#'
#' @return updated classification matrix
#' @keywords internal
#'
run.mahala <- function(z.,y.,x.,family,max.it = 1000){
  #modified from MoEClust
  init.exp <- TRUE
  stop.crit <- FALSE
  cnt <- 1
  M <- matrix(NA,nrow(z.),ncol(z.))
  pred <- mahala <- list()
 # df <- data.frame(y., x.)

  while(!stop.crit){
    for(g in 1:ncol(z.)){
      sub <- which(z.[,g]==1)
      mod <- tryCatch(lm(y.~x., subset = sub))
      if(inherits(mod, "try-error")){
        init.exp <- FALSE
        break
      } else {
        pred[[g]] <- cbind(rep(1,nrow(z.)),x.) %*% mod$coefficients
        res <- y. - pred[[g]]
        mahala[[g]] <- MoE_mahala(mod, res, squared=TRUE, identity=TRUE) ##! write my own function
        M[,g] <- mahala[[g]]
      }
    }
    if(anyNA(M)){
      init.exp <- FALSE
      break
    } else {
      new.z <- rep(0,nrow(z.))
      for(i in 1:nrow(z.)){
        new.z[i] <- which(M[i,]==min(M[i,]))
      }
      if(identical(map(z.), new.z) | cnt == max.it){ stop.crit <- TRUE}
      z. <- unmap(new.z)
      cnt <- cnt + 1
    }
  }
  return(z.)
}


