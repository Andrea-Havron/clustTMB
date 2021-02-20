#' genInit: generates initial values for cluster model
#'
#' @param init.method Method as string used to generate initial values. Options include: "random", "hc"
#' @param Data List containing TMB data objects
#' @param family Distribution family
#' @param dim Integer
#' @param data.trans String identifying any data transformation
#' @param hc.options Hierarchical clustering settings. Position[1] is model name -see mclust.options("hcModelName"). Position[2] is data transformation -see mclust.options("hcUse").
#' @param mix.method String stating initialization method
#' @param r.tol Relative convergence tolerance
#' @param true.class Vector of true classification integers for testing models
#'
#' @importFrom stats cutree gaussian hclust runif rmultinom nlminb cor
#' @importFrom mclust unmap hclass hc hcVVV hcE
#' @importFrom e1071 skewness
#' @importFrom cluster daisy pam
#' @importFrom geoR boxcoxfit
#' @importFrom clustMixType kproto
#'
#' @return list
#' @keywords internal
#' @noRd
#'
#' @examples
genInit <- function(init.method = NULL, Data, family = NULL, dim.list, data.trans = NA, hc.options = c("VVV", "SVD"), #default matches mclust
                     mix.method = NULL, r.tol = 1e-10, true.class = NULL){

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


  if(n.j == 1) hc.options[1] <- 'E'

  ## Apply any data transformations
  y <- Data$Y
  if(!is.na(data.trans)){
    if(data.trans == "log+1"){
      for(j in 1:n.j){
        u <- runif(n.i, 0.1,1)
        y[,j] <- log(Data$Y[,j]+u)
      }
    }
    if(data.trans == "log"){
      y <- log(Data$Y)
    }
    if(data.trans == "Yeojin Boxcox"){
      for(j in 1:n.j){
        u <- runif(n.i, 0.1,1)
        boxlambda <- geoR::boxcoxfit(Data$Y[,j] + u, add.to.data = 0)$lambda
        y[,j] <- ((Data$Y[,j] + u)^boxlambda - 1)/boxlambda
      }
    }
  }

  #set default init.methods based on dim(y)
  if(is.null(init.method)){
    if(n.j == 1){
      init.method <- 'mc.qclass'
    }
    if(n.j >1){
      if(Data$family == 0){
        init.method <- 'hc'
        hc.options = c("VVV", "SVD")
      }
      if(Data$family == 700){
        init.method <- 'mixed'
        mix.method <- 'Gower kmeans'
      }
    }
  }


  X <- subset(Data$Xg, select = colnames(Data$Xg) != ('(Intercept)'))
  if(ncol(X)>0){
    y <- cbind(y, X)
    if(init.method == 'mc.qclass'){
      init.method <- 'hc'
    }
    hc.options = c("VVV", "VARS")
  }

    #Apply classification method
    if(init.method =="random"){
      pi.init <- rep(1/n.g,n.g)
      Class <- t(rmultinom(n.i, 1, pi.init)) ## FIXME: write my own rmultinom function to avoid dependency here
    }

    if(init.method == 'mc.qclass'){
      classify <- mc.qclass(y, as.numeric(n.g))
      Class <- matrix(0, n.i, n.g)
      for(i in 1:n.i){
        Class[i,classify[i]] <- 1
      }
      pi.init <- apply(Class,2,sum)/n.i
    }

    if(init.method == "hc"){
      #Class <- unmap(hclass(hc(y),n.g))
      Class <- unmap(hclass(hc(y, modelName = hc.options[1], use = hc.options[2]),n.g))
      pi.init <-  apply(Class,2,function(x) sum(x)/nrow(Data$Y))
      classify <- apply(Class,1,function(x) which(x==1))

    }

    if(init.method == "kmeans"){
      classify <- cluster::pam(diss, k = n.g)$clustering
      Class <- matrix(0, n.i, n.g)
      for(i in 1:n.i){
        Class[i,classify[i]] <- 1
      }
      pi.init <- apply(Class,2,function(x) sum(x)/nrow(Data$Y))
    }

    if(init.method == "mixed"){
      tmp.pa <- matrix(NA, n.i, n.j)
      y1 <- as.data.frame(matrix(NA, n.i, n.j) )
      for(j in 1:n.j){
        tmp.pa[,j] <- ifelse(y[,j]==0, "A", "P")
        y1[,j] <- as.factor(tmp.pa[,j])
      }
      y <- data.frame(cbind(y1, y))


      if(mix.method != "kproto"){
        diss <- cluster::daisy(y, metric = "gower")
        if(mix.method == "Gower kmeans"){
          classify <- cluster::pam(diss, k = n.g)$clustering
        }
        if(mix.method == "Gower hclust"){
          classify <- cutree(hclust(diss),n.g)
        }
      } else {
        classify <- kproto(y, n.g, iter.max = 1000, nstart = 100, verbose = FALSE)$cluster
      }
      Class <- matrix(0, n.i, n.g)
      for(i in 1:n.i){
        Class[i,classify[i]] <- 1
      }
      pi.init <- apply(Class,2,function(x) sum(x)/nrow(Data$Y))
    }

    if(init.method == "true"){
      if(!is.null(true.class)){
        classify <- true.class
      } else {
        classify <- Data$classification + 1
      }
      Class <- matrix(0, n.i, n.g)
      for(i in 1:n.i){
        Class[i,classify[i]] <- 1
      }
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
  if(sum(Data$reStruct[1,])==0){
    ParList$betag[1,] <- matrix(log(pi.init[1:(n.g-1)]/(1 - sum(pi.init[1:(n.g-1)]))), nrow = 1)
  }

  #Set initial values for kappa and tau if n.r provided
  if(Data$reStruct[1,1] > 2){
    if(is.null(dim.list$n.r)){
      ParList$ln_kappag = rep(log(sqrt(8)/(n.r/2)), (n.g-1))
    }
  }
  if(Data$reStruct[2,1] > 2){
    if(is.null(dim.list$n.r)){
      ParList$ln_kappad = matrix(log(sqrt(8)/(n.r/2)), n.j,n.g)
      ParList$ln_taud = matrix( 1/(2*sqrt(pi)*sqrt(8)/(n.r/2)) ,n.j,n.g)
    }
  }


  #Update initial values based on classification
  for(g in 1:n.g){
    for(j in 1:n.j){
      y.sub <- Data$Y[Class[,g] == 1,j]
      if(sum(y.sub) == 0){
        mu <- 0.01
        var <- 0.01
        power.est <- 1.05
      } else {
        mu <- mean(y.sub)
        var <- var(y.sub)
        power.est <-  e1071::skewness(y.sub)*mu/sqrt(var)  # Clark and Thayer, 2004
      }

      if(Data$family == 700){ #Tweedie
        ParList$betad[1,j,g] <- log(mu) ##! ideally this will be based on link function not family
        if(power.est >= 2) power.est <- 1.95
        if(power.est <= 1) power.est <- 1.05
        ParList$thetaf[j,g] <- log((1-power.est)/(power.est - 2))
        ##! adjust varaince when random effects -  ParList$theta[j,g] <- log(var/mu^power.est/exp(1)) / 10 #exp(1) accounts for var=1 attributed to spatial
        ParList$theta[j,g] <- log(var/mu^power.est)
      } else {
        ParList$betad[1,j,g] <- mu
        ParList$theta[j,g] <- log(var)
      }
    }

    y.mat <- Data$Y[Class[,g] == 1,]
    if(Data$fixStruct != 10){
      cor.mat <- cor(y.mat)
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

#' mc.qclass: qclass function from mclust. Defaults used to initate 'E' or 'V' models
#'
#' @param x A numeric vector of observations for classification
#' @param k Integer specifying the number of mixtures
#'
#' @return vclassification vector
#'
#' @examples
#
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
