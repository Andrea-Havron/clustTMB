#' genInit: generates starting values for cluster model
#'
#' @param init.method method as string used to generate initial values. Options include: "random", "hc"
#' @parma Data list containing TMB data objects
#' @param dim integer
#' @param data.trans string identifying any data transformation
#' @param hc.options hierarchical clustering settings. Position[1] is model name -see mclust.options("hcModelName"). Position[2] is data transformation -see mclust.options("hcUse").
#' @param mix.method string stating initialization method
#' @param r.tol relative convergence tolerance
#' @param true.class vector of true classification integers for testing models
#'
#' @importFrom stats cutree gaussian hclust runif rmultinom nlminb cor
#' @importFrom mclust unmap hclass hc hcVVV hcE
#' @importFrom e1071 skewness
#' @importFrom cluster daisy pam
#' @importFrom geoR boxcoxfit
#' @importFrom clustMixType kproto
#'
#' @return list
#' @export
#'
#' @examples
genInit <- function(init.method, Data, dim.list, data.trans = NA, hc.options = c("VVV", "SVD"), #original default was VARS, changed to match mclust default 7/20/2020
                     mix.method = NULL, r.tol = 1e-10, true.class = NULL){

  #list2env(dim.list, environment(genStart))
  n.i <- dim.list$n.i
  n.j <- dim.list$n.j
  n.t <- dim.list$n.t
  n.g <- dim.list$n.g
  n.f <- dim.list$n.f
  n.v <- dim.list$n.v

  if(n.j == 1) hc.options[1] <- 'E'

  # #Apply any data transformations
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

  #Apply classification method
  if(init.method =="random"){
    pi.init <- rep(1/n.g,n.g)
    Class <- t(rmultinom(n.i, 1, pi.init)) ## FIXME: write my own rmultinom function to avoid dependency here
  }

  if(init.method == "hc"){
    #Class <- unmap(hclass(hc(y),n.g))
    Class <- unmap(hclass(hc(y, modelName = hc.options[1], use = hc.options[2]),n.g))
    pi.init <-  rep(1/n.g,n.g)#apply(Class,2,function(x) sum(x)/nrow(Data$Y))
    classify = apply(Class,1,function(x) which(x==1))

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
      tmp.pa[,j] = ifelse(y[,j]==0, "A", "P")
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

  if(Data$fixStruct == 36){
    nl <- n.j*n.f - (n.f*(n.f-1))/2
  } else {
    nl <- (n.j^2-n.j)/2
  }

  #setup ParList
  ParList <- list(
    betag = matrix(0,1,(n.g-1)),
    betad = array(0, dim = c(1,n.j,n.g)),
    betapz = array(numeric(0)),
    theta = matrix(0,n.j,n.g),
    thetaf = matrix(0,n.j,n.g),
    ld = matrix(0,nl,n.g),
    Hg_input = matrix(0,2,(n.g-1)),
    Hd_input = array(0, dim = c(2,n.j,n.g)), #fix in cpp - should be n.f, not n.j
    ln_kappag = rep(0, (n.g-1)),
    ln_taug = rep(0,(n.g-1)),
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
    Omega_vfg = array(0, dim = c(n.v,n.f,n.g))
  )
  if(sum(Data$reStruct[1,])==0){#no random effects in gating, turn on beta, o.w. beta in gating is mapped out
    ParList$betag <- matrix(log(pi.init[1:n.g-1]/(1 - sum(pi.init[1:n.g-1]))), nrow = 1)
  }

  #Update initial values based on classification
  for(g in 1:n.g){
    for(j in 1:n.j){
      y.sub <- Data$Y[Class[,g] == 1,j]
      #  if(init.model == "empirical"){
      if(sum(y.sub) == 0){
        mu <- 0.01
        var <- 0.01
        power.est <- 1.05
      } else {
        mu <- mean(y.sub)
        var <- var(y.sub)
        power.est <-  e1071::skewness(y.sub)*mu/sqrt(var)  # Clark and Thayer, 2004
      }
      #ideally this will be based on link function not family!!
      if(Data$family == 700){
        ParList$betad[,j,g] <- log(mu)
        if(power.est >= 2) power.est <- 1.95
        if(power.est <= 1) power.est <- 1.05
        ParList$thetaf[j,g] <- log((1-power.est)/(power.est - 2))
        ParList$theta[j,g] <- log(var/mu^power.est/exp(1)) / 10 #exp(1) accounts for var=1 attributed to spatial
      } else {
        ParList$betad[,j,g] <- mu
        ParList$theta[j,g] <- log(var)
      }


      y.mat <- Data$Y[Class[,g] == 1,]
      if( Data$fixStruct == 36 ){
        cor.mat <- cor(y.mat)
        if(sum(is.na(cor.mat))>0){
          idx.na <- which(is.na(cor.mat), arr.ind = TRUE)
          tmp.pa <- matrix(0, nrow(y.mat), n.j)
          for(j in 1:n.j){
            tmp.pa[y.mat[,j]>0,j] <- 1
          }

          for(n in 1:nrow(idx.na)){
            tmp.ctab <- cbind(
              c(nrow(tmp.pa[tmp.pa[,idx.na[n,1]]==1 & tmp.pa[,idx.na[n,2]]==1,]),
                nrow(tmp.pa[tmp.pa[,idx.na[n,1]]==0 & tmp.pa[,idx.na[n,2]]==1,])),
              c(nrow(tmp.pa[tmp.pa[,idx.na[n,1]]==1 & tmp.pa[,idx.na[n,2]]==0,]),
                nrow(tmp.pa[tmp.pa[,idx.na[n,1]]==0 & tmp.pa[,idx.na[n,2]]==0,]))
            )
            tmp.ctab[tmp.ctab==0] <- 1
            cor.mat[idx.na[n,2], idx.na[n,1]] <-
              (tmp.ctab[1,1]*tmp.ctab[2,2] - tmp.ctab[1,2]*tmp.ctab[2,1]) /
              sqrt(sum(tmp.ctab[1,])*sum(tmp.ctab[2,])*sum(tmp.ctab[,1])*sum(tmp.ctab[,2]))
          }
        }
        L.mat <- t(chol(cor.mat))[,1:n.f]
        cnt <- start <- 1
        for(f. in 1:n.f){
          for(j. in start:n.j){
            ParList$ld[cnt,g] <- L.mat[j.,f.]
            cnt <- cnt + 1
          }
          start <- start + 1
        }

      } else{
        if( Data$fixStruct == 20 ){
          ParList$ld[,g] <- 0
        } else {
          corr.mat <- cor(y.mat)
          corr.mat[is.na(corr.mat)] <- 0
          ParList$ld[,g] <- corr.mat[which(lower.tri(corr.mat)==TRUE, arr.ind = TRUE)]
        }
      }
      # if(init.model == "opt"){
      #   p.init <- e1071::skewness(y.sub)*mean(y.sub)/sd(y.sub)
      #   mu.init <- mean(y.sub)
      #   phi.init <- var(y.sub)/mu.init^p.init
      #   mod <- nlminb(c(log(mu.init),log(phi.init), log((1-p.init)/(p.init - 2))), tweed.fun, dat = y.sub)
      #   ParList$betad[1,j,g] <- mod$par[1]
      #   ParList$theta[j,g] <- mod$par[2]
      #   ParList$thetaf[j,g] <- mod$par[3]
      # }
      # }
    }
  }
  gen.init <- list(parms = ParList, class = classify)
  return(gen.init)

}

