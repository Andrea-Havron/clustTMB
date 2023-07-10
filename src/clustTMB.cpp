#define TMB_LIB_INIT R_init_clustTMB
#include <TMB.hpp>

using namespace density;
using namespace Eigen;
using namespace tmbutils;
using namespace R_inla;


template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

template<class Type>
Type lddirichlet(vector<Type> q, vector<Type> alpha){
  Type ans = 0;
  Type n_a = alpha.size();
  ans += lgamma(alpha.sum());
  for(int a=0; a<n_a; a++){
    ans -= lgamma(alpha(a));
    ans -= (alpha(a)-Type(1)) * log(q(a));
  }
  return ans;
}

//adapted from STAN - need to test
template<class Type>
vector<Type> rdirichlet(Type alpha){
  vector<Type> ans(alpha.size());
  if(alpha.minCoeff() < 1){
    vector<Type> log_y(alpha.size());
    Type log_u;
    for (int i=0; i<alpha.size(); ++i) {
      log_u = log(runif(Type(0), Type(1)));
      log_y(i) = log(rgamma(alpha(i) + Type(1), Type(1))) + log_u / alpha(i);
    }
    Type log_sum_y = log_y(0);
    for (int i=1; i<alpha.size(); ++i) {
      log_sum_y = logspace_add(log_sum_y, log_y(i));
    }
    for (int i=0; i<alpha.size(); ++i) {
      ans(i) = exp(log_y(i) - log_sum_y);
    }
  }
  if(alpha.minCoeff() == 1 | alpha.minCoeff() > 1){
    vector<Type> y(alpha.size());
    y = rgamma(alpha, 1e-7);
    for (int i=0; i<alpha.size(); ++i) {
      ans(i) = y(i)/y.sum();
    }
  }
  return(ans);
}

template<class Type>
Type log_sum_exp(vector<Type> x){
  //x is in log space
  Type ans;
  int K = x.size();
  Type tmax = max(x);

  ans = x(0) - tmax;
  for(int k=1; k<K; k++){
    ans = logspace_add(ans, x(k) - tmax);
  }
  ans += tmax;
  return(ans);
}

template<class Type>
vector<Type> norm_exp(vector<Type> x){
  //x is in log space
  vector<Type> ans(x.size());
  ans = exp(x - log_sum_exp(x));
  return(ans);
}

/* Define distributional families *///from glmmTMB
enum valid_family{
  gaussian_family = 0,
  Gamma_family = 300,
  lognormal_family = 600,
  tweedie_family = 700
};


enum valid_reStruct {
  // st covariance
  na_reStruct = 0,
  norm_reStruct = 1,
  ar1_reStruct = 2,
  gmrf_reStruct = 3,
  gmrf_speedup_reStruct = 4
};


enum valid_fixStruct{
  Univariate_fixStruct = 10,
  Diag_fixStruct = 20,
  General_fixStruct = 30
};

enum valid_rrStruct{
  full_rrStruct = 0,
  reduce_rrStruct = 1
};

/* Define link functions *///from glmmTMB
enum  valid_link{
  log_link = 0,
  logit_link = 1,
  probit_link = 2,
  inverse_link = 3,
  cloglog_link = 4,
  identity_link = 5,
  sqrt_link = 6,
  yoejin_boxcox_link = 7
};

enum valid_loglike{
  postMarginal_loglike = 0,
  marginal_loglike = 1,
  conditional_loglike = 2
};

template<class Type>
Type reNll(array<Type> reVec, vector<Type> parmVec, int reStruct, bool do_simulate = false){
  Type ans = 0;
  Type ar1sd;
  switch(reStruct){
      case na_reStruct:
        ans += 0;
        break;
      case norm_reStruct:
        for(int i=0; i<reVec.size(); i++){
          ans -= dnorm(reVec(i),Type(0), parmVec(0), true);
          if(do_simulate){
            reVec(i) = rnorm(Type(0), parmVec(0));
          }
        }
        break;
      case ar1_reStruct:
        ar1sd = sqrt(parmVec(1) * 1/(1-pow(parmVec(0),2)));
        ans += SCALE(AR1(parmVec(0)), ar1sd)(reVec);
         if(do_simulate){
          vector<Type>simVec(reVec.size());
          AR1(parmVec(0)).simulate(simVec);
          reVec = simVec * ar1sd;
        }
        break;
      default:
        error("reNLL method not implemented");
    }
    return ans;
}

template<class Type>
Type spNll(array<Type> reVec, vector<Type> parmVec, spde_aniso_t<Type> Spde, int reStruct, bool do_simulate = false){
  Type ans = 0;
  matrix<Type> H(2,2);
  H(0,0) = exp(parmVec(2));
  H(1,0) = parmVec(3);
  H(0,1) = parmVec(3);
  H(1,1) = ( 1 + pow(parmVec(3),2) )/exp(parmVec(2));
  int nq = reVec.size();

  SparseMatrix<Type> Q(nq,nq);
  switch(reStruct){
      case na_reStruct:
        ans += 0;
        break;
      case gmrf_reStruct:
        Q = Q_spde(Spde, parmVec(0), H);
        ans += SCALE(GMRF(Q),1/parmVec(1))( reVec );
        if(do_simulate){
           vector<Type>simVec(reVec.size());
           GMRF(Q).simulate(simVec);
           reVec = simVec/parmVec(1);
        }
        break;
      case gmrf_speedup_reStruct:
        Q = Q_spde(Spde, parmVec(0), H);
        ans += SCALE(GMRF(Q, false),1/parmVec(1))( reVec );
        if(do_simulate){
           vector<Type>simVec(reVec.size());
           GMRF(Q, false).simulate(simVec); //drop normalizing constant
           reVec = simVec/parmVec(1);
        }
        break;
      default:
       error("spNll method not implemented");
   }
   return ans;
}

//inverse links *///from glmmTMB
template<class Type>
Type inverse_linkfun(Type eta, int link){
  Type ans;
  switch(link){
    case log_link:
      ans = exp(eta);
      break;
    case identity_link:
      ans = eta;
      break;
    case logit_link:
      ans = invlogit(eta);
      break;
    case probit_link:
      ans = pnorm(eta);
      break;
    default:
      error("Link not implemented");
  } // end switch
  return ans;
}

template<class Type>
vector<Type> inverse_mlogit(vector<Type> eta){
  vector<Type> ans(eta.size()+1);
  Type pi_sum = 0;
  Type b = max(eta);
  for(int g=0; g<eta.size(); g++){
    ans(g) = exp(eta(g)-b)/(exp(-b) + exp(eta-b).sum());
    pi_sum += ans(g);
  }
  ans(eta.size()) = Type(1) - pi_sum;
  return ans;
}


template<class Type>
matrix<Type> corrmat_fun(vector<Type> l, int nj ){
  matrix<Type> ans(nj,nj);
  ans.setIdentity();
  vector<Type> l_input = invlogit(l)*Type(2)-Type(1);
  int cnt = 0;
  //matrix<Type> L(nj,nj);
  //L.setIdentity();
  //for(int j1=1; j1<nj; j1++){
  //  Type Norm2 = L(j1,j1);
  //  for(int j2=0; j2<=j1-1; j2++){
  //    L(j1,j2) = l(cnt);
  //    Norm2 += L(j1,j2)*L(j1*j2);
  //  }
  //  for(int j2=0; j2<=j1; j2++){
  //    L(j1,j2) /= sqrt(Norm2);
  //  }
  //}
  //matrix<Type> ans = L * L.transpose();
  for(int j1=0; j1<(nj-1); j1++){
    for(int j2=(j1+1); j2<nj; j2++){
      ans(j1,j2) = l_input(cnt);
      ans(j2,j1) = l_input(cnt);
      cnt += 1;
    }
  }
  return ans;
}

template<class Type>
vector<Type> rmultinom(vector<Type> prob){
  int n_g = prob.size();
  vector<Type> ans(n_g);

  //initiate random generation
  Type N_cnt = Type(1);
  ans(0) = rbinom(N_cnt, prob(0));
  Type ans_cnt = ans(0);
  Type prob_cnt = prob(0);
  N_cnt -= ans_cnt;
  for(int g=1; g<(n_g-1); g++){
  	ans(g) = rbinom(N_cnt, prob(g)/(1-prob_cnt));
  	N_cnt -= ans(g);
 	  prob_cnt += prob(g);
  }
  ans(n_g) = N_cnt;
  return ans;
}

template<class Type>
Type prec_fun(Type kappa, Type tau, int rrStruct){
  Type ans = 0;
  switch( rrStruct ){
    case reduce_rrStruct:
      ans = 1 / (2 * M_PI * kappa);
      break;
    case full_rrStruct:
      ans = tau;
      break;
  }
  return ans;
}

template<class Type>
array<Type> rr_fun(array<Type> x, matrix<Type>l, int nj, int rrStruct){
  //factor analysis from adapted from Thorson, James T. VAST, https://github.com/James-Thorson-NOAA/VAST 
  vector<int> d = x.dim;
  int ni = d(0);
  int nf = d(1);
  int ng = d(2);
  array<Type> ans(ni,nj,ng);
  int cnt;
  switch(rrStruct){
    case reduce_rrStruct:
      for(int g=0; g<ng; g++){
        cnt = 0;
        matrix<Type> L(nj, nf);
        for(int f=0; f<nf; f++){
        for(int j=0; j<nj; j++){
          if(j>=f){
            L(j,f) = l(cnt,g);
            cnt ++;
          } else {
            L(j,f) = 0.0;
          }
        }
        } //!!!! length of ld may not be equal for each g - need to fix code to account for this -- read i as list?
        for(int i=0; i<ni; i++){
          for(int j=0; j<nj; j++){
            for(int f=0; f<nf; f++){
               ans(i,j,g) += x(i,f,g) * L(j,f);
            }
          }
        }
      }
      break;
    case full_rrStruct:
      ans = x;
      break;
  }
  return ans;
}


template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_ARRAY( Y ); // Data matrix: i obs x j species
  DATA_IVECTOR( t ); //temporal index
  DATA_MATRIX( Xd ); //Static Covariates that drive distribution means: i x k for k covariates
  DATA_MATRIX( Xg ); //Static Covariates that drive group membership: i x l for l covariates
  DATA_MATRIX( Xpz ); //Static Covariates that drive zero inflation: i x m for m covariates
  DATA_SPARSE_MATRIX( A ); // Map vertex v to site s
  DATA_SPARSE_MATRIX( A_proj ); // Map vertex v to projection grid n
  DATA_INTEGER( doProj );
  DATA_MATRIX( Xd_proj ); //Distribution covariates used for prediction
  DATA_MATRIX( Xg_proj ); //Cluster covariates used for prediction

  DATA_STRUCT( spde, spde_aniso_t );
  DATA_INTEGER( family );
  DATA_INTEGER( link );
  DATA_INTEGER( loglike );
  DATA_INTEGER( fixStruct  );
  DATA_IVECTOR( rrStruct ); //1: random, 2:spatial
  DATA_IMATRIX( reStruct  ); //row 1: gating; row 2: expert; col 1: spatial; col 2: temporal; col3: overdispersion

  //DATA_ARRAY_INDICATOR( keep, Y );

  vector<int> y_dim = Y.dim; //dims: i*t,j (i: no. obs; t: no. time events; j: no. spp)

  // Parameters
  PARAMETER_MATRIX( betag ); // Group membership covariate coefficients; dims: kg,(g-1)
  PARAMETER_ARRAY( betad ); // Distribution covariate coefficients; dims: k,j,g
  PARAMETER_ARRAY( betapz ); //Zero-inflation covariate coefficients; dims: m,j,g
  PARAMETER_MATRIX( theta );  //variance parameters; dims: j,g
  PARAMETER_MATRIX( thetaf ); // If NOT tweedie, these are mapped to 0, o.w. tweedie power parameter
  PARAMETER_MATRIX( logit_corr_fix );
  PARAMETER_MATRIX( ld_rand ); //!!!! length of ld may not be equal for each g - need to fix code to account for this -- read in as vector with separate id vector
  PARAMETER_MATRIX( ld_sp ); //!!!! length of ld may not be equal for each g - need to fix code to account for this -- read in as vector with separate id vector
  PARAMETER_MATRIX( Hg_input ); //dims: 2, n_g-1
  PARAMETER_ARRAY( Hd_input ); //dims: 2, (Multivariate:n_j, RR:n_f), n_g
  PARAMETER_VECTOR( ln_kappag );   //dims: n_g-1
  PARAMETER_MATRIX( ln_kappad );   //dims: (Multivariate:n_j, RR:n_f), n_g
  PARAMETER_MATRIX( ln_taud );     //dims: (Multivariate:n_j, RR:n_f), n_g
  PARAMETER_VECTOR( logit_rhog );    //dims: n_g-1
  PARAMETER_MATRIX( logit_rhod ); //dims: n_j, n_g
  PARAMETER_VECTOR( ln_sigmaup );   //dims: ng-1
  PARAMETER_MATRIX( ln_sigmaep ); //dims: n_j, n_g
  //PARAMETER_VECTOR( ln_sigmau );    //dims: n_g-1 //!not implemented yet
  PARAMETER_MATRIX( ln_sigmav );   //dims: (Multivariate:n_j, RR:n_f), n_g
  PARAMETER_ARRAY( upsilon_tg ); //dims: n_t, n_g-1
  PARAMETER_ARRAY( epsilon_tjg ); //dims: n_t, n_j, n_g
  //PARAMETER_ARRAY( u_ig ); //dims: n_i, n_g-1 //!not implemented yet
  PARAMETER_ARRAY( v_ifg ); //dims: n_i, (Multivariate:n_j, RR:n_f), n_g
  PARAMETER_ARRAY( Gamma_vg ); //dims: n_v, n_g-1
  PARAMETER_ARRAY( Omega_vfg ); //dims: n_v, (Multivariate:n_j, RR:n_f), n_g


  //Type nll = 0;
  Type nll_re = 0;
  Type nll_data = 0;
 // parallel_accumulator<Type> nll(this);

  vector<int> omega_dim = Omega_vfg.dim;
  vector<int> v_dim = v_ifg.dim;

  int n_g = omega_dim(2);
  int n_i = y_dim(0);
  int n_j = y_dim(1);
  int n_k = Xd.row(0).size(); //no. covariates in distribution regression
  int n_m = Xpz.row(0).size(); //no. covariate in zero-inflation
  int n_v = omega_dim(0);
  int n_x = A_proj.col(0).size();

  bool pz_flag = (betapz.size() > 0);


  vector<Type> kappag = exp(ln_kappag);
  vector<Type> taug = 1 / (sqrt(4*M_PI)*kappag);
  vector<Type> rhog = invlogit(logit_rhog);
  vector<Type> sigmaup = exp(ln_sigmaup);
  //vector<Type> sigmau = exp(ln_sigmau); //!not implemented yet
  matrix<Type> kappad(n_j, n_g);
  matrix<Type> taud(n_j, n_g);
  matrix<Type> rhod(n_j, n_g);
  matrix<Type> sigmaep(n_j, n_g);
  matrix<Type> sigmav(n_j, n_g);
  for(int g=0; g<n_g; g++){
    for(int j=0; j<n_j; j++){
      rhod(j,g) = invlogit(logit_rhod(j,g));
      sigmaep(j,g) = exp(ln_sigmaep(j,g));
    }
    for(int j=0; j<ln_kappad.col(0).size(); j++){
      kappad(j,g) = exp(ln_kappad(j,g));
      taud(j,g) = exp(ln_taud(j,g));
    }
    for(int j=0; j<ln_sigmav.col(0).size(); j++){
      sigmav(j,g) = exp(ln_sigmav(j,g));
    }
  }

  //// tweedie priors
  //switch(family){
  //   case tweedie_family:
  //   for(int j=0; j<n_j; j++){
  //    for(int g=0; g<n_g; g++){
  //      nll -= log(Type(1)/(Type(100)-Type(0))); //U(0,100) prior on phi
  //      nll -= log(Type(1)/(Type(1.99)-Type(1.01))); //U(1.01,1.99) prior on power
  //      nll -= dnorm(betad(0,j,g), Type(0), Type(1000), true); //prior on betag
  //    }
  //   }
  //   break;
  //  default:
  //    nll = 0;
  //}




  //// Random Effects likelihood ===================================================

  ////// Cluster probability ==================================================
  //////// Gating
  vector<Type> tParmg(2);
  tParmg.setZero();
  vector<Type> spParmg(4);
  spParmg.setZero();
  //vector<Type> ovParmg(1);
  //ovParmg.setZero();


  for(int g=0; g<(n_g-1); g++){
    spParmg(0) = kappag(g);
    spParmg(1) = taug(g);
    spParmg(2) = Hg_input(0,g);
    spParmg(3) = Hg_input(1,g);
    tParmg(0) = rhog(g);
    tParmg(1) = sigmaup(g);
    //ovParmg(0) = sigmau(g);
    //Spatial effect

    nll_re += spNll(Gamma_vg.col(g), spParmg, spde, reStruct(0,0), this->do_simulate);

    //Temporal effect
    nll_re += reNll(upsilon_tg.col(g), tParmg, reStruct(0,1), this->do_simulate);

    //Overdispersion
    //nll_re += reNll(u_ig.col(g), ovParmg, reStruct(0,2), this->do_simulate);
  }

  //////// Expert
  vector<Type> tParmd(2);
  tParmd.setZero();
  vector<Type> spParmd(4);
  spParmd.setZero();
  vector<Type> ovParmd(1);
  ovParmd.setZero();

  for(int g=0; g<n_g; g++){
    for(int f=0; f<omega_dim(1); f++){
      spParmd(0) = kappad(f,g);
      spParmd(1) = prec_fun(kappad(f,g), taud(f,g), rrStruct(1));
      spParmd(2) = Hd_input(0,f,g);
      spParmd(3) = Hd_input(1,f,g);
      //Spatial effect
      nll_re += spNll((Omega_vfg.col(g)).col(f), spParmd, spde, reStruct(1,0), this->do_simulate);
    }
    for(int j=0; j<n_j; j++){
      tParmd(0) = rhod(j,g);
      tParmd(1) = sigmaep(j,g);
      ovParmd(0) = sigmav(j,g);
      //Temporal effect
      nll_re += reNll((epsilon_tjg.col(g)).col(j), tParmd, reStruct(1,1), this->do_simulate);

      //Overdispersion
      nll_re += reNll((v_ifg.col(g)).col(j), ovParmd, reStruct(1,2), this->do_simulate);
    }
  }
  matrix<Type> gamma_vg(n_v, (n_g-1));
  for(int g=0; g<(n_g-1); g++){
    gamma_vg.col(g) = Gamma_vg.col(g);
  }

  matrix<Type> Gamma_ig = A * gamma_vg;

  //Apply any rank reduction
  array<Type> Omega_vjg(n_v, n_j, n_g);
  array<Type> Omega_ijg (n_i, n_j, n_g);
  Omega_ijg.setZero();
  Omega_vjg = rr_fun(Omega_vfg, ld_sp, n_j, rrStruct(1));

  for(int g=0; g<n_g; g++){
    for(int i=0; i<n_i; i++){
      for(int j=0; j<n_j; j++){
        for(int v=0; v<n_v; v++){
          Omega_ijg(i,j,g) += A.coeffRef(i,v) * Omega_vjg(v,j,g);  //-check to make sure indexing array correctly
        }
      }
    }
  }

  array<Type> v_ijg(n_i, n_j, n_g);
  v_ijg = rr_fun(v_ifg, ld_rand, n_j, rrStruct(0));



  ////Linear predictor and Links ========================================

  //Cluster probability linear predictor and link
  matrix<Type> etag = Xg*betag;
  for(int g=0; g<(n_g-1); g++){
    for(int i=0; i<n_i; i++){
      etag(i,g) += Gamma_ig(i,g) + upsilon_tg(t(i),g);// + u_ig(i,g);
    }
  }
  matrix<Type> pi(n_i, n_g);
  for(int i=0; i<n_i; i++){
    pi.row(i) = inverse_mlogit(vector<Type>(etag.row(i)));
  }

  //Observation linear predictor and link
  array<Type> etad(n_i, n_j, n_g);
  array<Type> etapz(n_i, n_j, n_g);
  array<Type> mu(n_i, n_j, n_g);
  array<Type> pz(n_i, n_j, n_g);
  for(int g=0; g<n_g; g++){
    for(int i=0; i<n_i; i++){
      for(int j=0; j<n_j; j++){
        for(int k=0; k<n_k; k++){
          etad(i,j,g) += Xd(i,k) * betad(k,j,g);
        }
        etad(i,j,g) += Omega_ijg(i,j,g) + epsilon_tjg(t(i),j,g) + v_ijg(i,j,g);
        if(pz_flag){
          for(int m=0; m<n_m; m++){
            etapz(i,j,g) += Xpz(i,m) * betapz(m,j,g);
          }
        }
        mu(i,j,g) = inverse_linkfun(etad(i,j,g), link);
        pz(i,j,g) = invlogit(etapz(i,j,g));
      }
    }
  }


  //Variance link
  matrix<Type> var(n_j, n_g);
  for(int g=0; g<n_g; g++){
    var.col(g) = exp(vector<Type>(theta.col(g)));
  }

  ////Observation Likelihood ===============================================
  Type s1;
  Type s2;
  Type s3;
  //tweedie power link

  matrix<Type> tmp_ll(n_i, n_g); //dim: i,g
  tmp_ll.setZero();
  array<Type>Corr_mat_g(n_j,n_j,n_g);
  Corr_mat_g.setZero();
  vector<Type> residual(n_j);
  residual.setZero();
    //  if(!isNA(Y(i,j))){
  switch(family){
    case gaussian_family:
      switch(fixStruct){
        case Univariate_fixStruct:
        case Diag_fixStruct:
          for(int i=0; i<n_i; i++){
            for(int g=0; g<n_g; g++){
              for(int j=0; j<n_j; j++){
                tmp_ll(i,g) += dnorm(Y(i,j), mu(i,j,g), sqrt(var(j,g)), true);
              }
            }
          }
          break;
        case General_fixStruct:
          for(int g=0; g<n_g; g++){
            vector<Type> sds(n_j);
            matrix<Type> Corr_mat = corrmat_fun(vector<Type>(logit_corr_fix.col(g)), n_j );
            Corr_mat_g.col(g) = Corr_mat;
            for(int j=0; j<n_j; j++){
              sds(j) = sqrt(var(j,g));
            }
            MVNORM_t<Type> neg_log_dmvnorm(Corr_mat);
            for(int i=0; i<n_i; i++){
              for(int j=0; j<n_j; j++){
               residual(j) = Y(i,j) - mu(i,j,g);
              }
              tmp_ll(i,g) -= VECSCALE(neg_log_dmvnorm, sds)(residual); //record positve likelihood
            }
          }
          break;
        default:
        error("Covariance structure not implemented");
      } //end method switch
      break;
    case Gamma_family:
      switch(fixStruct){
        case Univariate_fixStruct:
        case Diag_fixStruct:
	      for(int i=0; i<n_i; i++){
	        for(int g=0; g<n_g; g++){
	          for(int j=0; j<n_j;j++){
	            s1 = var(j,g);
	            s2 = mu(i,j,g) / var(j,g);
	            if(pz_flag){
	              //Type logit_pz = etapz(i,j,g) ;
	              //Type log_pz   = -logspace_add( Type(0) , -logit_pz );
	              //Type log_1mpz = -logspace_add( Type(0) ,  logit_pz );
	              if(Y(i,j) == Type(0)){
	                //tmp_ll(i,g) += logpz;
	                tmp_ll(i,g) += log(pz(i,j,g));
	              } else {
	                // Is this correct? - tmp_ll(i,g) += logspace_add( log_1mpz, dgamma(Y(i,j), s1, s2, true) );
	                tmp_ll(i,g) += log(Type(1) - pz(i,j,g)) + dgamma(Y(i,j), s1, s2, true);
	              }
	            } else {
	              tmp_ll(i,g) += dgamma(Y(i,j), s1, s2, true);
	            }
	          }
	        }
	      }
	      break;
	    default:
	      	error("General covariance structure not implemented for Gamma");
	      }
    case lognormal_family:
      switch(fixStruct){
        case Univariate_fixStruct:
        case Diag_fixStruct:
          for(int i=0; i<n_i; i++){
            for(int g=0; g<n_g; g++){
              for(int j=0; j<n_j; j++){
                if(pz_flag){
                //Type logit_pz = etapz(i,j,g) ;
                //Type log_pz   = -logspace_add( Type(0) , -logit_pz );
                //Type log_1mpz = -logspace_add( Type(0) ,  logit_pz );
                  if(Y(i,j) == Type(0)){
                    //tmp_ll(i,g) += logpz;
                    tmp_ll(i,g) += log(pz(i,j,g));
                  } else {
                    // Is this correct? - tmp_ll(i,g) += logspace_add( log_1mpz, dgamma(Y(i,j), s1, s2, true) );
                    tmp_ll(i,g) += log(Type(1) - pz(i,j,g)) + dnorm(log(Y(i,j)), mu(i,j,g), sqrt(var(j,g)), true) - log(Y(i,j));
                  }
                } else {
                  tmp_ll(i,g) += dnorm(log(Y(i,j)), mu(i,j,g), sqrt(var(j,g)), true) - log(Y(i,j));
                }
              }
            }
          }
          break;
        case General_fixStruct:
        for(int g=0; g<n_g; g++){
          vector<Type> sds(n_j);
          matrix<Type> Corr_mat = corrmat_fun(vector<Type>(logit_corr_fix.col(g)), n_j );
          Corr_mat_g.col(g) = Corr_mat;
          for(int j=0; j<n_j; j++){
            sds(j) = sqrt(var(j,g));
          }
          MVNORM_t<Type> neg_log_dmvnorm(Corr_mat);
          for(int i=0; i<n_i; i++){
            for(int j=0; j<n_j; j++){
             residual(j) = log(Y(i,j)) - mu(i,j,g);
            }
            tmp_ll(i,g) -= VECSCALE(neg_log_dmvnorm, sds)(residual); //record positve likelihood
            //lognormal correction
            for(int j=0; j<n_j; j++){
              tmp_ll(i,g) -= log(Y(i,j));
            }
          }
        }
        break;
        default:
          error("Covariance structure not implemented for lognormal");
        } // end method switch
        break;
    case tweedie_family:
      switch(fixStruct){
        case Univariate_fixStruct:
        case Diag_fixStruct:
	      for(int g=0; g<n_g; g++){
	        for(int j=0; j<n_j;j++){
	          s2 = var(j,g);
	          s3 = invlogit(thetaf(j,g)) + Type(1); //p, 1<p<2
	          for(int i=0; i<n_i; i++){
	            s1 = mu(i,j,g);
	            tmp_ll(i,g) += dtweedie(Y(i,j), s1, s2, s3, true);
	          }
	        }
	      }
	      break;
	    default:
	      error("General covariance structure not implemented for tweedie");
	  }
    } // end switch
    //Add zero inflation -- from glmmTMB - is this correct for delta models?


  ////Marginalize z by hand   ==================================================================
  matrix<Type> z_ig(n_i,n_g);
  z_ig.setZero();
  matrix<Type> zhat(n_i,n_g);
  zhat.setZero();
  matrix<Type> ln_fy(n_i,n_g);
  vector<Type> newpi(n_g);
  newpi.setZero();
  //vector<Type> ln_fy_sum(n_i);
  //Type test = 0;
  //vector<Type> denom(n_i);
  //denom.setZero();
  for(int i=0; i<n_i; i++){
  //  denom(i) = log(pi(i,0)) + tmp_ll(i,0);
  //  for(int g=1; g<n_g; g++){
  //    denom(i) = logspace_add(denom(i), log(pi(i,g)) + tmp_ll(i,g));
  //  } //end g loop - complete denom sum before next g loop

  //  for(int g=0; g<n_g; g++){
  //    zhat(i,g) = exp(log(pi(i,g)) + tmp_ll(i,g) - denom(i));
  //  } // end g loop - fill in zhat row before z_ig assignment
  //  Type M = max(vector<Type>(zhat.row(i)));
  //  for(int g=0; g<n_g; g++){
  //    z_ig(i,g) = CppAD::CondExpEq(zhat(i,g), M, Type(1), Type(0));
  //  } //end g loop

    for(int g=0; g<n_g; g++){
      ln_fy(i,g) = log(pi(i,g)) + tmp_ll(i,g);
    }
   // nll -= log_sum_exp(vector<Type>(ln_fy.row(i)));
    for(int g=0; g<n_g; g++){
      zhat(i,g) = exp( ln_fy(i,g) - log_sum_exp(vector<Type>(ln_fy.row(i))) );
      newpi(g) += zhat(i,g);
    }
    Type zmax = max(vector<Type>(zhat.row(i)));
    for(int g=0; g<n_g; g++){
      z_ig(i,g) = CppAD::CondExpEq(zhat(i,g), zmax, Type(1), Type(0));
    }

    //for(int g=0; g<n_g; g++){
    //  ln_fy(i,g) = z_ig(i,g) * ln_fy(i,g);
    //}

  }
  newpi = newpi/n_i;

  switch(loglike){
    case postMarginal_loglike: //0
      for(int i=0; i<n_i; i++){
        for(int g=0; g<n_g; g++){
          ln_fy(i,g) = log(newpi(g)) + tmp_ll(i,g);
        }
        nll_data -= log_sum_exp(vector<Type>(ln_fy.row(i)));
       // rep2 = ln_fy;
       // REPORT(rep2);
      }
      break;
    case marginal_loglike: //1
      for(int i=0; i<n_i; i++){
        nll_data -= log_sum_exp(vector<Type>(ln_fy.row(i)));
      }
      break;
    case conditional_loglike: //5
      for(int i=0; i<n_i; i++){
        for(int g=0; g<n_g;g++){
          nll_data -= z_ig(i,g) * ln_fy(i,g);
        }
      }
      break;
    default:
      error("ll not specified");
  } //end switch



  //Type test2 = nll;

  //REPORT(test);
  //REPORT(test2);
  REPORT(newpi);





  //// Cluster Likelihood =======================================================================
 // for(int i=0; i<n_i; i++){
 // 	for(int g=0; g<n_g; g++){
 //       nll -= z_ig(i,g) * (log(pi(i,g)) + tmp_ll(i,g));
 //   } // end g loop
 // } // end i loop

  SIMULATE{
    vector<Type> res_sim(n_j);
    matrix<Type> Corr_mat(n_j,n_j);

    for(int i=0; i<n_i; i++){
      for(int g=0; g<n_g; g++){
        if(pi(i,g) == max(vector<Type>(pi.row(i)))){
          switch(family){
            case gaussian_family:
              switch(fixStruct){
                case Univariate_fixStruct:
                case Diag_fixStruct:
                  for(int j=0; j<n_j; j++){
                    Y(i,j) = rnorm(mu(i,j,g), sqrt(var(j,g)));
                  }
                  break;
                case General_fixStruct:
                  for(int j1=0; j1<n_j; j1++){
                    for(int j2=0; j2<n_j; j2++){
                      Corr_mat(j1,j2) = Corr_mat_g(j1,j2,g);
                    }
                  }
                  MVNORM(Corr_mat).simulate(res_sim);
                  for(int j=0; j<n_j; j++){
                    res_sim(j) = res_sim(j)*sqrt(var(j,g));
                    Y(i,j) = res_sim(j) + mu(i,j,g);
                  }
                  break;
                default:
                  error("Method not implemented");
              } //end switch
              break;
             case lognormal_family:
              switch(fixStruct){
                case Univariate_fixStruct:
                case Diag_fixStruct:
                  for(int j=0; j<n_j; j++){
                    Y(i,j) = exp(rnorm(mu(i,j,g), sqrt(var(j,g))) );
                  }
                  break;
                case General_fixStruct:
                  for(int j1=0; j1<n_j; j1++){
                    for(int j2=0; j2<n_j; j2++){
                      Corr_mat(j1,j2) = Corr_mat_g(j1,j2,g);
                    }
                  }
                  MVNORM(Corr_mat).simulate(res_sim);
                  for(int j=0; j<n_j; j++){
                    res_sim(j) = res_sim(j)*sqrt(var(j,g));
                    Y(i,j) = exp ( res_sim(j) + mu(i,j,g) );
                  }
                  break;
                default:
                  error("Method not implemented");
              } //end switch
              break;
            case tweedie_family:
              switch(fixStruct){
                case Univariate_fixStruct:
                case Diag_fixStruct:
	              for(int j=0; j<n_j; j++){
	                s1 = mu(i,j,g);
	                s2 = var(j,g);
	                s3 = invlogit(thetaf(j,g)) + Type(1); //p, 1<p<2
	                Y(i,j) = rtweedie(s1, s2, s3);
	              }
	              break;
                default:
                  error("Method not implemented");
              } //end switch
            default:
              error("Family not implmented");
          } // end switch
        } // end conditional
        if(pz_flag == 1){
          for(int j=0; j<n_j; j++){
            Y(i,j) = Y(i,j)*rbinom(Type(1), Type(1)-pz(i,j,g));
          }
        } //end pz flag
      } // end g loop
    } //end i loop
  } // end simulate



 //// Classification ========================================================================
  vector<int> classification(n_i);
  for(int i=0; i<n_i; i++){
    for(int g=0; g<n_g; g++){
      //if(pi(i,g) == max(vector<Type>(pi.row(i)))){ //CHECK! pi_ig and z_ig produce slightly different results! - why?
      if(z_ig(i,g) == 1){
          classification(i) = g;
      }
    }
  }

  //// Prediction ===========================================================================
  if(doProj){
    matrix<Type> Gamma_xg = A_proj * gamma_vg;


    matrix<Type> etapi_pred = Xg_proj * betag;
    matrix<Type> pi_pred(n_x,n_g);
    vector<int> Class_pred(n_x);
    for(int g=0; g<(n_g-1); g++){
      etapi_pred.col(g) += Gamma_xg.col(g);
    }
    for(int n=0; n<n_x; n++){
      pi_pred.row(n) = inverse_mlogit(vector<Type>(etapi_pred.row(n)));
      for(int g=0; g<n_g; g++){
        if(pi_pred(n,g) == max(vector<Type>(pi_pred.row(n)))) Class_pred(n) = g;
      }
    }

    //density prediction only implemented for spatial
    array<Type> Omega_xjg(n_x, n_j, n_g);
    matrix<Type> tmp_omega2(n_v, n_j);
    for(int g=0; g<n_g; g++){
      tmp_omega2.setZero();
      for(int v=0; v<n_v; v++){
        for(int j=0; j<n_j; j++){
          tmp_omega2(v,j) = Omega_vjg(v,j,g);
        }
      }
      Omega_xjg.col(g) = A_proj * tmp_omega2;
    }

    matrix<Type> etad_pred(n_x, n_j);
    etad_pred.setZero();
    matrix<Type> mu_pred(n_x, n_j);
    for(int n=0; n<n_x; n++){
      for(int j=0; j<n_j; j++){
        for(int k=0; k<n_k; k++){
          etad_pred(n,j) += Xd_proj(n,k) * betad(k,j,Class_pred(n));
        }
        etad_pred(n,j) += Omega_xjg(n,j,Class_pred(n));
        mu_pred(n,j) = inverse_linkfun(etad_pred(n,j), link);
      }
    }
    REPORT( Gamma_xg );
    REPORT( pi_pred );
    REPORT( Class_pred );
    REPORT( etad_pred );
    REPORT( mu_pred );
  }




  //// NLL =====================================================================================

  Type nll = nll_re + nll_data;

  //// Report ==================================================================================

  REPORT( var );
  REPORT( pi );
  REPORT( z_ig );
  REPORT( mu );
  REPORT( etad );
  REPORT( nll_re );
  REPORT( nll_data );
  //REPORT( nll );
  REPORT( zhat );
  REPORT( classification );

  REPORT( tmp_ll );
  REPORT( Gamma_ig );
  REPORT( Omega_ijg );
  REPORT( Corr_mat_g );
  REPORT( residual );
  REPORT( Omega_vjg );
  REPORT( v_ijg );

  SIMULATE{
    REPORT( Y );
    REPORT( upsilon_tg );
    REPORT( epsilon_tjg );
  //  REPORT( u_ig );
    REPORT( v_ifg );
    REPORT( Gamma_vg );
    REPORT( Omega_vfg );
    REPORT( Omega_vjg );
  }

 //// ===================================================================================================
  return nll;
}
