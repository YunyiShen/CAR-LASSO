// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h> // to use sparse matrix
#include <tgmath.h>
using namespace Rcpp;
using namespace arma;
static const double pi = 3.141592653589793238462643383280;

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
#include "helper.h"

/*
 * We would like to develope a Simulteneous Regressive Graphical LASSO, 
 * Basic idea was to embed Graphical LASSO into a normal LASSO using the hirechical structure
 *   described by Wang (2012) and Park and Casella 2008
 * 
 * 
 */


/*
 * TODO: 
 *  test Gibbs sampler for SAR-lasso
 */


/*Full conditional dist'n of beta was normal with:
 *   
 * Sigma_beta = sigma^2 * (\sum_i [(B-I)X_i]^T[(B-I)X_i] + D_{\tau}^{-1} )^{-1}
 * 
 * mu_beta = (\sum_i [(B-I)X_i]^T[(B-I)X_i] + D_{\tau}^{-1} )^{-1} (\sum_i [(B-I)X_i]^T[(B-I)(Z_i-\mu)]  )
 * 
 */

// tested for dimension compatibility 20200528
arma::vec update_beta_helper(const arma::mat & data,
                             const arma::mat & design,
                             const arma::vec & mu,
                             const arma::vec & tau2,
                             const arma::mat & Omega,
                             int k, int p, int n){
  
  
  arma::vec mu_beta(k*p,fill::zeros);
  arma::mat Q_beta(k*p,k*p,fill::zeros);// percision matrix up to sigma^2 scaling
  Q_beta.diag() += 1/tau2;
  arma::vec res;
  
  arma::sp_mat X_i;
  
  for(int i = 0 ; i < n ; ++i){
    X_i = getDesign_i_helper(design.row(i),k);
    Q_beta +=  X_i.t() * Omega * X_i ;
    mu_beta +=  X_i.t() * Omega * (data.col(i)-mu);
  }
  
  mu_beta = solve(Q_beta,mu_beta);
  res = mvnrnd(mu_beta, inv(Q_beta));
  return(res);
}




// return Omega (first col), eta2 (last col)
arma::mat update_Omega_helper(const arma::mat & data,
                          const arma::mat & design,
                          const arma::vec & mu,
                          const arma::vec & beta,
                          const double & lambda_curr,
                          int k, int p,int n){
  arma::mat beta_mat = reshape(beta,p,k);
  arma::mat Omega;
  arma::mat Y_tilde;
  
  arma::mat expectation = trans(design * beta);
  expectation.each_col() += mu;
  
  Y_tilde = trans( data - expectation);
  
  arma::mat S = Y_tilde.t() * Y_tilde;
  arma::mat Sigma = cov(Y_tilde);
  
  // Concentration matrix and it's dimension:
  arma::mat Omega = pinv(Sigma); // Moore-Penrose inverse
  
  arma::uvec pertub_vec = linspace<uvec>(0,k-1,k); 
  
  arma::uvec Omega_upper_tri = trimatu_ind(size(Omega),1);
  int n_upper_tri = Omega_upper_tri.n_elem;
  
  arma::mat tau_curr(k,k,fill::zeros);
  
  arma::uvec perms_j;
  arma::uvec ind_j(1,fill::zeros);
  arma::vec tauI;
  arma::mat Sigma11;
  arma::mat Sigma12;
  arma::mat S21;
  arma::mat Omega11inv;
  arma::mat Ci;
  arma::mat CiChol;
  arma::mat S_temp;
  arma::mat mui;
  arma::mat gamma;
  
  double gamm_rn;
  arma::mat OmegaInvTemp;
  
  for(int j = 0 ; j < n_upper_tri ; ++j){
    tau_curr(Omega_upper_tri(j)) = 
      rinvGau(sqrt(lambda_curr*lambda_curr/(Omega(Omega_upper_tri(j))*Omega(Omega_upper_tri(j)))),
              lambda_curr*lambda_curr);
  }
  tau_curr = tau_curr + tau_curr.t(); // use symmertric to update lower tri
  
  for(int j = 0 ; j < k ; ++j){
    perms_j = find(pertub_vec!=j);
    ind_j = ind_j.zeros();
    ind_j += j;
    tauI = tau_curr.col(j);
    tauI = tauI(perms_j);
    
    Sigma11 = Sigma(perms_j,perms_j);
    Sigma12 = Sigma.col(j);
    Sigma12 = Sigma12(perms_j);
    
    S21 = S.row(j);
    S21 = S21(perms_j);
    
    Omega11inv = Sigma11-Sigma12 * Sigma12.t() / Sigma(j,j);
    Ci = (S(j,j)+tau_curr) * Omega11inv;
    Ci.diag() += 1/tauI;
    CiChol = chol(Ci);
    
    S_temp = S.col(j);
    S_temp = S_temp(perms_j);
    mui = solve(-Ci,S_temp);
    
    gamma = mui+solve(CiChol,randn(k-1));
    
    // Replacing omega entries
    Omega.submat(perms_j,ind_j) = gamma;
    Omega.submat(ind_j,perms_j) = gamma.t();
    
    gamm_rn = R::rgamma(n/2+1,2/( as_scalar( S(0,0) )+tau_curr));
    Omega(j,j) = gamm_rn + as_scalar( gamma.t() * Omega11inv * gamma);
  }
  
  arma::mat Res(k*k,2,fill::zeros);
  Res.col(0) = vectorise(Omega);
  Res.col(1) = vectorise(tau_curr);
  
  return(Res);
  
}



arma::vec tau2_eta2_update_helper(const arma::vec & beta,
                                  const double & lambda2,
                                  const double & sigma2){
  int n = beta.n_elem;
  arma::vec invtau2(n);
  arma::vec mu_prime = sqrt(lambda2*sigma2/(beta%beta));
  for(int i = 0 ; i < n ; ++i){
    invtau2(i) =  rinvGau(mu_prime(i),lambda2);
  }
  return(1/invtau2);
}

/* 
 * Convention:
 * mcmc obects were vectorized column first, e.g. each row: looks like
 * beta_11,beta_12,beta_13,...,beta_kp
 * 
 */

// [[Rcpp::export]]
List SRG_LASSO_Cpp(const arma::mat & data, // raw composition data, column as a sample
                   const arma::mat & design, // design matrix, each ROW as a sample
                   const int n_iter, // how many iteractions?
                   const int n_burn_in, // burn in
                   const int thin_by, // thinning?
                   const double r_beta, // prior on lambda of beta
                   const double delta_beta,
                   const double r_B,
                   const double delta_B,
                   bool progress){
  int k = data.n_rows; // number of nodes
  int p = design.n_cols; //number of predictors
  int n = data.n_cols; // number of samples
  
  int n_save = floor(n_iter/thin_by); //
  int i_save = 0;  
  
  // mcmc matrices:
  arma::mat beta_mcmc(n_save,k * p); // beta mcmc
  beta_mcmc += NA_REAL; 
  
  arma::mat B_mcmc(n_save , k * (k - 1)); // vectorized column first, but had no diagnol
  B_mcmc += NA_REAL;
  
  arma::mat sigma_mcmc(n_save , 1); //noise standard deviation
  sigma_mcmc += NA_REAL;
  
  arma::mat mu_mcmc(n_save , k); // mean for node 1 to k
  mu_mcmc += NA_REAL;
  
  arma::mat lambda_mcmc(n_save , 2); // LASSO parameter for beta and B
  lambda_mcmc += NA_REAL;
  
  arma::vec tau2_curr = randg<arma::vec> (k*p,distr_param(r_beta,1/delta_beta)); // current latent variable tau^2, for prior of beta
  arma::vec eta2_curr = randg<arma::vec> (k*(k-1),distr_param(r_B,1/delta_B)); // current latent eta^2, for prior of B
  //
  arma::vec mu_curr = mean(data,1); // current value of mean
  
  arma::vec beta_curr(k*p , fill::randu); // current value of beta
  arma::vec B_curr(k*(k-1) , fill::randu); // current value of B
  double sigma_curr = sqrt(1/ R::rgamma(1,1)); // current value of sigma
  
  arma::vec mean_uncertain(k); // for sampling mu
  
  double lambda2_beta = R::rgamma(r_beta,1/delta_beta); // current value of squared LASSO parameter of \beta
  double lambda2_B = R::rgamma(r_B,1/delta_B); // current value of squared LASSO parameter of B
  
  Progress prog((n_iter+n_burn_in), progress); // progress bar
  //Rcout << "flag" <<endl;
  // main loop
  for(int i = 0 ; i<(n_iter+n_burn_in) ; ++i){
    if (Progress::check_abort()){
      Rcerr << "keyboard abort\n";
      return(Rcpp::List::create(
          Rcpp::Named("beta") = beta_mcmc,
          Rcpp::Named("mu") = mu_mcmc,
          Rcpp::Named("B") = B_mcmc,
          Rcpp::Named("sigma") = sigma_mcmc,
          Rcpp::Named("lambda") = lambda_mcmc
      ));
    }
    // block update start:
    
    // Update betas:
    beta_curr = update_beta_helper(data,design,mu_curr,
                                   sigma_curr * sigma_curr,
                                   tau2_curr,B_curr,
                                   k,p,n);
    
    // Update Bs
    B_curr = update_B_helper(data,design,mu_curr,beta_curr,
                             sigma_curr * sigma_curr,
                             eta2_curr,
                             k,n);
      
    // Update mu
    mu_curr = mean(data,1) + sigma_curr/sqrt(n*k) * mean_uncertain.randu();
    
    // Update sigma
    sigma_curr = update_sigma_helper(data, design, mu_curr,
                                     beta_curr, B_curr,
                                     eta2_curr, tau2_curr,
                                     k,p,n);
      
    // Update tau
    tau2_curr = tau2_eta2_update_helper(beta_curr,lambda2_beta,
                                          sigma_curr * sigma_curr);
    
    
    
    // Update eta
    eta2_curr = tau2_eta2_update_helper(B_curr,lambda2_B,
                                        sigma_curr * sigma_curr);
    
    
    
    // Update lambda_beta
    
    lambda2_beta = R::rgamma(r_beta+k*p,1/(delta_beta+sum(tau2_curr)/2));
    
    
    // Update lambda_B
    lambda2_B = R::rgamma(r_B+k*(k-1),1/(delta_B+sum(eta2_curr)/2));
    
    // saving the state
    if( (i-n_burn_in)>=0 && (i+1-n_burn_in)%thin_by ==0 ){
      beta_mcmc.row(i_save) = beta_curr.t();
      B_mcmc.row(i_save) = B_curr.t();
      mu_mcmc.row(i_save) = mu_curr.t();
      sigma_mcmc(i_save) = sigma_curr;
      
      lambda_mcmc(i_save,0) = sqrt( lambda2_beta);
      lambda_mcmc(i_save,1) = sqrt( lambda2_B );
      
      i_save++;
    }
    
    
    prog.increment();
  }
  return(Rcpp::List::create(
      Rcpp::Named("beta") = beta_mcmc,
      Rcpp::Named("mu") = mu_mcmc,
      Rcpp::Named("B") = B_mcmc,
      Rcpp::Named("sigma") = sigma_mcmc,
      Rcpp::Named("lambda") = lambda_mcmc
  ));
}

