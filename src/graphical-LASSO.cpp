// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include <RcppArmadillo.h> // to use sparse matrix
#include <tgmath.h>
using namespace Rcpp;
using namespace arma;



#include <progress.hpp>
#include <progress_bar.hpp>
#include "helper.h"


/*
 * Implements the block Gibbs sampler for the Bayesian graphical lasso
 * introduced in Wang (2012). Samples from the conditional distribution of a 
 * permuted column/row for simulating the posterior distribution for the concentration 
 * matrix specifying a Gaussian Graphical Model
 */

// [[Rcpp::export]]
Rcpp::List Graphical_LASSO_Cpp(const arma::mat & data,
                           const int n_iter,
                           const int n_burn_in,
                           const int thin_by,
                           const double lambda_a, // Gamma prior for LASSO parameter
                           const double lambda_b,
                           bool progress){// Gamma prior for LASSO parameter
  // Sum of product matrix, covariance matrix, n
  int n = data.n_rows;
  arma::mat S = data.t() * data;
  arma::mat Sigma = cov(data);
  
  // Concentration matrix and it's dimension:
  arma::mat Omega = pinv(Sigma); // Moore-Penrose inverse
  int k = Omega.n_rows;
  
  
  // indicator matrix and permutation matrix for looping through columns & rows ("blocks")
  arma::uvec pertub_vec = linspace<uvec>(0,k-1,k); 
  
  // mcmc storage
  int n_store = floor(n_iter/thin_by);
  //arma::mat Sigma_mcmc(n_store,k*k,fill::zeros);
  //Sigma_mcmc += NA_REAL;
  
  arma::mat Omega_mcmc(n_store,floor( k * (k+1)/2 ) ,fill::zeros);
  Omega_mcmc += NA_REAL;
  
  arma::vec lambda_mcmc(n_store,fill::zeros);
  lambda_mcmc + NA_REAL;
  
  int i_save = 0;
  
  // latent tau
  arma::mat tau_curr(k,k,fill::zeros);
  
  double lambda_a_post = (lambda_a+(k*(k+1)/2));
  double lambda_b_post;
  
  double lambda_curr = 0;
  
  arma::uvec Omega_upper_tri = trimatu_ind(size(Omega),1);
  int n_upper_tri = Omega_upper_tri.n_elem;
  
  // some objects needed in block update
  
  arma::uvec perms_j;
  arma::uvec ind_j(1,fill::zeros);
  arma::vec tauI;
  arma::mat Omega11;
  arma::mat Omega11inv;
  arma::mat Omega12;
  
  arma::mat S11;
  arma::mat S12;
  
  arma::mat Ci;
  arma::mat invCi;
  
  arma::mat CiChol;
  arma::mat S_temp;
  arma::mat mui;
  arma::mat gamma;
  double gamm_rn = 0;
  arma::mat OmegaInvTemp;
  
  
  // progress bar
  Progress prog((n_iter+n_burn_in), progress); 
  
  // main iterations
  for(int i = 0 ; i < (n_iter+n_burn_in) ; ++i){
    if (Progress::check_abort()){
      Rcerr << "keyboard abort\n"; // abort using esc
      return(Rcpp::List::create(
          //Rcpp::Named("Sigma") = Sigma_mcmc,
          Rcpp::Named("Omega") = Omega_mcmc,
          Rcpp::Named("lambda") = lambda_mcmc));
    }
    
    // update LASSO parameter lambda
    lambda_b_post = (lambda_b+sum(sum(abs(Omega)))/2);
    lambda_curr = R::rgamma(lambda_a_post,1/lambda_b_post);
    
    // update latent tau (it is symmertric so we only work on upper tri)
    tau_curr.zeros();
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
      
      Omega11 = Omega(perms_j,perms_j);
      Omega12 = Omega(perms_j,ind_j);
      S11 = S(perms_j,perms_j);
      
      S12 = S(perms_j,ind_j);
      //S12 = S21(perms_j);
      
      Omega11inv = inv(Omega11);
      
      Ci = (S(j,j)+lambda_curr) * Omega11inv;
      Ci.diag() += tauI;
      invCi = inv(Ci);
      //CiChol = chol(Ci);
      
      //S_temp = S.col(j);
      //S_temp = S_temp(perms_j);
      mui = -invCi*S12;
      
      gamma = mvnrnd(mui, invCi);
      
      // Replacing omega entries
      Omega.submat(perms_j,ind_j) = gamma;
      Omega.submat(ind_j,perms_j) = gamma.t();
      
      
      
      gamm_rn = R::rgamma(n/2+1,2/( as_scalar( S(j,j) )+lambda_curr));
      Omega(j,j) = gamm_rn + as_scalar( gamma.t() * Omega11inv * gamma);
      
    }
    
    
    
    
    
    // saving current state
    if( (i-n_burn_in)>=0 && (i+1-n_burn_in)%thin_by == 0 ){
      lambda_mcmc(i_save) = lambda_curr;
      //Sigma_mcmc.row(i_save) = Sigma(trimatu_ind(size(Sigma)));
      Omega_mcmc.row(i_save) = trans(Omega(trimatu_ind(size(Omega))));
      i_save++ ;
      
    }
    
    prog.increment();
  }
  
  return(Rcpp::List::create(
      //Rcpp::Named("Sigma") = Sigma_mcmc,
      Rcpp::Named("Omega") = Omega_mcmc,
      Rcpp::Named("lambda") = lambda_mcmc));
  
}

