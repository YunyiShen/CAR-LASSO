// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include <RcppArmadillo.h> 
#include <tgmath.h>
using namespace Rcpp;
using namespace arma;

#include <progress.hpp>
#include <progress_bar.hpp>
#include "CAR_LASSO_helper.h"
#include "Probit_helper.h"

/*
 * We would like to develope a Probit CAR LASSO for binary response, 
 * Basic idea was to embed a modified version of Graphical LASSO into a normal LASSO using the hirechical structure
 *   described by Wang (2012) and Park and Casella 2008
 * 
 * 
 */

// [[Rcpp::export]]
List Probit_CAR_LASSO_Cpp(const arma::mat & data, // col composition data, ROW as a sample
                   const arma::mat & design, // design matrix, each ROW as a sample
                   const int n_iter, // how many iterations?
                   const int n_burn_in, // burn in
                   const int thin_by, // thinning?
                   const double r_beta, // prior on lambda of beta
                   const double delta_beta,
                   const double r_Omega,
                   const double delta_Omega,
                   bool progress){
  int k = data.n_cols; // number of nodes
  int p = design.n_cols; //number of predictors
  int n = data.n_rows; // number of samples
  
  int n_save = floor(n_iter/thin_by); //
  int i_save = 0;  
  
  // mcmc matrices:
  arma::mat beta_mcmc(n_save,k * p); // beta mcmc
  beta_mcmc += NA_REAL; 
  
  arma::mat Omega_mcmc(n_save , floor( k * (k+1)/2 )  ); // vectorized column first, but had no diagnol
  Omega_mcmc += NA_REAL;
  
  
  arma::mat mu_mcmc(n_save , k); // mean for node 1 to k
  mu_mcmc += NA_REAL;
  
  arma::mat lambda_mcmc(n_save , 2); // LASSO parameter for beta and B
  lambda_mcmc += NA_REAL;
  
  //arma::mat Z_mcmc(n_save , k*n); // latent normal 
  //Z_mcmc += NA_REAL; // we may not have enough memory to save all latent variables for now.
  
  arma::vec tau2_curr = randg<arma::vec> (k*p,distr_param(r_beta,delta_beta)); // current latent variable tau^2, for prior of beta

  arma::vec mu_curr = trans( mean(data) ); // current value of mean
  arma::mat Omega_curr(k,k); // current value of Omega
  Omega_curr = pinv(cov(data));
  arma::mat beta_curr(p,k,fill::zeros); // current value of beta
  
  arma::mat Z_curr = 0 * data; // latent normal variable
  
  //arma::vec mean_uncertain(k); // for sampling mu
  
  double lambda2_beta = R::rgamma(r_beta,1/delta_beta); // current value of squared LASSO parameter of \beta
  double lambda_Omega = 0;//R::rgamma(r_Omega,1/delta_Omega); // current value of squared LASSO parameter of B
  
  double Omega_r_post = (r_Omega+(k*(k+1)/2));
  double Omega_delta_post;
  
  Progress prog((n_iter+n_burn_in), progress); // progress bar
  
  // main loop
  for(int i = 0 ; i<(n_iter+n_burn_in) ; ++i){
    
    
    if (Progress::check_abort()){
      Rcerr << "keyboard abort\n";
      return(Rcpp::List::create(
          Rcpp::Named("beta") = beta_mcmc,
          Rcpp::Named("mu") = mu_mcmc,
          Rcpp::Named("Omega") = Omega_mcmc,
          Rcpp::Named("lambda") = lambda_mcmc//,
          //Rcpp::Named("Z") = Z_mcmc
      ));
    }
    // block update start:
    Omega_delta_post = (delta_Omega+sum(sum(abs(Omega_curr)))/2);
    lambda_Omega = R::rgamma(Omega_r_post,1/Omega_delta_post);
    
    
    // Update Latent Zs as truncated normal
    
    update_Z_helper_CAR(Z_curr, data,design,mu_curr,beta_curr,
                             Omega_curr,k,p,n);
    
    //Update betas:
    beta_curr = update_car_beta_helper(Z_curr,design,mu_curr,
                                   tau2_curr,Omega_curr, 
                                   k,p,n);
    
    
    // update Omega
    update_car_Omega_helper(Omega_curr,Z_curr, design, 
                                     mu_curr,beta_curr,
                                     lambda_Omega,
                                     k,p,n);
    
    
    
    // Update mu
    
    mu_curr = update_car_mu_helper(Z_curr,design,beta_curr,
                               Omega_curr, 
                               k,p,n);
    
    
    
      
    //Rcout << "detOmega curr in main loop:" << det(Omega_curr) << endl;
    //Rcout << "sum beta in main loop:" <<sum(sum(beta_curr)) <<endl;
    //Rcout << "mean of mu in main loop:" << mean(mu_curr) <<endl;
  
    // Update tau
    tau2_curr = update_car_tau2_helper(beta_curr,lambda2_beta,
                                   Omega_curr,k,p,n);
    
    
    // Update lambda_beta
    
    lambda2_beta = R::rgamma(r_beta+k*p,1/(delta_beta+sum(tau2_curr)/2));
    
    
    // Update lambda_Omega
    
    
    
    // saving the state
    if( (i-n_burn_in)>=0 && (i+1-n_burn_in)%thin_by ==0 ){
      
      //Z_mcmc.row(i_save) = trans(vectorise(Z_curr));
      beta_mcmc.row(i_save) = trans(vectorise(beta_curr));
      Omega_mcmc.row(i_save) = trans( Omega_curr(trimatu_ind( size(Omega_curr) )));
      mu_mcmc.row(i_save) = mu_curr.t();
      
      lambda_mcmc(i_save,0) = sqrt( lambda2_beta);
      lambda_mcmc(i_save,1) = lambda_Omega;
      
      i_save++;
    }
    
    
    prog.increment();
  }
  return(Rcpp::List::create(
      Rcpp::Named("beta") = beta_mcmc,
      Rcpp::Named("mu") = mu_mcmc,
      Rcpp::Named("Omega") = Omega_mcmc,
      Rcpp::Named("lambda") = lambda_mcmc//,
      //Rcpp::Named("Z") = Z_mcmc // it is not a good idea to sace all latent normal
  ));
}

