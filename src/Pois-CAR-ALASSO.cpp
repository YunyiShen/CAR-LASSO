// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h> 
#include <tgmath.h>
using namespace Rcpp;
using namespace arma;

#include <progress.hpp>
#include <progress_bar.hpp>
#include "CAR_LASSO_helper.h"
#include "CAR_ALASSO_helper.h"
#include "ars_pois_helper.h"

/*
 * We would like to develope a Poisson Simulteneous Regressive Graphical LASSO for counting response, 
 * Basic idea was to embed Graphical LASSO into a normal LASSO using the hirechical structure
 *   described by Wang (2012) and Park and Casella 2008
 * 
 * And add another layer of Poisson to account for counting nature of the data
 */

// [[Rcpp::export]]
List Pois_CAR_ALASSO_Cpp(const arma::mat & data, // col composition data, ROW as a sample
                   const arma::mat & design, // design matrix, each ROW as a sample
                   const int n_iter, // how many iterations?
                   const int n_burn_in, // burn in
                   const int thin_by, // thinning?
                   const arma::mat r_beta, // prior on lambda of beta
                   const arma::mat delta_beta,
                   const arma::vec r_Omega, // prior on lambda of Omega
                   const arma::vec delta_Omega,
                   const arma::vec & lambda_diag,// penalty for diagonals 
                   const int ns, const int m, const double emax, // ars parameters
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
  
  arma::mat lambda_beta_mcmc(n_save , k*p); // LASSO parameter for beta and B
  lambda_beta_mcmc += NA_REAL;
  arma::mat lambda_Omega_mcmc(n_save , r_Omega.n_elem); // LASSO parameter for beta and B
  lambda_beta_mcmc += NA_REAL;
  
  arma::vec tau2_curr = randg<arma::vec> (k*p,distr_param(r_beta(0),delta_beta(0))); // current latent variable tau^2, for prior of beta

  arma::vec mu_curr = trans( mean(data) ); // current value of mean
  arma::mat Omega_curr(k,k); // current value of Omega
  Omega_curr = inv(cov(log(data+.1)));
  arma::mat beta_curr(p,k,fill::zeros); // current value of beta
  
  arma::mat Z_curr = log(data+.1); // latent normal variable
  
  
  
  arma::vec lambda2_beta = randg<arma::vec> (k*p,distr_param(1.0,0.01)); // current lambda, for prior of beta
  arma::vec lambda_Omega = randg<arma::vec> (.5*k*(k-1),distr_param(r_Omega(0),delta_Omega(0)));
  
  Progress prog((n_iter+n_burn_in), progress); // progress bar
  
  // main loop
  for(int i = 0 ; i<(n_iter+n_burn_in) ; ++i){
    
    
    if (Progress::check_abort()){
      Rcerr << "keyboard abort\n";
      return(Rcpp::List::create(
          Rcpp::Named("beta") = beta_mcmc,
          Rcpp::Named("mu") = mu_mcmc,
          Rcpp::Named("Omega") = Omega_mcmc,
          Rcpp::Named("lambda_beta") = lambda_beta_mcmc,
          Rcpp::Named("lambda_Omega") = lambda_Omega_mcmc
      ));
    }
    // block update start:
    update_car_lambda_Omega_adp_helper(lambda_Omega,
                                       Omega_curr,
                                       r_Omega,
                                       delta_Omega);    
    
    // Update Latent Zs using ars
    
    update_Z_helper_Pois_CAR(Z_curr,
                             data, design,mu_curr, beta_curr, Omega_curr,
                             k,p,n,ns,m,emax);
    //Rcout << "updated:\n" << Z_curr <<endl;
    //Update betas:
    beta_curr = update_car_beta_helper(Z_curr, design, mu_curr,
                                   tau2_curr, Omega_curr, 
                                   k, p, n);
    
    
    
    // update Omega
    //Rcout<<Z_curr<<endl;
    update_car_Omega_adp_helper2(Omega_curr, Z_curr, design, 
                                     mu_curr, beta_curr,
                                     lambda_Omega,
                                     lambda_diag,
                                     k, p, n);
    
    
    
    // Update mu
    
    mu_curr = update_car_mu_helper(Z_curr,design,beta_curr,
                               Omega_curr, 
                               k, p, n);
    
    
    
  
    // Update tau
    tau2_curr = update_car_tau2_adp_helper(beta_curr,lambda2_beta,
                                   Omega_curr,k,p,n);
    
    
    // Update lambda_beta
    
    update_car_lambda2_beta_adp_helper(lambda2_beta,tau2_curr,
                                        r_beta,delta_beta,k, p);
        
    
    
    // saving the state
    if( (i-n_burn_in)>=0 && (i+1-n_burn_in)%thin_by ==0 ){
      
      
      beta_mcmc.row(i_save) = trans(vectorise(beta_curr));
      Omega_mcmc.row(i_save) = trans( Omega_curr(trimatu_ind( size(Omega_curr) )));
      mu_mcmc.row(i_save) = mu_curr.t();
      
      lambda_beta_mcmc.row(i_save) = trans( sqrt( lambda2_beta));
      lambda_Omega_mcmc.row(i_save) =  lambda_Omega.t();
      
      i_save++;
    }
    
    
    prog.increment();
  }
  return(Rcpp::List::create(
      Rcpp::Named("beta") = beta_mcmc,
      Rcpp::Named("mu") = mu_mcmc,
      Rcpp::Named("Omega") = Omega_mcmc,
      Rcpp::Named("lambda_beta") = lambda_beta_mcmc,
      Rcpp::Named("lambda_Omega") = lambda_Omega_mcmc
  ));
}

