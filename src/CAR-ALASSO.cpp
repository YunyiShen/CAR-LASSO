// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h> 
#include <tgmath.h>
using namespace Rcpp;
using namespace arma;


// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
//#include "helper.h"
#include "CAR_LASSO_helper.h"
#include "CAR_ALASSO_helper.h"
/*
 * Helper functions for Conditional Auto Regression LASSO, 
 * Basic idea was to embed Graphical LASSO into a normal LASSO using the hirechical structure
 *   described by Wang (2012) and Park and Casella 2008
 * In this model average structure offer some extra information of conditional correlation
 * 
 * A CAR can be reparameterize into a model s.t.: 
 * Y~N(Sigma (Xbeta+mu),Sigma)
 */

/* 
Input:
  @ data: a matrix with column as nodes and row as samples
  @ design: a design matrix of common input to the network, should have same # of rows as data
  @ n_iter: number of saved sampling iterations in the Gibbs sampler
  @ n_burn_in: number of burn in before sampling
  @ thin_by: subsampling steps, integer
  @ r_beta, r_Omega: shape parameter for shrinkage parameter lambda of beta and Omega, 
    for Omega it should be a vector sise as upper.tri, for beta, a matrix as size of beta
  @ delta_beta, delta_Omega: RATE parameter for lambda prior
  @ progress: whether to show a progress bar from C++

Output:
  A list with component:
  @ beta: a matrix with each row as an MCMC sample, 
    columns are the vectorization of beta, 
    while beta matrix has p row and k columns
  @ mu: a matrix with each row as an MCMC sample, columns are intercept vectors
  @ Omega: a matrix with each row as an MCMC sample, 
    columns are the upper diagnol entries of precision matrix Omega
  @ lambda: a matrix with only row columns, first was for beta, second was for Omega
    each row was an MCMC sample of shrinkage parameter lambda


*/



// [[Rcpp::export]]
List CAR_ALASSO_Cpp(const arma::mat & data, // col composition data, ROW as a sample
                   const arma::mat & design, // design matrix, each ROW as a sample
                   const int n_iter, // how many iterations?
                   const int n_burn_in, // burn in
                   const int thin_by, // thinning?
                   const arma::mat r_beta, // prior on lambda of beta
                   const arma::mat delta_beta,
                   const arma::vec r_Omega, // prior on lambda of Omega
                   const arma::vec delta_Omega,
                   bool progress) {// whether to report progress
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
  arma::mat centered_data = data;
  centered_data.each_row() -= mu_curr.t();
  arma::mat Omega_curr(k,k); // current value of Omega
  Omega_curr = inv(cov(data));
  arma::mat beta_curr = solve( design.t()*design,design.t()*(centered_data*Omega_curr)); // current value of beta
  
  
  arma::vec lambda2_beta = randg<arma::vec> (k*p,distr_param(r_beta(0),delta_beta(0))); // current lambda, for prior of beta
  arma::vec lambda_Omega = randg<arma::vec> (.5*k*(k+1),distr_param(r_Omega(0),delta_Omega(0)));
  
  Progress prog((n_iter+n_burn_in), progress); // progress bar
  
  // main loop
  for(int i = 0 ; i<(n_iter+n_burn_in) ; ++i){
    
    
    if (Progress::check_abort()){
      Rcerr << "keyboard abort\n";
      return(
        Rcpp::List::create(
          Rcpp::Named("beta") = beta_mcmc,
          Rcpp::Named("mu") = mu_mcmc,
          Rcpp::Named("Omega") = Omega_mcmc,
          Rcpp::Named("lambda_beta") = lambda_beta_mcmc,
          Rcpp::Named("lambda_Omega") = lambda_Omega_mcmc
        )
      );
    }
    // block update start:
	  
	// Update lambda_Omega
    update_car_lambda_Omega_adp_helper(lambda_Omega,
                                       Omega_curr,
                                       r_Omega,
                                       delta_Omega);
    
    //Update betas:
    beta_curr = update_car_beta_helper(data, design, mu_curr,
                                   tau2_curr, Omega_curr, 
                                   k, p, n);
    
    
    // update Omega
    
    update_car_Omega_adp_helper(Omega_curr, data, design, 
                                     mu_curr, beta_curr,
                                     lambda_Omega,
                                     k, p, n);
    
    
    
    // Update mu
    
    mu_curr = update_car_mu_helper(data,design,beta_curr,
                               Omega_curr, 
                               k, p, n);
    
    
  
    // Update tau2 for beta
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
  return(
    Rcpp::List::create(
          Rcpp::Named("beta") = beta_mcmc,
          Rcpp::Named("mu") = mu_mcmc,
          Rcpp::Named("Omega") = Omega_mcmc,
          Rcpp::Named("lambda_beta") = lambda_beta_mcmc,
          Rcpp::Named("lambda_Omega") = lambda_Omega_mcmc
    )
  );
}

