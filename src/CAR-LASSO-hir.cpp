// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h> 
#include <tgmath.h>
using namespace Rcpp;
using namespace arma;


// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
//#include "helper.h"
#include "CAR_LASSO_helper.h"
#include "ars_pois_helper.h"
#include "ars_multinomial_helper.h"
#include "Probit_helper.h"

/*
 * Main sampling functions for Conditional Auto Regression LASSO, 
 * Basic idea was to embed Graphical LASSO into a normal LASSO using the hirechical structure
 *   described by Wang (2012) and Park and Casella 2008
 * In this model average structure offer some extra information of conditional correlation
 * 
 * This one was not adaptive, i.e. lambda_Omega was a fixed for all entries of Omega and 
 *   similar for betas
 * 
 * A CAR can be reparameterize into a model s.t.: 
 * Z~N(Sigma (Xbeta+mu),Sigma)
 * Y~Pois, multinomial, probit
 */

/* 
Input:
  @ data: a matrix with column as nodes and row as samples
  @ design: a design matrix of common input to the network, should have same # of rows as data
  @ link: code for the link function, 0: identity, 1: log for pois, 2: logit for multinomial, 3: probit for binary
  @ n_iter: number of saved sampling iterations in the Gibbs sampler
  @ n_burn_in: number of burn in before sampling
  @ thin_by: subsampling steps, integer
  @ r_beta, r_Omega: shape parameter for shrinkage parameter lambda of beta and Omega
  @ delta_beta, delta_Omega: RATE parameter for lambda prior
  @ int ns,int m,int emax: parameters for ars
  @ progress: whether to show a progress bar from C++

Output:
  A list with component:
  @ beta: a matrix with each row as an MCMC sample, 
    columns are the vectorization of beta, 
    while beta matrix has p row and k columns
  @ mu: a matrix with each row as an MCMC sample, columns are intercept vectors
  @ Omega: a matrix with each row as an MCMC sample, 
    columns are the upper diagnol entries of precision matrix Omega
  @ lambda: a matrix with only two columns, first was for beta, second was for Omega
    each row was an MCMC sample of shrinkage parameter lambda


*/



// [[Rcpp::export]]
List CAR_LASSO_hir_Cpp(const arma::mat & data, // col composition data, ROW as a sample
                   const arma::mat & design, // design matrix, each ROW as a sample
                   const int link, 
                   const int n_iter, // how many iterations?
                   const int n_burn_in, // burn in
                   const int thin_by, // thinning?
                   const double r_beta, // prior on lambda of beta
                   const double delta_beta,
                   const double r_Omega, // prior on lambda of Omega
                   const double delta_Omega,
                   int ns,int m,int emax,
                   bool progress) {// whether to report progress
  int k = data.n_cols; // number of nodes
  int p = design.n_cols; //number of predictors
  int n = data.n_rows; // number of samples

  bool flag = true;
  arma::mat Z;
  if(link==0){
    Z = data;
    flag = false;
  }

  if(link==1){
    Z = log(data + .1);
    flag = false;
  }
  
  if(link==2){
    Z = data; // latent normal variable
    Z.each_col() /= (Z.col(k-1)+.1);
    Z.shed_col(k-1);
    Z = log(Z+.1);
    k = k-1;// composition has less df
    flag = false; 
  }
  
  if(link==3){
    Z = data;
    flag = false;
  }

  if(flag){
    Rcerr << "Unknown link code, this should not happen in R interface, if you see it, please file an issue.\n" << endl;
  }

  
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
  
  //arma::vec tau2_curr = randg<arma::vec> (k*p,distr_param(r_beta,delta_beta)); // current latent variable tau^2, for prior of beta
  arma::vec tau2_curr = randg<arma::vec> (k*p,distr_param(1.0,.01)); // current latent variable tau^2, for prior of beta
  
  //Rcout << tau2_curr <<endl;
  arma::vec mu_curr = trans( mean(Z) ); // current value of mean
  arma::mat centered_Z = Z;
  centered_Z.each_row() -= mu_curr.t();
  arma::mat Omega_curr(k,k); // current value of Omega
  Omega_curr = inv(cov(Z));
  arma::mat beta_curr = solve( design.t()*design,design.t()*(centered_Z*Omega_curr)); // current value of beta
  
  
  
  
  double lambda2_beta = randg<double>(distr_param(r_beta,1/delta_beta)); // current value of squared LASSO parameter of \beta
  double lambda_Omega = randg<double>(distr_param(r_Omega,1/delta_Omega)); // current value of squared LASSO parameter of B
  
  double Omega_r_post = (r_Omega+(k*(k+1)/2));
  double Omega_delta_post;
  
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
          Rcpp::Named("lambda") = lambda_mcmc
        )
      );
    }
    // update Z
        if(link == 1){
      update_Z_helper_Pois_CAR(Z,data,design,mu_curr, beta_curr,Omega_curr, k,p,n,ns,m,emax);
    }

    if(link == 2){
      update_Z_helper_multinomial_CAR(Z,data,design,mu_curr, beta_curr,Omega_curr, k,p,n,ns,m,emax);
    }

    if(link == 3){
      update_Z_helper_CAR(Z, data,design,mu_curr,beta_curr,
                             Omega_curr,k,p,n);
    }

    // block update start:
	  
	// Update lambda_Omega
    Omega_delta_post = (delta_Omega + sum(sum(abs(Omega_curr)))/2);
    lambda_Omega = randg<double>(distr_param(Omega_r_post, 1/Omega_delta_post));
    
    //Update betas:
    beta_curr = update_car_beta_helper(Z, design, mu_curr,
                                   tau2_curr, Omega_curr, 
                                   k, p, n);
    
    
    // update Omega
    
    update_car_Omega_helper(Omega_curr, Z, design, 
                                     mu_curr, beta_curr,
                                     lambda_Omega,
                                     k, p, n);
    
    
    
    // Update mu
    
    mu_curr = update_car_mu_helper(Z,design,beta_curr,
                               Omega_curr, 
                               k, p, n);
    
    
  
    // Update tau2 for beta
    tau2_curr = update_car_tau2_helper(beta_curr,lambda2_beta,
                                   Omega_curr,k,p,n);
    
    
    // Update lambda_beta
    
    lambda2_beta = randg<double>(distr_param(r_beta+k*p,1/(delta_beta+sum(tau2_curr)/2)));
    
    
    
    
    
    
    // saving the state
    if( (i-n_burn_in)>=0 && (i+1-n_burn_in)%thin_by ==0 ){
      
      beta_mcmc.row(i_save) = trans(vectorise(beta_curr));
      Omega_mcmc.row(i_save) = trans( Omega_curr(trimatu_ind( size(Omega_curr) )));
      mu_mcmc.row(i_save) = mu_curr.t();
      
      lambda_mcmc(i_save,0) = sqrt( lambda2_beta);
      lambda_mcmc(i_save,1) = lambda_Omega;
      
      i_save++;
    }
    
    
    prog.increment();
  }
  return(
    Rcpp::List::create(
      Rcpp::Named("beta") = beta_mcmc,
      Rcpp::Named("mu") = mu_mcmc,
      Rcpp::Named("Omega") = Omega_mcmc,
      Rcpp::Named("lambda") = lambda_mcmc
    )
  );
}

