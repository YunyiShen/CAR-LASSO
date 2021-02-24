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
#include "ars_multinomial_helper.h"

/*
 * We would like to develope a multinomial CAR LASSO for compositional response, 
 * Basic idea was to embed Graphical LASSO into a normal LASSO using the hirechical structure
 *   described by Wang (2012) and Park and Casella 2008
 * 
 * And add another layer of multinomial to account for compositional nature of the data
 */

// [[Rcpp::export]]
List Multinomial_CAR_LASSO_Cpp(const arma::mat & data, // col composition data, ROW as a sample
                   const arma::mat & design, // design matrix, each ROW as a sample
                   const int n_iter, // how many iterations?
                   const int n_burn_in, // burn in
                   const int thin_by, // thinning?
                   const double r_beta, // prior on lambda of beta
                   const double delta_beta,
                   const double r_Omega, // prior on lambda of Omega
                   const double delta_Omega,
                   const int ns, const int m, const double emax, // ars parameters
                   bool progress){
  int k = data.n_cols-1; // number of nodes
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
  
  arma::vec tau2_curr = randg<arma::vec> (k*p,distr_param(1.0,0.01)); // current latent variable tau^2, for prior of beta

  
  arma::mat beta_curr(p,k,fill::zeros); // current value of beta
  
  arma::mat Z_curr = data; // latent normal variable
  Z_curr.each_col() /= (Z_curr.col(k)+.1);
  Z_curr.shed_col(k);
  Z_curr = log(Z_curr+.1);

  arma::vec mu_curr(k,fill::zeros); // current value of mean
  mu_curr = trans(mean(Z_curr));
  arma::mat Omega_curr(k,k); // current value of Omega
  Omega_curr = inv(cov(Z_curr));
  
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
          Rcpp::Named("lambda") = lambda_mcmc
      ));
    }
    
    
    // block update start:
    Omega_delta_post = (delta_Omega+sum(sum(abs(Omega_curr)))/2);
    lambda_Omega = R::rgamma(Omega_r_post,1/Omega_delta_post);
                                          
    
    
    
    //Update betas:
    beta_curr = update_car_beta_helper(Z_curr, design, mu_curr,
                                   tau2_curr, Omega_curr, 
                                   k, p, n);
    
    
    
    // update Omega
    //Rcout<<Z_curr<<endl;
    update_car_Omega_helper(Omega_curr, Z_curr, design, 
                                     mu_curr, beta_curr,
                                     lambda_Omega,
                                     k, p, n);
    
    
    
    // Update mu
    
    mu_curr = update_car_mu_helper(Z_curr,design,beta_curr,
                               Omega_curr, 
                               k, p, n);
    
    
    
  
    // Update tau
    tau2_curr = update_car_tau2_helper(beta_curr,lambda2_beta,
                                   Omega_curr,k,p,n);
    
    
    // Update lambda_beta
    
    lambda2_beta = R::rgamma(r_beta+k*p,1/(delta_beta+sum(tau2_curr)/2));

    //Rcout << i <<endl;    
    //Rcout << Z_curr <<endl;
    update_Z_helper_multinomial_CAR(Z_curr,
                             data, design,mu_curr, beta_curr, Omega_curr,
                             k,p,n,ns,m,emax);


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
  return(Rcpp::List::create(
      Rcpp::Named("beta") = beta_mcmc,
      Rcpp::Named("mu") = mu_mcmc,
      Rcpp::Named("Omega") = Omega_mcmc,
      Rcpp::Named("lambda") = lambda_mcmc
  ));
}

