// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include <RcppArmadillo.h> // to use sparse matrix
#include <tgmath.h>
using namespace Rcpp;
using namespace arma;

#include <progress.hpp>
#include <progress_bar.hpp>
#include "SRG_LASSO_helper.h"
#include "ars_pois_helper.h"
#include "ZIP_helper.h"

// [[Rcpp::export]]
List ZIP_SRG_LASSO_Cpp(const arma::mat & data, // col composition data, ROW as a sample
                        const arma::mat & design_P, // design matrix, each ROW as a sample
                        const arma::mat & design_ZI
                        const int n_iter, // how many iterations?
                        const int n_burn_in, // burn in
                        const int thin_by, // thinning?
                        const double r_beta_P, // prior on lambda of beta
                        const double delta_beta_P,
                        const double r_Omega_P,
                        const double delta_Omega_P,
                        const double r_beta_ZI, // prior on lambda of beta
                        const double delta_beta_ZI,
                        const double r_Omega_ZI,
                        const double delta_Omega_ZI,
                        const int ns, const int m, const double emax, // ars parameters
                        bool progress){
  int k = data.n_cols; // number of nodes
  int p = design_P.n_cols; //number of predictors
  int n = data_P.n_rows; // number of samples
  
  int n_save = floor(n_iter/thin_by); //
  int i_save = 0;  
  
  // mcmc matrices:
  arma::mat beta_P_mcmc(n_save,k * p); // beta mcmc
  beta_P_mcmc += NA_REAL; 
  
  arma::mat beta_ZI_mcmc(n_save,k * p); // beta mcmc
  beta_ZI_mcmc += NA_REAL; 
  
  arma::mat Omega_P_mcmc(n_save , floor( k * (k+1)/2 )  ); // vectorized column first, but had no diagnol
  Omega_P_mcmc += NA_REAL;
  
  
  arma::mat Omega_ZI_mcmc(n_save , floor( k * (k+1)/2 )  ); // vectorized column first, but had no diagnol
  Omega_ZI_mcmc += NA_REAL;
  
  arma::mat mu_P_mcmc(n_save , k); // mean for node 1 to k
  mu_P_mcmc += NA_REAL;
  
  arma::mat mu_P_mcmc(n_save , k); // mean for node 1 to k
  mu_P_mcmc += NA_REAL;
  
  arma::mat lambda_mcmc(n_save , 4); //mu_P, S_P, beta_ZI,S_ZI // LASSO parameter for beta and B
  lambda_mcmc += NA_REAL;
  
  //arma::mat Z_mcmc(n_save , k*n); // latent normal 
  //Z_mcmc += NA_REAL;
  
  arma::vec tau2_P_curr = randg<arma::vec> (k*p,distr_param(r_beta_P,delta_beta_P)); // current latent variable tau^2, for prior of beta
  arma::vec tau2_ZI_curr = randg<arma::vec> (k*p,distr_param(r_beta_ZI,delta_beta_ZI)); // current latent variable tau^2, for prior of beta
  
  
  
  arma::vec mu_curr = trans( mean(data) ); // current value of mean
  arma::mat Omega_curr(k,k); // current value of Omega
  Omega_curr = inv(cov(log(data+.1)));
  arma::mat beta_curr(p,k,fill::zeros); // current value of beta
  
  arma::mat Z_curr = log(data+.1); // latent normal variable
  
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
          Rcpp::Named("lambda") = lambda_mcmc,
          //Rcpp::Named("Z") = Z_mcmc
      ));
    }
    // block update start:
    Omega_delta_post = (delta_Omega+sum(sum(abs(Omega_curr)))/2);
    lambda_Omega = R::rgamma(Omega_r_post,1/Omega_delta_post);
    
    
    // Update Latent Zs using ars
    //Rcout << "current Z : \n" << Z_curr <<endl;
    update_Z_helper_Pois_reg(Z_curr,
                             data, design,mu_curr, beta_curr, Omega_curr,
                             k,p,n,ns,m,emax);
    //Rcout << "updated:\n" << Z_curr <<endl;
    //Update betas:
    beta_curr = update_beta_helper(Z_curr,design,mu_curr,
                                   tau2_curr,Omega_curr, 
                                   k,p,n);
    
    
    // update Omega
    //Rcout<<Z_curr<<endl;
    Omega_curr = update_Omega_helper(Z_curr, design, 
                                     mu_curr,beta_curr,
                                     lambda_Omega,
                                     k,p,n);
    
    
    
    // Update mu
    
    mu_curr = update_mu_helper(Z_curr,design,beta_curr,
                               Omega_curr, 
                               k,p,n);
    
    
    
    
    //Rcout << "detOmega curr in main loop:" << det(Omega_curr) << endl;
    //Rcout << "sum beta in main loop:" <<sum(sum(beta_curr)) <<endl;
    //Rcout << "mean of mu in main loop:" << mean(mu_curr) <<endl;
    
    // Update tau
    tau2_curr = update_tau2_helper(beta_curr,lambda2_beta,
                                   Omega_curr,k,p,n);
    
    
    // Update lambda_beta
    
    lambda2_beta = R::rgamma(r_beta+k*p,1/(delta_beta+sum(tau2_curr)/2));
    
    
    // Update lambda_Omega
    
    
    
    // saving the state
    if( (i-n_burn_in)>=0 && (i+1-n_burn_in)%thin_by ==0 ){
      
      Z_mcmc.row(i_save) = trans(vectorise(Z_curr));
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
      Rcpp::Named("lambda") = lambda_mcmc,
      //Rcpp::Named("Z") = Z_mcmc
  ));
}




