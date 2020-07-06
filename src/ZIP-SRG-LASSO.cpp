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
  int p_P = design_P.n_cols; //number of predictors
  int p_ZI = design_ZI.n_cols; //number of predictors for ZI
  int n = data.n_rows; // number of samples
  
  int n_save = floor(n_iter/thin_by); //
  int i_save = 0;  
  
  // mcmc matrices:
  arma::mat beta_P_mcmc(n_save,k * p_P); // beta mcmc
  beta_P_mcmc += NA_REAL; 
  
  arma::mat beta_ZI_mcmc(n_save,k * p_ZI); // beta mcmc
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
  
  arma::vec tau2_P_curr = randg<arma::vec> (k*p_P,distr_param(r_beta_P,delta_beta_P)); // current latent variable tau^2, for prior of beta
  arma::vec tau2_ZI_curr = randg<arma::vec> (k*p_ZI,distr_param(r_beta_ZI,delta_beta_ZI)); // current latent variable tau^2, for prior of beta
  
  
  
  arma::vec mu_P_curr = trans( mean(data) ); // current value of mean in Pois
  arma::mat Omega_P_curr(k,k); // current value of Omega
  
  arma::vec mu_ZI_curr = trans( mean(data) ); // current value of mean in ZI
  arma::mat Omega_ZI_curr(k,k); // current value of Omega
  
  Omega_P_curr = inv(cov(log(data+.1)));
  Omega_ZI_curr = Omega_P_curr;
  
  arma::mat beta_P_curr(p,k,fill::zeros); // current value of beta for Pois
  arma::mat beta_ZI_curr(p,k,fill::zeros); // current value of beta for ZI
  
  arma::mat Z_P_curr = log(data+.1); // latent normal variable
  arma::mat Pois_latent = data;
  arma::mat Z_ZI_curr = 0 * data + 1 ;
  
  //arma::vec mean_uncertain(k); // for sampling mu
  
  double lambda2_beta_P = R::rgamma(r_beta,1/delta_beta); // current value of squared LASSO parameter of \beta
  double lambda_Omega_P = 0;//R::rgamma(r_Omega,1/delta_Omega); // current value of squared LASSO parameter of B
  
  double lambda2_beta_ZI = R::rgamma(r_beta,1/delta_beta); // current value of squared LASSO parameter of \beta
  double lambda_Omega_ZI = 0;//R::rgamma(r_Omega,1/delta_Omega); // current value of squared LASSO parameter of B
  
  double Omega_P_r_post = (r_Omega_P+(k*(k+1)/2));
  double Omega_P_delta_post;
  
  double Omega_ZI_r_post = (r_Omega_ZI+(k*(k+1)/2));
  double Omega_ZI_delta_post;
  
  Progress prog((n_iter+n_burn_in), progress); // progress bar
  
  // main loop
  for(int i = 0 ; i<(n_iter+n_burn_in) ; ++i){
    
    
    if (Progress::check_abort()){
      Rcerr << "keyboard abort\n";
      return(Rcpp::List::create(
          Rcpp::Named("beta_P") = beta_P_mcmc,
          Rcpp::Named("mu_P") = mu_P_mcmc,
          Rcpp::Named("Omega_P") = Omega_P_mcmc,
          Rcpp::Named("beta_ZI") = beta_ZI_mcmc,
          Rcpp::Named("mu_ZI") = mu_ZI_mcmc,
          Rcpp::Named("Omega_ZI") = Omega_ZI_mcmc,
          Rcpp::Named("lambda") = lambda_mcmc,
          //Rcpp::Named("Z") = Z_mcmc
      ));
    }
    // block update start:
    Omega_P_delta_post = (delta_Omega_P+sum(sum(abs(Omega_P_curr)))/2);
    lambda_P_Omega = R::rgamma(Omega_P_r_post,1/Omega_P_delta_post);
    Omega_ZI_delta_post = (delta_Omega_ZI+sum(sum(abs(Omega_ZI_curr)))/2);
    lambda_ZI_Omega = R::rgamma(Omega_ZI_r_post,1/Omega_ZI_delta_post);
    
    ZIP_update_Pois_latent_helper(Pois_latent, data,
                                  Z_P_curr, Z_ZI_curr,
                                  k, p_P, n)
    
    // Update Latent Zs using ars
    //Rcout << "current Z : \n" << Z_curr <<endl;
    update_Z_helper_Pois_reg(Z_P_curr,
                             Pois_latent, design_P,
                             mu_P_curr, beta_P_curr, Omega_P_curr,
                             k,p_P,n,ns,m,emax);
    //Rcout << "updated:\n" << Z_curr <<endl;
    //Update betas:
    beta_P_curr = update_beta_helper(Z_P_curr,design_P,mu_P_curr,
                                   tau2_P_curr,Omega_P_curr, 
                                   k,p_P,n);
    
    
    // update Omega
    //Rcout<<Z_curr<<endl;
    Omega_P_curr = update_Omega_helper(Z_P_curr, design_P, 
                                     mu_P_curr,beta_P_curr,
                                     lambda_Omega_P,
                                     k,p_P,n);
    
    
    
    // Update mu
    
    mu_P_curr = update_mu_helper(Z_P_curr,design_P,beta_P_curr,
                               Omega_P_curr, 
                               k,p_P,n);
    
    
    
    
    //Rcout << "detOmega curr in main loop:" << det(Omega_curr) << endl;
    //Rcout << "sum beta in main loop:" <<sum(sum(beta_curr)) <<endl;
    //Rcout << "mean of mu in main loop:" << mean(mu_curr) <<endl;
    
    // Update tau
    tau2_P_curr = update_tau2_helper(beta_P_curr,lambda2_beta_P,
                                   Omega_P_curr,k,p_P,n);
    
    
    // Update lambda_beta
    
    lambda2_P_beta = R::rgamma(r_beta_P+k*p_P,1/(delta_beta_P+sum(tau2_P_curr)/2));
    
    
    // Start the ZI part //
    
    ZIP_update_Z_ZI_helper(Z_ZI_curr,
                           data, design_ZI,mu_ZI_curr,beta_ZI_curr,
                           Omega_ZI_curr,Pois_latent,
                           k, p_ZI, n);
    
    //Update betas:
    beta_ZI_curr = update_beta_helper(Z_ZI_curr,design_ZI,mu_ZI_curr,
                                     tau2_ZI_curr,Omega_ZI_curr, 
                                     k,p_ZI,n);
    
    
    // update Omega
    //Rcout<<Z_curr<<endl;
    Omega_ZI_curr = update_Omega_helper(Z_ZI_curr, design_ZI, 
                                       mu_ZI_curr,beta_ZI_curr,
                                       lambda_Omega_ZI,
                                       k,p_ZI,n);
    
    
    
    // Update mu
    
    mu_ZI_curr = update_mu_helper(Z_ZI_curr,design_ZI,beta_ZI_curr,
                                 Omega_ZI_curr, 
                                 k,p_ZI,n);
    
    
    
    // Update tau
    tau2_ZI_curr = update_tau2_helper(beta_ZI_curr,lambda2_beta_ZI,
                                     Omega_ZI_curr,k,p_ZI,n);
    
    
    // Update lambda_beta
    
    lambda2_ZI_beta = R::rgamma(r_beta_ZI+k*p_ZI,1/(delta_beta_ZI+sum(tau2_ZI_curr)/2));
    
    
    
    
    // saving the state
    if( (i-n_burn_in)>=0 && (i+1-n_burn_in)%thin_by ==0 ){
      
      beta_P_mcmc.row(i_save) = trans(vectorise(beta_P_curr));
      Omega_P_mcmc.row(i_save) = trans( Omega_P_curr(trimatu_ind( size(Omega_P_curr) )));
      mu_P_mcmc.row(i_save) = mu_P_curr.t();
      
      beta_ZI_mcmc.row(i_save) = trans(vectorise(beta_ZI_curr));
      Omega_ZI_mcmc.row(i_save) = trans( Omega_ZI_curr(trimatu_ind( size(Omega_ZI_curr) )));
      mu_ZI_mcmc.row(i_save) = mu_ZI_curr.t();
      
      lambda_mcmc(i_save,0) = sqrt( lambda2_beta_P);
      lambda_mcmc(i_save,1) = lambda_Omega_P;
      lambda_mcmc(i_save,2) = sqrt( lambda2_beta_ZI);
      lambda_mcmc(i_save,3) = lambda_Omega_ZI;
      
      i_save++;
    }
    
    
    prog.increment();
  }
  return(Rcpp::List::create(
      Rcpp::Named("beta_P") = beta_P_mcmc,
      Rcpp::Named("mu_P") = mu_P_mcmc,
      Rcpp::Named("Omega_P") = Omega_P_mcmc,
      Rcpp::Named("beta_ZI") = beta_ZI_mcmc,
      Rcpp::Named("mu_ZI") = mu_ZI_mcmc,
      Rcpp::Named("Omega_ZI") = Omega_ZI_mcmc,
      Rcpp::Named("lambda") = lambda_mcmc,
  ));
}




