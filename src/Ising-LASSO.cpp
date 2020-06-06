// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include <RcppArmadillo.h> // to use sparse matrix
#include <tgmath.h>
using namespace Rcpp;
using namespace arma;

#include <progress.hpp>
#include <progress_bar.hpp>
#include "helper.h"
#include "Ising_helper.h"
//#include "Ising_helper_sparse.h"

// log MH ratio for inverse Ising 
double logMH_Ising(const arma::mat & data,
                   const arma::mat & thr_prop,
                   const arma::mat & thr_curr,
                   const arma::mat & J_prop,
                   const arma::mat & J_curr,
                   const arma::vec & responses,
                   double log_prior_beta_curr,
                   double log_prior_J_curr,
                   double log_prior_beta_prop,
                   double log_prior_J_prop,             
                   bool exact,
                   int nIter,
                   int n,int k){
  
  //arma::vec Z_prime;
  arma::vec Z_prime_i;
  double log_q_theta_curr_data = 0;
  double log_q_theta_prop_data = 0;
  double log_q_theta_curr_Z_prime = 0;
  double log_q_theta_prop_Z_prime = 0;
  
  //arma::vec thrvec_prop = vectorise(thr_prop);
  //arma::sp_mat bigJ_prop = block_diag(J_prop,k,n);
  //Z_prime = IsingSampler_sparse_Cpp(1,bigJ_prop,thrvec_prop,nIter,responses,exact);
  //Rcout<<Z_prime.n_elem << endl;
  for(int i = 0 ; i < n ; ++i){
    // Optimization todo: change this to joint sample at once instead of run CFTP for n times
    //Z_prime_i = Z_prime(span(i*k,(i+1)*k-1));
    Z_prime_i = IsingSamplerCpp(1,J_prop,thr_prop.col(i),nIter,responses,exact);
    
    log_q_theta_curr_data -= H(J_curr,data.col(i),thr_curr.col(i));
    log_q_theta_prop_data -= H(J_prop,data.col(i),thr_prop.col(i));
    
    log_q_theta_curr_Z_prime -= H(J_curr,Z_prime_i,thr_curr.col(i));
    log_q_theta_prop_Z_prime -= H(J_prop,Z_prime_i,thr_prop.col(i));
  }
  
  double log_MH = log_prior_beta_prop + 
    log_prior_J_prop + 
    log_q_theta_prop_data + 
    log_q_theta_curr_Z_prime -
    log_prior_beta_curr - 
    log_prior_J_curr - 
    log_q_theta_curr_data - 
    log_q_theta_prop_Z_prime;
  
  return(log_MH);
  
}


// main sampler

// [[Rcpp::export]]
Rcpp::List Ising_LASSO_Cpp(const arma::mat & data_,
                           const arma::mat & design,
                           const int n_iter, // how many iteractions?
                           const int n_burn_in, // burn in
                           const int thin_by, // thinning?
                           const double r_beta, // prior on lambda of beta
                           const double delta_beta,
                           const double r_J,
                           const double delta_J,
                           const arma::vec & propsd_mu,
                           const arma::vec & propsd_beta,
                           const arma::vec & propsd_J,
                           const arma::vec & propsd_lambda,// beta first then J
                           bool exact, // whether use CFTP for the sampler in sampling Z*
                           int nauxIter, // still, interactions in Ising sampling 
                           bool progress,
                           bool verbos,int reportby){
  arma::vec responses(2);
  responses(0) = min(min(data));
  responses(1) = max(max(data));
  arma::mat data = data_.t(); // convert to column vector as sample
  int k = data.n_rows; // number of nodes
  int p = design.n_cols; //number of predictors
  int n = data.n_cols; // number of samples
  
  int accept_Ising = 0;
  double accept_lambda = 0;
  
  int n_save = floor(n_iter/thin_by); //
  int i_save = 0;  
  
  // mcmc matrices:
  arma::mat beta_mcmc(n_save,k * p); // beta mcmc
  beta_mcmc += NA_REAL; // initial with NAs
  
  arma::mat J_mcmc(n_save,floor(.5*(k-1)*k));// adjacency matrix
  J_mcmc += NA_REAL;
  
  arma::mat mu_mcmc(n_save , k); // mean for node 1 to k
  mu_mcmc += NA_REAL;
  
  arma::mat lambda_mcmc(n_save , 2); // LASSO parameter for beta and J
  lambda_mcmc += NA_REAL;
  
  // some initial values
  arma::mat J_mat_curr(k,k,fill::zeros);
  arma::uvec upperdiag = trimatu_ind(size(J_mat_curr),1);
  
  arma::vec beta_curr(k*p,fill::zeros);
  arma::mat beta_mat_curr = reshape(beta_curr,p,k);
  arma::vec mu_curr(k,fill::zeros);
  arma::vec J_curr(upperdiag.n_elem,fill::zeros);
  
  
  // obejects for proposal state
  arma::vec beta_prop(k*p,fill::zeros);
  arma::mat beta_mat_prop(p,k,fill::zeros);
  arma::vec J_prop(J_curr.n_elem);
  arma::mat J_mat_prop(k,k,fill::zeros);
  arma::vec mu_prop(k,fill::zeros);
  
  double lambda_beta_curr = 1;
  double lambda_J_curr = 1;
  
  double lambda_beta_prop;
  double lambda_J_prop;
  
  double log_prior_beta_curr = sum( dLaplace_Cpp(beta_curr,0,lambda_beta_curr,true));
  double log_prior_J_curr = sum( dLaplace_Cpp(J_curr,0,lambda_beta_curr,true));
  
  double log_prior_beta_prop;
  double log_prior_J_prop;
  
  double logMH;
  
  
  double logPostlambda_beta_curr = sum(dLaplace_Cpp(beta_curr,0,lambda_beta_curr,true)) + 
    R::dgamma(lambda_beta_curr,r_beta,1/delta_beta,true);
  double logPostlambda_J_curr = sum(dLaplace_Cpp(J_curr,0,lambda_J_curr,true)) + 
    R::dgamma(lambda_J_curr,r_J,1/delta_J,true);
  
  double logPostlambda_beta_prop;
  double logPostlambda_J_prop;
  
  
  arma::mat thr_curr(size(data),fill::zeros);
  arma::mat thr_prop = thr_curr;
  
  Progress prog((n_iter+n_burn_in), progress); 
  // main loop
  
  
  for(int i = 0 ; i < (n_iter+n_burn_in) ; ++i){
    
    if (Progress::check_abort()){
      Rcerr << "keyboard abort\n";
      return(Rcpp::List::create(
          Rcpp::Named("beta") = beta_mcmc,
          Rcpp::Named("mu") = mu_mcmc,
          Rcpp::Named("J") = J_mcmc,
          Rcpp::Named("lambda") = lambda_mcmc
      ));
    }
    
    // update Ising parameters
    beta_prop = propsd_beta % randn(size(beta_curr)) + beta_curr;
    beta_mat_prop = reshape(beta_prop,p,k);
    J_prop = propsd_J % randn(size(J_curr)) + J_curr;
    J_mat_prop = J_mat_prop.zeros();
    
    J_mat_prop(upperdiag) = J_prop;
    
    J_mat_prop += J_mat_prop.t();
    mu_prop = propsd_mu % randn(size(mu_curr)) + mu_curr;
    
    thr_prop = trans(design * beta_mat_prop);
    thr_prop.each_col() += mu_prop;
    
    log_prior_beta_prop = sum( dLaplace_Cpp(beta_prop,0,lambda_beta_curr,true));
    log_prior_J_prop = sum( dLaplace_Cpp(J_prop,0,lambda_J_curr,true));
    
    // MH ratio
    
    logMH = logMH_Ising(data,
                        thr_prop,thr_curr,
                        J_mat_prop,J_mat_curr,
                        responses,
                        log_prior_beta_curr,log_prior_J_curr,
                        log_prior_beta_prop,log_prior_J_prop,             
                        exact,nauxIter,n,k);
    //Rcout<<logMH<<endl;
    if(log(R::runif(0,1)) <= logMH){ // accept
      beta_curr = beta_prop;
      beta_mat_curr = beta_mat_prop;
      
      J_curr = J_prop;
      J_mat_curr = J_mat_prop;
      
      mu_curr = mu_prop;
      
      thr_curr = thr_prop;
      
      log_prior_beta_curr = log_prior_beta_prop;
      
      log_prior_J_curr = log_prior_J_prop;
      
      accept_Ising++;
      
    }
    // updating lambdas
    lambda_beta_prop = lambda_beta_curr + propsd_lambda(0) * R::rnorm(0,1);
    lambda_J_prop = lambda_J_curr + propsd_lambda(1) * R::rnorm(0,1);
    
    if(lambda_beta_prop>0){
      logPostlambda_beta_prop = sum(dLaplace_Cpp(beta_curr,0,lambda_beta_prop,true)) + 
        R::dgamma(lambda_beta_prop,r_beta,1/delta_beta,true);
      if(log(R::runif(0,1)) <= logPostlambda_beta_prop-logPostlambda_beta_curr){
        lambda_beta_curr = lambda_beta_prop;
        logPostlambda_beta_curr = logPostlambda_beta_prop;
        accept_lambda += .5 ;
      }
    }
    
    if(lambda_J_prop>0){
      logPostlambda_J_prop = sum(dLaplace_Cpp(J_curr,0,lambda_J_prop,true)) + 
        R::dgamma(lambda_J_prop,r_J,1/delta_J,true);
      if(log(R::runif(0,1)) <= logPostlambda_J_prop-logPostlambda_J_curr){
        lambda_J_curr = lambda_J_prop;
        logPostlambda_J_curr = logPostlambda_J_prop;
        accept_lambda += .5 ;
      }
    }
    
    if(verbos && ((i+1) % reportby ==0)){
      Rcout << "Iteration: " << (i+1-reportby) << " to " << "iteration " << i+1 << "\n" << endl;
      Rcout << "    Grand acceptance rate for Ising:" << 100.0 * accept_Ising/(i+1.0) << "%" <<endl;
      Rcout << "    Grand acceptance rate for LASSO:" << 100.0 * accept_lambda/(i+1.0) << "%\n" <<endl;
    }
    
    // 
    
    if( (i-n_burn_in)>=0 && (i+1-n_burn_in)%thin_by == 0 ){
      lambda_mcmc(i_save,0) = lambda_beta_curr;
      lambda_mcmc(i_save,1) = lambda_J_curr;
      
      beta_mcmc.row(i_save) = beta_curr.t();
      J_mcmc.row(i_save) = J_curr.t();
      
      mu_mcmc.row(i_save) = mu_curr.t();

      i_save++ ;
    }
    
    
    
    prog.increment(); // progress bar
  }
  
  return(Rcpp::List::create(
      Rcpp::Named("beta") = beta_mcmc,
      Rcpp::Named("mu") = mu_mcmc,
      Rcpp::Named("J") = J_mcmc,
      Rcpp::Named("lambda") = lambda_mcmc
  ));
  
  
}
