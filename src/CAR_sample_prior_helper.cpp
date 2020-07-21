// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h> 
#include <tgmath.h>
using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
#include "helper.h"


void sample_Omega_prior_helper_naive(arma::mat & Omega,
                         const double & lambda_curr,
                         int k) {
  arma::uvec pertub_vec = linspace<uvec>(0,k-1,k); 
  
  arma::uvec Omega_upper_tri = trimatu_ind(size(Omega),1);
  int n_upper_tri = Omega_upper_tri.n_elem;
  
  arma::mat tau_curr(k,k,fill::zeros);
  
  arma::uvec perms_j;
  arma::uvec ind_j(1,fill::zeros);
  arma::vec tauI;
  arma::mat Omega11;
  arma::mat Omega11inv;
  arma::mat Omega12;
  
  arma::mat S12(k-1,1,fill::zeros);
  
  arma::mat Ci;
  arma::mat invCi;
  
  arma::mat CiChol;

  arma::mat mui;
  arma::mat gamma;
  double gamm_rn = 0;
  arma::mat OmegaInvTemp;
  
  
  tau_curr.zeros();
  for(int j = 0 ; j < n_upper_tri ; ++j){
    tau_curr(Omega_upper_tri(j)) = 
      rinvGau(sqrt(lambda_curr*lambda_curr/(Omega(Omega_upper_tri(j))*Omega(Omega_upper_tri(j)))),
              lambda_curr*lambda_curr);
  }
  
  tau_curr += tau_curr.t(); // use symmertric to update lower tri
  
  for(int j = 0 ; j < k ; ++j){
    perms_j = find(pertub_vec!=j);
    ind_j = ind_j.zeros();
    ind_j += j;
    tauI = tau_curr.col(j);
    tauI = tauI(perms_j);
    
    Omega11 = Omega(perms_j,perms_j);
    Omega12 = Omega(perms_j,ind_j);

    
    Omega11inv = inv(Omega11);
    
    Ci = (lambda_curr) * Omega11inv;
    Ci.diag() += tauI;
    invCi = inv(Ci);

    mui = -invCi*S12;
    
    gamma = mvnrnd(mui, invCi);
    
    // Replacing omega entries
    Omega.submat(perms_j,ind_j) = gamma;
    Omega.submat(ind_j,perms_j) = gamma.t();
    
    
    
    gamm_rn = R::rgamma(1,2/(lambda_curr));
    Omega(j,j) = gamm_rn + as_scalar( gamma.t() * Omega11inv * gamma);
    
  }
  
  //return(Omega);
}

/*
This function was for sampling Omega directly from the prior distribution,
which is useful when evaluating the optimal experimental design numerically, 
prior to any experiment 
*/


// [[Rcpp::export]]
Rcpp::List sample_Omega_prior_cpp(int k,
                                  const int n_iter,
                                  const int n_burn_in,
                                  const int thin_by,
                                  const double lambda_a, // Gamma prior for LASSO parameter
                                  const double lambda_b,
                                  bool progress) { 
  int n_store = floor(n_iter/thin_by);
  //arma::mat Sigma_mcmc(n_store,k*k,fill::zeros);
  //Sigma_mcmc += NA_REAL;
  
  arma::mat Omega(k,k,fill::eye);

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
  Progress prog((n_iter+n_burn_in), progress); 
  for(int i = 0 ; i < (n_iter + n_burn_in) ; ++i){
        if (Progress::check_abort()){
            Rcerr << "keyboard abort\n"; // abort using esc
            return(Rcpp::List::create(
                    //Rcpp::Named("Sigma") = Sigma_mcmc,
                    Rcpp::Named("Omega") = Omega_mcmc,
                    Rcpp::Named("lambda") = lambda_mcmc)
                   );
        }
        lambda_b_post = (lambda_b+sum(sum(abs(Omega)))/2);
        lambda_curr = R::rgamma(lambda_a_post,1/lambda_b_post);
        sample_Omega_prior_helper_naive(Omega, lambda_curr, k);
        if( (i-n_burn_in)>=0 && (i+1-n_burn_in)%thin_by == 0 ){
            lambda_mcmc(i_save) = lambda_curr;
            //Sigma_mcmc.row(i_save) = Sigma(trimatu_ind(size(Sigma)));
            Omega_mcmc.row(i_save) = trans(Omega(trimatu_ind(size(Omega))));
            i_save++ ;
        }
        prog.increment();
    }  
  return(Rcpp::List::create(
            Rcpp::Named("Omega") = Omega_mcmc,
            Rcpp::Named("lambda") = lambda_mcmc
        )
    );

}

