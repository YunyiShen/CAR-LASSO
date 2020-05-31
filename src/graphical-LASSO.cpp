// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h> // to use sparse matrix
#include <tgmath.h>
using namespace Rcpp;
static const double pi = 3.141592653589793238462643383280;

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

// [[Rcpp::export]]
double rinvGau_trunc(const double & mu, const double & lambda){
  mu = mu < 1e-12 ? 1e-12 : mu;
  double b = 0.5 * mu / lambda;
  double a = mu * b;
  double c = 4.0 * mu * lambda;
  double d = mu * mu;
  double res;
  double u = R::runif(0,1);
  double chisqr = R::rchisq(1);
  
  double x = mu + a * chisqr - b * sqrt( c * chisqr + d * chisqr * chisqr); // solve the smaller root
  res = (u < (mu / (mu + x))) ? x : d/x; // accept the small one with prob mu/(mu+x), otherwise the larger one
  return(res);
}



/*
 * Implements the block Gibbs sampler for the Bayesian graphical lasso
 * introduced in Wang (2012). Samples from the conditional distribution of a 
 * permuted column/row for simulating the posterior distribution for the concentration 
 * matrix specifying a Gaussian Graphical Model
 */

Rcpp::List Graphical_LASSO_Cpp(const arma::mat & data,
                           const int nIter,
                           const int burn_in,
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
  int p = Omega.n_rows;
  
  
  // indicator matrix and permutation matrix for looping through columns & rows ("blocks")
  //arma::vec indMat_vec = linspace<vec>(1, p*p , p*p);
  //arma::mat indMat = reshape(indMat_vec,p,p);
  arma::uvec pertub_vec = linspace<uvec>(0,p-1,p); // we will not convert this to matrix
  
  // mcmc storage
  n_store = floor(nIter/thin_by);
  arma::mat Sigma_mcmc(n_store,p*p,fill::zeros);
  Sigma_mcmc += NA_REAL;
  
  arma::mat Omega_mcmc(n_store,p*p,fill::zeros);
  Omega_mcmc += NA_REAL;
  
  arma::vec lambda_mcmc(n_store,fill::zeors);
  lambda_mcmc + NA_REAL;
  
  int i_save = 0;
  
  // latent tau
  arma::mat tau_curr(p,p,fill::zeros);
  
  double lambda_a_post = (lambda_a+(p*(p+1)/2));
  double lambda_b_post;
  
  double lambda_curr = 0;
  
  arma::uvec Omega_upper_tri = trimatu_ind(size(Omega),1);
  int n_upper_tri = Omega_upper_tri.n_elem;
  
  // some objects needed in block update
  arma::uvec perms_j;
  arma::uvec ind_j(1,fill::zeros);
  arma::vec tauI;
  arma::mat Sigma11;
  arma::mat Sigma12;
  arma::mat S21;
  arma::mat Omega11inv;
  arma::mat Ci;
  arma::mat CiChol;
  arma::mat S_temp;
  arma::mat mui;
  arma::mat beta;
  double gamm;
  arma::mat OmegaInvTemp;
  
  // progress bar
  Progress p((n_iter+n_burn_in), progress); 
  
  // main iterations
  for(int i = 0 ; i < (nIter+burn_in) ; ++i){
    if (Progress::check_abort()){
      Rcerr << "keyboard abort\n"; // abort using esc
      return(Rcpp::List::create(
          Rcpp::Named("Sigma") = Sigma_mcmc,
          Rcpp::Named("Omega") = Omega_mcmc,
          Rcpp::Named("lambda") = lambda_mcmc);
    }
    
    // update LASSO parameter lambda
    lambda_b_post = (lambda_b+sum(sum(abs(Omega)))/2);
    lambda_curr = R::rgamma(lambda_a_post,1/lambda_b_post);
    
    // update latent tau (it is symmertric so we only work on upper tri)
    tau_curr.zeros();
    for(int j = 0 ; j < n_upper_tri ; ++n){
      tau_curr(Omega_upper_tri(j)) = 
        rinvGau_trunc(abs(lambda_curr/Omega(Omega_upper_tri(j))),
                      lambda_curr*lambda_curr)
    }
    tau_curr = tau_curr + t_curr.t(); // use symmertric to update lower tri
    
    for(int j = 0 ; j < p ; ++j){
      perms_j = find(pertub_vec!=j);
      ind_j = ind_j.zeros();
      ind_j += j;
      tauI = tau_curr.col(j);
      tauI = tauI(perms_j);
      
      Sigma11 = Sigma(perms_j,perms_j);
      Sigma12 = Sigma.col(j);
      Sigma12 = Sigma12(perms_j);
      
      S21 = S.row(i);
      S21 = S21(perms_j);
      
      Omega11inv = Sigma11-Sigma12 * Sigma12.t() / Sigma(j,j);
      Ci = (S(j,j)+lambda_curr) * Omega11inv;
      Ci.diag() += 1/tauI;
      CiChol = chol(Ci)
      
      S_temp = S.col(i);
      S_temp = S_temp(perms_j);
      mui = solve(-Ci,S_temp);
      
      // Sampling:
      beta = mui+solve(CiChol,randn(p-1));
        
      // Replacing omega entries
      Omega.submat(perms_j,ind_j) = beta;
      Omega.submat(ind_j,perms_j) = beta;
      gamm = R::rgamma(n/2+1,(S(0,0)+lambda_curr)/2);
      Omega(j,j) = gamm + beta.t() * Omega11inv * beta;
          
      // Replacing sigma entries
      OmegaInvTemp = Omega11inv * beta;
      Sigma.submat(perms_j,perms_j) = Omega11inv+(OmegaInvTemp * OmegaInvTemp.t())/gamm ;
      Sigma[perms_j,ind_j]<-Sigma[ind_j,perms_j] = -OmegaInvTemp/gamm;
      Sigma(j,j) = 1/gamm ;
      
    }
    
    // saving current state
    if( (i-burn_in)>=0 && (i+1-burn_in)%thin_by == 0 ){
      lambda_mcmc(i_save) = lambda_curr;
      Sigma_mcmc.row(i_save) = vectorise(Sigma,1);
      Omega_mcmc.row(i_save) = vectorise(Omega,1);
      i_save++ ;
    }
    
    p.increment();
  }
  
  return(Rcpp::List::create(
      Rcpp::Named("Sigma") = Sigma_mcmc,
      Rcpp::Named("Omega") = Omega_mcmc,
      Rcpp::Named("lambda") = lambda_mcmc);
  
}

