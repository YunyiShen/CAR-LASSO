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
#include "GIG_helper.h"
#include "CAR_LASSO_helper.h"
/*
 * Helper functions for Conditional Auto Regression LASSO, 
 * Basic idea was to embed Graphical LASSO into a normal LASSO using the hirechical structure
 *   described by Wang (2012) and Park and Casella 2008
 * In this model average structure offer some extra information of conditional correlation
 * 
 * A CAR can be reparameterize into a model s.t.: 
 * Y~N(Sigma (Xbeta+mu),Sigma)
 */


// [[Rcpp::export]]
arma::mat update_car_beta_helper(const arma::mat & data,
                             const arma::mat & design,
                             const arma::vec & mu,
                             const arma::vec & tau2,
                             const arma::mat & Omega,
                             int k, int p, int n){
  

  arma::mat Q_beta(k*p,k*p,fill::zeros);
  
  arma::mat XtX = design.t() * design;
  
  
  arma::mat res;
  
  arma::uvec ind_para = linspace<uvec>(0,k-1,k);
  
  arma::mat Sigma = inv_sympd(Omega);

  arma::mat Y_tilde = data;
  Y_tilde.each_row() -= mu.t() * Sigma;
  arma::uvec ind_p = linspace<uvec>(0,p-1,p);
  
  
  arma::mat mu_beta_mat = design.t() * Y_tilde;
  arma::vec mu_beta = vectorise(mu_beta_mat);
  
  Q_beta = kron(Sigma,XtX); // precision matrix of beta, more intuitive way was sum_i X_i^TSigmaX_i, but kron is faster

  Q_beta.diag() += 1/tau2; // penalty due to Laplace prior
  
  
  arma::mat Sigma_beta = inv(Q_beta);
  
  mu_beta = Sigma_beta*mu_beta;
  res = mvnrnd(mu_beta, Sigma_beta);
  res = reshape(res,p,k);
  
  return(res);
}


// tested 20200712
// [[Rcpp::export]]
arma::vec update_car_mu_helper(const arma::mat & data,
                           const arma::mat & design,
                           const arma::mat & beta,
                           const arma::mat & Omega, 
                           int k,int p,int n){
  arma::vec mu(k);
  
  // mu was normal with mean Omega*data-Xbeta
  arma::mat YminusXbeta = data*Omega - design * beta;
  arma::vec mu_mu = trans(mean(YminusXbeta));
  
  
  mu = mvnrnd(mu_mu, Omega/n);
  
  return(mu);
  
}


// update Omega matrix using block method
// [[Rcpp::export]]
void update_car_Omega_helper(arma::mat & Omega,
                             const arma::mat & data,
                             const arma::mat & design,
                             const arma::vec & mu,
                             const arma::mat & beta,
                             const double & lambda_curr,
                             int k, int p,int n){
  //arma::mat Omega;
  //arma::mat Y_tilde;
  
  arma::mat expectation = design * beta;
  expectation.each_row() += mu.t();
  
  
  arma::mat S = data.t() * data;
  arma::mat U = expectation.t() * expectation;
  
  arma::uvec pertub_vec = linspace<uvec>(0,k-1,k); 
  
  arma::uvec Omega_upper_tri = trimatu_ind(size(Omega),1);
  int n_upper_tri = Omega_upper_tri.n_elem;
  
  arma::mat tau_curr(k,k,fill::zeros);
  
  arma::uvec perms_j;
  arma::uvec just_j;
  arma::uvec ind_j(1,fill::zeros);
  arma::vec tauI;
  
  arma::mat S11;
  arma::mat S12;
  double S22;
  
  arma::mat U11;
  arma::mat U12;
  double U22;
  
  arma::mat Omega_11;
  arma::mat inv_Omega_11;
  
  arma::mat omega_12;
  arma::mat Sigma_omega_12(k-1,k-1,fill::zeros);
  arma::mat Omega_omega_12(k-1,k-1,fill::zeros);
  arma::mat mu_omega_12;
  
  double gamma;
  double lambda_gamma,psi_gamma,chi_gamma;
  
  
  
  tau_curr.zeros();
  // update tau, using old Omega
  for(int j = 0 ; j < n_upper_tri ; ++j){
    tau_curr(Omega_upper_tri(j)) = 
      rinvGau(sqrt(lambda_curr*lambda_curr/(Omega(Omega_upper_tri(j))*Omega(Omega_upper_tri(j)))),
              lambda_curr*lambda_curr);
  }
  
  tau_curr += tau_curr.t(); // use symmertric to update lower tri
  
  for(int j = 0 ; j < k ; ++j){
    perms_j = find(pertub_vec!=j);
    just_j = find(pertub_vec==j);
    tauI = tau_curr(perms_j,just_j); // tau for this column
    
    
    // partitioning:
    S11 = S(perms_j,perms_j);
    S12 = S(perms_j,just_j);
    U11 = U(perms_j,perms_j);
    U12 = U(perms_j,just_j);
    Omega_11 = Omega(perms_j,perms_j);

    inv_Omega_11 = inv(Omega_11);
    
    // the current gamma=Omega_22-omega_12^T * Omega_{11}^{-1} * omega_{12}
    gamma = as_scalar( Omega(j,j)-Omega(just_j,perms_j)*inv_Omega_11*Omega(perms_j,just_j));
    
    // update omega_12, which is normal
    Omega_omega_12 = (S(j,j)+lambda_curr)*inv_Omega_11+(1/gamma)*inv_Omega_11*U11*inv_Omega_11;
    Omega_omega_12.diag() += tauI;

    Sigma_omega_12 = inv(Omega_omega_12);
    
    mu_omega_12 = (S12-(1/gamma) * inv_Omega_11*U12);
    mu_omega_12 = - Sigma_omega_12 * mu_omega_12;
    
    omega_12 = mvnrnd(mu_omega_12, Sigma_omega_12);
    
    Omega(perms_j,just_j) = omega_12;
    Omega(just_j,perms_j) = omega_12.t();
    
    // update gamma, follow GIG 
    lambda_gamma = n/2+1;
    psi_gamma = lambda_curr + S(j,j);
    chi_gamma = U(j,j) - 
      2*as_scalar(U12.t()*inv_Omega_11*omega_12)+
      as_scalar(omega_12.t()*inv_Omega_11*U11*inv_Omega_11*omega_12);
    
    gamma = rgig(lambda_gamma, chi_gamma, psi_gamma); // function in GIG_helper.cpp
    
    // update diagonal using determinant
    Omega(j,j) = gamma + as_scalar( omega_12.t()*inv_Omega_11*omega_12);
    
  }
  return;
}

// update 1/tau^2, this is for beta, 
// tau in Omega was updated within update_Omega
// [[Rcpp::export]]
arma::vec update_car_tau2_helper(const arma::mat & beta,
                             const double & lambda2,
                             const arma::mat & Omega,
                             int k, int p, int n){
  arma::vec betavec = vectorise(beta);
  arma::vec invtau2(k*p);
  
  
  arma::vec mu_prime = sqrt(lambda2/(betavec%betavec));
  for(int i = 0 ; i < k*p ; ++i){
    invtau2(i) =  rinvGau(mu_prime(i),lambda2); 
  }
  return(1/invtau2);
}
