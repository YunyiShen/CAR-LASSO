// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h> // to use sparse matrix
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
 * We would like to develope a Conditional Auto Regression LASSO, 
 * Basic idea was to embed Graphical LASSO into a normal LASSO using the hirechical structure
 *   described by Wang (2012) and Park and Casella 2008
 * 
 * Y~N(Sigma (Xbeta+mu),Sigma)
 */



// tested for dimension compatibility 20200603
arma::mat update_car_beta_helper(const arma::mat & data,
                             const arma::mat & design,
                             const arma::vec & mu,
                             const arma::vec & tau2,
                             const arma::mat & Omega,
                             int k, int p, int n){
  
  arma::mat Y = data.t(); // convert to col vectors
  arma::vec mu_beta(k*p,fill::zeros);
  arma::mat Q_beta(k*p,k*p,fill::zeros);// percision matrix up to sigma^2 scaling
  Q_beta.diag() += 1/tau2;
  arma::mat D_i(k,k,fill::zeros);
  arma::mat res;
  
  arma::uvec ind_para = linspace<uvec>(0,k-1,k);
  
  arma::mat Sigma = inv_sympd(Omega);
  //arma::mat sqrt_Omega =  sqrtmat_sympd(Omega);
  arma::mat chol_Omega =  chol(Omega);
  //for(int i = 0 ; i < p ; ++i){
  //  D_i.zeros();
  //  D_i.diag() = 1/tau2(ind_para * p + i);
    //D_i = sqrt_Omega * D_i * sqrt_Omega;
    
  //  D_i = chol_Omega.t() * D_i * chol_Omega;
  //  Q_beta(ind_para * p + i,ind_para * p + i) = D_i;
  //}
  
  arma::sp_mat X_i;
  
  for(int i = 0 ; i < n ; ++i){
    X_i = getDesign_i_helper(design.row(i),k);// this function was in helper.cpp
    Q_beta +=  X_i.t() * Sigma * X_i ;
    mu_beta +=  X_i.t() * (Y.col(i)-Sigma * mu);
  }
  
  
  //Rcout << "beta" <<endl;
  arma::mat Sigma_beta = inv_sympd(Q_beta);
  mu_beta = Sigma_beta*mu_beta;
  res = mvnrnd(mu_beta, Sigma_beta);
  res = reshape(res,p,k);
  
  return(res);
}

// tested dimension 20200603
arma::vec update_car_mu_helper(const arma::mat & data,
                           const arma::mat & design,
                           const arma::mat & beta,
                           const arma::mat & Omega, 
                           int k,int p,int n){
  arma::vec mu(k);
  //Rcout << "mu" <<endl;
  //arma::mat Sigma_mu = inv_sympd(Omega);
  arma::mat YminusXbeta = data*Omega - design * beta;
  arma::vec mu_mu = trans(mean(YminusXbeta));
  
  
  mu = mvnrnd(mu_mu, Omega/n);
  
  return(mu);
  
}


// return Omega matrix
// tested dimension 20200603
void update_car_Omega_helper(arma::mat Omega,
                             const arma::mat & data,
                             const arma::mat & design,
                             const arma::vec & mu,
                             const arma::mat & beta,
                             const double & lambda_curr,
                             int k, int p,int n){
  //arma::mat Omega;
  arma::mat Y_tilde;
  
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
    ind_j = ind_j.zeros();
    ind_j += j;
    tauI = tau_curr(perms_j,just_j); // tau for this column
    
    S11 = S(perms_j,perms_j);
    S12 = S(perms_j,just_j);
    U11 = U(perms_j,perms_j);
    U12 = U(perms_j,just_j);
    Omega_11 = Omega(perms_j,perms_j);
    inv_Omega_11 = inv_sympd(Omega_11);
    
    // the current gamma
    gamma = as_scalar( Omega(j,j)-Omega(just_j,perms_j)*inv_Omega_11*Omega(perms_j,just_j));
    
    // update omega_12
    Omega_omega_12 = (S(j,j)+lambda_curr)*inv_Omega_11+(1/gamma)*inv_Omega_11*U11*inv_Omega_11;
    Omega_omega_12.diag() += 1/tauI;
    Sigma_omega_12 = inv_sympd(Omega_omega_12);
    //Rcout << "flag" <<endl;
    mu_omega_12 = S12+(1/gamma) * inv_Omega_11*U12;
    mu_omega_12 = Sigma_omega_12 * mu_omega_12;
    
    omega_12 = mvnrnd(mu_omega_12, Sigma_omega_12);
    
    Omega(perms_j,just_j) = omega_12;
    Omega(just_j,perms_j) = omega_12.t();
    // flagged, after error
    // update gamma
    lambda_gamma = k/2+1;
    psi_gamma = lambda_curr + S(j,j);
    chi_gamma = U(j,j) + 
      2*as_scalar(U12.t()*inv_Omega_11*omega_12)+
      as_scalar(omega_12.t()*inv_Omega_11*U11*inv_Omega_11*omega_12);
    
    gamma = rgig(lambda_gamma, chi_gamma, psi_gamma);
    
    // update diagonal
    Omega(j,j) = gamma + as_scalar( omega_12.t()*inv_Omega_11*omega_12);
    
  }
}


arma::vec update_car_tau2_helper(const arma::mat & beta,
                             const double & lambda2,
                             const arma::mat & Omega,
                             int k, int p, int n){
  //arma::mat sqrt_Omega =  sqrtmat_sympd(Omega);
  arma::mat chol_Omega = chol(Omega);
  arma::mat beta_temp = beta;
  arma::vec betavec = vectorise(beta_temp);
  arma::vec invtau2(k*p);
  
  //double detSigma = 1/det(Omega);
  //Rcout << "detOmega:" << 1/detSigma << endl;
  //detSigma = pow(detSigma,1/(k));
  
  
  
  arma::vec mu_prime = sqrt(lambda2/(betavec%betavec));
  for(int i = 0 ; i < k*p ; ++i){
    invtau2(i) =  rinvGau(mu_prime(i),lambda2);
  }
  return(1/invtau2);
}





