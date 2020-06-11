// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h> // to use sparse matrix
#include <tgmath.h>
using namespace Rcpp;
using namespace arma;


// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
#include "helper.h"

/*
 * We would like to develope a Simulteneous Regressive Graphical LASSO, 
 * Basic idea was to embed Graphical LASSO into a normal LASSO using the hirechical structure
 *   described by Wang (2012) and Park and Casella 2008
 * 
 * 
 */



/*Full conditional dist'n of beta was normal with:
 *   
 * Sigma_beta = (\sum_i [X_i]^T\Omega[X_i] + D_{\tau}^{-1} )^{-1}
 * 
 * mu_beta = (\sum_i [X_i]^T\Omega[X_i] + D_{\tau}^{-1} )^{-1} (\sum_i [X_i]^T\Omega[(Z_i-\mu)]  )
 * 
 */

/*
 * Test in 0610, we assume Sigma * beta_{i,*}|Sigma,tau_{i,*} \sim N(0,D_{i})
 * 
 */

// tested for dimension compatibility 20200603
arma::mat update_beta_helper(const arma::mat & data,
                             const arma::mat & design,
                             const arma::vec & mu,
                             const arma::vec & tau2,
                             const arma::mat & Omega,
                             int k, int p, int n){
  
  arma::mat Y = data.t(); // convert to col vectors
  arma::vec mu_beta(k*p,fill::zeros);
  arma::mat Q_beta(k*p,k*p,fill::zeros);// percision matrix up to sigma^2 scaling
  //Q_beta.diag() += 1/tau2;
  arma::mat D_i(k,k,fill::zeros);
  arma::mat res;
  
  arma::uvec ind_para = linspace<uvec>(0,k-1,k);
  
  arma::mat chol_Omega = chol(Omega);
  //arma::mat sqrt_Omega =  sqrtmat_sympd(Omega);
  
  for(int i = 0 ; i < p ; ++i){
    D_i.zeros();
    D_i.diag() = 1/tau2(ind_para * p + i);
    //D_i = sqrt_Omega * D_i * sqrt_Omega;
    D_i = chol_Omega.t() * D_i * chol_Omega;
    Q_beta(ind_para * p + i,ind_para * p + i) = D_i;
  }
  
  arma::sp_mat X_i;
  
  for(int i = 0 ; i < n ; ++i){
    X_i = getDesign_i_helper(design.row(i),k);// this function was in helper.cpp
    Q_beta +=  X_i.t() * Omega * X_i ;
    mu_beta +=  X_i.t() * Omega * (Y.col(i)-mu);
  }
  
  
  //Rcout << "beta" <<endl;
  arma::mat Sigma_beta = inv_sympd(Q_beta);
  mu_beta = Sigma_beta*mu_beta;
  res = mvnrnd(mu_beta, Sigma_beta);
  res = reshape(res,p,k);
  
  return(res);
}

// tested dimension 20200603
arma::vec update_mu_helper(const arma::mat & data,
                           const arma::mat & design,
                           const arma::mat & beta,
                           const arma::mat & Omega, 
                           int k,int p,int n){
  arma::vec mu(k);
  arma::mat Sigma_mu = inv_sympd(Omega)/n;
  arma::mat YminusXbeta = data - design * beta;
  arma::vec mu_mu = trans(mean(YminusXbeta));
  
  //Rcout << "mu" <<endl;
  mu = mvnrnd(mu_mu, Sigma_mu);
  
  return(mu);
  
}


// return Omega matrix
// tested dimension 20200603
arma::mat update_Omega_helper(const arma::mat & data,
                              const arma::mat & design,
                              const arma::vec & mu,
                              const arma::mat & beta,
                              const double & lambda_curr,
                              int k, int p,int n){
  //arma::mat Omega;
  arma::mat Y_tilde;
  
  arma::mat expectation = design * beta;
  expectation.each_row() += mu.t();
  
  Y_tilde = data - expectation;
  
  arma::mat S = Y_tilde.t() * Y_tilde;
  arma::mat Sigma = cov(Y_tilde);
  
  //Rcout << "det cov of centered data: " << det(Sigma) << endl;
  //Rcout << "lambda of Omega: " << lambda_curr <<endl;
  
  // Concentration matrix and it's dimension:
  arma::mat Omega = pinv(Sigma); // Moore-Penrose inverse
  //int d = Omega.n_rows;
  arma::uvec pertub_vec = linspace<uvec>(0,k-1,k); 
  
  arma::uvec Omega_upper_tri = trimatu_ind(size(Omega),1);
  int n_upper_tri = Omega_upper_tri.n_elem;
  
  arma::mat tau_curr(k,k,fill::zeros);
  
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
    
    Sigma11 = Sigma(perms_j,perms_j);
    Sigma12 = Sigma.col(j);
    Sigma12 = Sigma12(perms_j);
    
    S21 = S.row(j);
    S21 = S21(perms_j);
    
    Omega11inv = Sigma11-Sigma12 * Sigma12.t() / Sigma(j,j);
    Ci = (S(j,j)+lambda_curr) * Omega11inv;
    Ci.diag() += 1/tauI;
    CiChol = chol(Ci);
    
    S_temp = S.col(j);
    S_temp = S_temp(perms_j);
    mui = solve(-Ci,S_temp);
    
    gamma = mui+solve(CiChol,randn(k-1));
    
    // Replacing omega entries
    Omega.submat(perms_j,ind_j) = gamma;
    Omega.submat(ind_j,perms_j) = gamma.t();
    
    
    
    gamm_rn = R::rgamma(n/2+1,2/( as_scalar( S(0,0) )+lambda_curr));
    Omega(j,j) = gamm_rn + as_scalar( gamma.t() * Omega11inv * gamma);
    
    
    // Replacing sigma entries
    OmegaInvTemp = Omega11inv * gamma;
    Sigma.submat(perms_j,perms_j) = Omega11inv+(OmegaInvTemp * OmegaInvTemp.t())/gamm_rn ;
    
    //flagged
    
    Sigma(perms_j,ind_j) = -OmegaInvTemp/gamm_rn;
    
    
    Sigma(ind_j,perms_j) = trans( -OmegaInvTemp/gamm_rn);
    Sigma(j,j) = 1/gamm_rn ;
    
    
    
  }
  
  return(Omega);
}


arma::vec update_tau2_helper(const arma::mat & beta,
                             const double & lambda2,
                             const arma::mat & Omega,
                             int k, int p, int n){
  //arma::mat sqrt_Omega =  sqrtmat_sympd(Omega);
  arma::mat chol_Omega = chol(Omega);
  arma::mat beta_temp = beta * chol_Omega;
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


// [[Rcpp::export]]
Rcpp::List Sigma_to_CAR_Cpp(const arma::mat & Sigma){
  
  arma::mat Q = inv_sympd(Sigma); // percision matrix
  arma::mat R = -Sigma; // R matrix, see Ver Hoef et al. 2018
  R.diag() *= 0;
  arma::vec M = 1/Q.diag();
  arma::mat invD(size(Q));
  invD.diag() = M ;
  return(Rcpp::List::create(
      Rcpp::Named("C") = invD * R, // C matrix per Ver Hoef et al. 2018
      Rcpp::Named("M") = M // diag of M matrix per Ver Hoef et al. 2018
  ));
}



