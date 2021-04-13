// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h> 
#include <tgmath.h>
using namespace Rcpp;
using namespace arma;
using namespace std;
#include "helper.h"
#include "GIG_helper.h"
#include "CAR_LASSO_helper.h"
#include "CAR_LASSO_randeff_helper.h"


// helper function for updating random effect matrix, so that the conditional mean was design_r * res

arma::mat update_car_nu_helper(const arma::mat & data,
                               const arma::mat & design, // design matrix
                               const arma::mat & design_r, // design mat for random effect
                               const arma::mat & membership, // membership matrix for random precision, this should be something similar to a design matrix, that tells the algorithm which precision correspond to which random effect, important for multiple memberships, then the precision will be membership * xi
                               const arma::mat & beta,
                               const arma::vec & mu, // grand mean
                               const arma::mat & xi, // the precision matrix of random effect, should be m rows k column, m is how many different distributions we have
                               const arma::mat & Omega,
                               int k, int pr, int n){

  arma::mat Q_nu(k*pr,k*pr,fill::zeros);
  
  arma::mat XtX = design_r.t() * design_r ;

  arma::mat Sigma = inv_sympd(Omega);

  arma::mat mu_nu_mat = data*Omega - design * beta;
  mu_nu_mat.each_row() -= mu.t();
  mu_nu_mat = design_r.t() * (mu_nu_mat * Sigma);

  arma::vec mu_nu = vectorise(mu_nu_mat);
  
  Q_nu = kron(Sigma,XtX); // precision matrix of beta, more intuitive way was sum_i X_i^TSigmaX_i, but kron is faster
  arma::mat perc_random = membership * xi;
  //Rcout << perc_random << endl;
  Q_nu.diag() += vectorise(perc_random);
  //Rcout << Q_nu << endl;
  //Rcout << mu_nu << endl;

  arma::mat res(size(mu_nu),fill::randn);

  arma::mat chol_Q = arma::chol(Q_nu);
  res = arma::solve(chol_Q,res) + arma::solve(Q_nu,mu_nu);

  res = reshape(res,pr,k);
  
  return(res);
  
}

// used in random effect using conjugate 
arma::mat update_srg_nu_helper(const arma::mat & data,
                               const arma::mat & design, // design matrix
                               const arma::mat & design_r, // design mat for random effect
                               const arma::mat & membership, // membership matrix for random precision, this should be something similar to a design matrix, that tells the algorithm which precision correspond to which random effect, important for multiple memberships, then the precision will be membership * xi
                               const arma::mat & beta,
                               const arma::vec & mu, // grand mean
                               const arma::mat & xi, // the precision matrix of random effect, should be m rows k column, m is how many different distributions we have
                               const arma::mat & Omega,
                               int k, int pr, int n){

  arma::mat Q_nu(k*pr,k*pr,fill::zeros);
  
  arma::mat XtX = design_r.t() * design_r ;

  //arma::mat Sigma = inv_sympd(Omega);

  arma::mat mu_nu_mat = data - design * beta;
  mu_nu_mat.each_row() -= mu.t();
  mu_nu_mat = design_r.t() * (mu_nu_mat);

  arma::vec mu_nu = vectorise(mu_nu_mat);
  
  Q_nu = kron(Omega,XtX); // precision matrix of beta, more intuitive way was sum_i X_i^TSigmaX_i, but kron is faster
  arma::mat perc_random = membership * xi;
  //Rcout << perc_random << endl;
  Q_nu.diag() += vectorise(perc_random);
  //Rcout << Q_nu << endl;
  //Rcout << mu_nu << endl;

  arma::mat res(size(mu_nu),fill::randn);

  arma::mat chol_Q = arma::chol(Q_nu);
  res = arma::solve(chol_Q,res) + arma::solve(Q_nu,mu_nu);

  res = reshape(res,pr,k);
  
  return(res);
  
}




void update_xi_helper(arma::mat & xi,
                      const arma::mat & nu,
                      const arma::mat & membership,
                      const double & alpha,
                      const double & beta, 
                      int k, int pr, int m){
    double alpha_post; 
    double beta_post;
    int nn; // number of latent variable for that member
    arma::vec latent_m; //vector of latent variable of mth membership
    for(int i = 0; i < m; ++i){
        for (int j = 0; j < k; ++j)
        {
            nn = as_scalar(sum(membership.col(i))); // number of latent var for node i, member j 
            alpha_post = alpha + nn/2;
            latent_m = membership.col(i) % nu.col(j);// this vector was 0 for non members and nu for members
            beta_post = beta + arma::as_scalar( arma::sum(  latent_m % latent_m) )/2;
            xi(i,j) = randg<double>( distr_param(alpha_post,1/beta_post) );
        }
    }
    return;
}


void get_data_centered(arma::mat & centered_data,
                       const arma::mat & data,
                       const arma::mat & design_r,
                       const arma::mat & nu,
                       const arma::mat & Omega){
    //arma::mat Sigma = inv(Omega);
    centered_data = data - arma::trans( arma::solve(Omega, arma::trans(design_r * nu))); //* Sigma;
    return;
}



void update_car_randeff_Omega_helper(arma::mat & Omega,
                             const arma::mat & data,
                             const arma::mat & design,
                             const arma::mat & design_r,
                             const arma::mat & nu,
                             const arma::vec & mu,
                             const arma::mat & beta,
                             const double & lambda_curr,
                             int k, int p,int n){
  //arma::mat Omega;
  //arma::mat Y_tilde;
  
  arma::mat expectation = design * beta + design_r * nu;
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
  //double S22;
  
  arma::mat U11;
  arma::mat U12;
  //double U22;
  
  arma::mat Omega_11;
  arma::mat inv_Omega_11;
  
  arma::mat omega_12;
  arma::mat Sigma_omega_12(k-1,k-1,fill::zeros);
  arma::mat chol_Omega_omega_12;
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
    //inv_Omega_11 += inv_Omega_11.t();
    //inv_Omega_11/=2;
    // the current gamma=Omega_22-omega_12^T * Omega_{11}^{-1} * omega_{12}
    gamma = as_scalar( Omega(j,j)-Omega(just_j,perms_j)*inv_Omega_11*Omega(perms_j,just_j));
    
    // update omega_12, which is normal
    Omega_omega_12 = (S(j,j)+lambda_curr)*inv_Omega_11+(1/gamma)*inv_Omega_11*U11*inv_Omega_11;
    Omega_omega_12.diag() += tauI;
    // substitute inv-mvnrnd to chol
    //Sigma_omega_12 = inv(Omega_omega_12);
    
    //mu_omega_12 = (S12-(1/gamma) * inv_Omega_11*U12);
    //mu_omega_12 = - Sigma_omega_12 * mu_omega_12;
    
    //omega_12 = mvnrnd(mu_omega_12, Sigma_omega_12);
    chol_Omega_omega_12 = chol(Omega_omega_12);
    mu_omega_12 = (S12-(1/gamma) * solve(Omega_11,U12));
    mu_omega_12 = - solve(Omega_omega_12,mu_omega_12);
    omega_12 = solve(chol_Omega_omega_12,randn(size(mu_omega_12))) + mu_omega_12;

    
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
