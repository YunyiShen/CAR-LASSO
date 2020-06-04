// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h> // to use sparse matrix
#include <tgmath.h>
using namespace Rcpp;
using namespace arma;
static const double pi = 3.141592653589793238462643383280;

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


/*
 * TODO: 
 *  test Gibbs sampler for SAR-lasso
 */


/*Full conditional dist'n of beta was normal with:
 *   
 * Sigma_beta = (\sum_i [X_i]^T\Omega[X_i] + D_{\tau}^{-1} )^{-1}
 * 
 * mu_beta = (\sum_i [X_i]^T\Omega[X_i] + D_{\tau}^{-1} )^{-1} (\sum_i [X_i]^T\Omega[(Z_i-\mu)]  )
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
  Q_beta.diag() += 1/tau2;
  arma::mat res;
  
  arma::sp_mat X_i;
  
  for(int i = 0 ; i < n ; ++i){
    X_i = getDesign_i_helper(design.row(i),k);
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
                              const double lambda_curr,
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
  int d = Omega.n_rows;
  arma::uvec pertub_vec = linspace<uvec>(0,d-1,d); 
  
  arma::uvec Omega_upper_tri = trimatu_ind(size(Omega),1);
  int n_upper_tri = Omega_upper_tri.n_elem;
  
  arma::mat tau_curr(d,d,fill::zeros);
  
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
  
  for(int j = 0 ; j < d ; ++j){
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
    
    gamma = mui+solve(CiChol,randn(d-1));
    
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
  arma::vec betavec = vectorise(beta);
  arma::vec invtau2(k*p);
  
  double detSigma = 1/det(Omega);
  //Rcout << "detOmega:" << 1/detSigma << endl;
  detSigma = pow(detSigma,1/(k));
  
  
  
  arma::vec mu_prime = sqrt(lambda2*detSigma/(betavec%betavec));
  for(int i = 0 ; i < k*p ; ++i){
    invtau2(i) =  rinvGau(mu_prime(i),lambda2);
  }
  return(1/invtau2);
}
/* 
 * Convention:
 * mcmc obects were vectorized column first, e.g. each row: looks like
 * beta_11,beta_12,beta_13,...,beta_kp
 * 
 */

// [[Rcpp::export]]
List SRG_LASSO_Cpp(const arma::mat & data, // col composition data, ROW as a sample
                   const arma::mat & design, // design matrix, each ROW as a sample
                   const int n_iter, // how many iterations?
                   const int n_burn_in, // burn in
                   const int thin_by, // thinning?
                   const double r_beta, // prior on lambda of beta
                   const double delta_beta,
                   const double r_Omega,
                   const double delta_Omega,
                   bool progress){
  int k = data.n_cols; // number of nodes
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
  
  arma::vec tau2_curr = randg<arma::vec> (k*p,distr_param(r_beta,delta_beta)); // current latent variable tau^2, for prior of beta

  arma::vec mu_curr = trans( mean(data) ); // current value of mean
  arma::mat centered_data = data;
  centered_data.each_row() -= mu_curr.t();
  arma::mat Omega_curr(k,k); // current value of Omega
  Omega_curr = pinv(cov(data));
  arma::mat beta_curr = solve( design.t()*design,design.t()*(centered_data)); // current value of beta
  
  
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
          Rcpp::Named("lambda") = lambda_mcmc
      ));
    }
    // block update start:
    Omega_delta_post = (delta_Omega+sum(sum(abs(Omega_curr)))/2);
    lambda_Omega = R::rgamma(Omega_r_post,1/Omega_delta_post);
    
    //Update betas:
    beta_curr = update_beta_helper(data,design,mu_curr,
                                   tau2_curr,Omega_curr, 
                                   k,p,n);
    
    
    // update Omega
    Omega_curr = update_Omega_helper(data, design, 
                                     mu_curr,beta_curr,
                                     lambda_Omega,
                                     k,p,n);
    
    
    
    // Update mu
    
    mu_curr = update_mu_helper(data,design,beta_curr,
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

