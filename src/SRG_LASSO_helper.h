#ifndef SRG_HELPER_H
#define SRG_HELPER_H
arma::mat update_beta_helper(const arma::mat & data,
                             const arma::mat & design,
                             const arma::vec & mu,
                             const arma::vec & tau2,
                             const arma::mat & Omega,
                             int k, int p, int n);
                             
arma::vec update_mu_helper(const arma::mat & data,
                           const arma::mat & design,
                           const arma::mat & beta,
                           const arma::mat & Omega, 
                           int k,int p,int n);

arma::mat update_Omega_helper(const arma::mat & data,
                              const arma::mat & design,
                              const arma::vec & mu,
                              const arma::mat & beta,
                              const double & lambda_curr,
                              int k, int p,int n);

arma::vec update_tau2_helper(const arma::mat & beta,
                             const double & lambda2,
                             const arma::mat & Omega,
                             int k, int p, int n);

Rcpp::List Sigma_to_CAR_Cpp(const arma::mat & Sigma);


#endif