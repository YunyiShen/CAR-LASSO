#ifndef PROBIT_HELPER_H
#define PROBIT_HELPER_H

// version for regression
void update_Z_helper(arma::mat & Z_curr, // persumably large, thus will not copy
                     const arma::mat & data, 
                     const arma::mat & design,
                     const arma::vec & mu_curr,
                     const arma::mat & beta_curr,
                     const arma::mat & Omega_curr,
                     int k, int p, int n);

void update_Z_helper_CAR(arma::mat & Z_curr, 
                     const arma::mat & data, 
                     const arma::mat & design,
                     const arma::vec & mu_curr,
                     const arma::mat & beta_curr,
                     const arma::mat & Omega_curr,
                     int k, int p, int n);

// version for graphical LASSO
void update_Z_graphical_helper(arma::mat & Z_curr, 
                               const arma::mat & data, 
                               const arma::vec & mu_curr,
                               const arma::mat & Omega_curr,
                               int k, int n);
                     
#endif

