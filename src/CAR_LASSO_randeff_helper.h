#ifndef CAR_LASSO_RANDEFF_HELPER_H
#define CAR_LASSO_RANDEFF_HELPER_H
arma::mat update_car_nu_helper(const arma::mat & data,
                               const arma::mat & design, // design matrix
                               const arma::mat & design_r, // design mat for random effect
                               const arma::mat & membership, // membership matrix for random precision
                               const arma::mat & beta,
                               const arma::vec & mu, // grand mean
                               const arma::mat & xi, // the precision vector of random effect, should be size of k
                               const arma::mat & Omega,
                               int k, int pr, int n);

void update_xi_helper(arma::mat & xi,
                      const arma::mat & nu,
                      const arma::mat & membership,
                      const double & alpha,
                      const double & beta, 
                      int k, int pr, int m);


void get_data_centered(arma::mat & centered_data,
                       const arma::mat & data,
                       const arma::mat & design_r,
                       const arma::mat & nu,
                       const arma::mat & Omega);

void update_car_randeff_Omega_helper(arma::mat & Omega,
                             const arma::mat & data,
                             const arma::mat & design,
                             const arma::mat & design_r,
                             const arma::mat & nu,
                             const arma::vec & mu,
                             const arma::mat & beta,
                             const double & lambda_curr,
                             int k, int p,int n);
#endif