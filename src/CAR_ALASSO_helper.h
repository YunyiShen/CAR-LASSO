#ifndef CAR_ALASSO_HELPER_H
#define CAR_ALASSO_HELPER_H

void update_car_Omega_adp_helper(arma::mat & Omega,
                             const arma::mat & data,
                             const arma::mat & design,
                             const arma::vec & mu,
                             const arma::mat & beta,
                             const arma::vec & lambda_curr,// different lambda for different entries
                             int k, int p,int n);


//the one allow panalty on digonal entries
void update_car_Omega_adp_helper2(arma::mat & Omega,
                             const arma::mat & data,
                             const arma::mat & design,
                             const arma::vec & mu,
                             const arma::mat & beta,
                             const arma::vec & lambda_curr,// different lambda for different entries
                             const arma::vec & lambda_diag, 
                             int k, int p,int n);
                             
arma::vec update_car_tau2_adp_helper(const arma::mat & beta,
                             const arma::vec & lambda2,
                             const arma::mat & Omega,
                             int k, int p, int n);

void update_car_lambda_Omega_adp_helper(arma::vec & lambda_curr,
                                        const arma::mat & Omega,
                                        const arma::vec & r,
                                        const arma::vec & delta
                                        );

void update_car_lambda2_beta_adp_helper(arma::vec & lambda2_curr,
                                        const arma::vec & tau2,
                                        const arma::mat & r,
                                        const arma::mat & delta,
                                        int k, int p
                                        );
#endif