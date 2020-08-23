#ifndef CAR_FI_HELPER_H
#define CAR_FI_HELPER_H

arma::mat sec_dev_betaij_beta( const arma::mat & design, // a row vector of design
                              const arma::mat & Sigma, // the cov mat
                              int i, int j);

arma::mat sec_dev_betaij_mu(const arma::mat & design, 
                            const arma::mat & Sigma,
                            int i, int j);

arma::mat sec_dev_betaij_Omega(const arma::mat & design, 
                               const arma::mat & Sigma,
                               const arma::mat & beta,
                               const arma::vec & mu,
                               int i, int j);

arma::mat sec_dev_muj_Omega(const arma::mat & design, 
                            const arma::mat & Sigma,
                            const arma::mat & beta,
                            const arma::vec & mu,
                            int j);

arma::mat sec_dev_Omegasl_Omega(const arma::mat & design, 
                                const arma::mat & Sigma,
                                const arma::mat & beta,
                                const arma::vec & mu,
                                int s, int l);

arma::mat CAR_FI(const arma::mat & design, 
                 const arma::mat & Omega,
                 const arma::mat & beta,
                 const arma::vec & mu,
                 int k, int p);

arma::mat CAR_FI_graph(const arma::mat & design, 
                 const arma::mat & Omega,
                 const arma::mat & beta,
                 const arma::vec & mu,
                 int k, int p);
                 
#endif