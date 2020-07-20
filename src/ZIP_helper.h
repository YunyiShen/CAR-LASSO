# ifndef ZIP_HELPER_H
# define ZIP_HELPER_H
void ZIP_update_Pois_latent_helper(mat & Pois_latent,
                                   const mat & data,
                                   const mat & Z_P,// Z for poiss
                                   const mat & Z_ZI,// Z for probit (ZI)
                                   int k, int p, int n);

void ZIP_update_Z_ZI_helper(mat & Z_ZI,
                            const mat & data,
                            const arma::mat & design,
                            const arma::vec & mu_curr,
                            const arma::mat & beta_curr,
                            const arma::mat & Omega_curr,
                            const mat & Pois_latent,
                            int k, int p, int n);



# endif
