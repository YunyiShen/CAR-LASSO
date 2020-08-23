/*
 *This part was intend to calculate the FI matrix for CAR model, 
 *for experimental design aand active learning of such model
*/

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h> 
#include <tgmath.h>
using namespace Rcpp;
using namespace arma;
using namespace std;


// For coders: ALL THE INDECES START FROM 0!!

// This function calculate the second derivitave: 
//   \partial{l}^2/\partial{\beta_{ij}}\partial{beta}
//   result was a p by k matrix of above partial derivitives, useful in FI
//   need to be vectorized when constructing the FI

arma::mat sec_dev_betaij_beta( const arma::mat & design, // a row vector of design
                              const arma::mat & Sigma, // the cov mat
                              int i, int j){
    arma::mat XtX = design.t() * design;
    return(vectorise(-XtX.col(i) * Sigma.row(j)));
    
}

// This function calculate the second derivitave: 
//   \partial{l}^2/\partial{\beta_{ij}}\partial{mu}
//   result was a k by 1 matrix of above partial derivitives, useful in FI
arma::mat sec_dev_betaij_mu(const arma::mat & design, 
                            const arma::mat & Sigma,
                            int i, int j){
    return(-Sigma.col(j) * design(0,i));
}

// There is no need for partial l^2 / partial mu partial mu since it is -Sigma

// This function calculate the second derivitave: 
//   \partial{l}^2/\partial{\beta_{ij}}\partial{Omega(upper.tri)}
//   result was a k by 1 matrix of above partial derivitives, useful in FI
arma::mat sec_dev_betaij_Omega(const arma::mat & design, 
                               const arma::mat & Sigma,
                               const arma::mat & beta,
                               const arma::vec & mu,
                               int i, int j){

    arma::mat XtX = design.t() * design;
    arma::mat temp1,temp;
    temp1 = Sigma * mu * design; 
    temp = temp1.col(i) * Sigma.row(j) + Sigma * beta.t() * XtX.col(i) * Sigma.row(j);

    // get only the upper tri by chain rule
    temp += temp.t();
    temp.diag() /= 2;
    return(temp(trimatu_ind(size(temp))));
}


// This function calculate the second derivitave: 
//   \partial{l}^2/\partial{\mu_{j}}\partial{Omega(upper.tri)}
//   result was a .5*k*(k-1) by 1 matrix of above partial derivitives, useful in FI
arma::mat sec_dev_muj_Omega(const arma::mat & design, 
                            const arma::mat & Sigma,
                            const arma::mat & beta,
                            const arma::vec & mu,
                            int j){

    arma::mat temp;

    temp = Sigma.col(j) * (design * beta + mu.t()) * Sigma;

    // get only the upper tri by chain rule
    temp += temp.t();
    temp.diag() /= 2;
    return(temp(trimatu_ind(size(temp))));
}

// This function calculate the second derivitave: 
//   \partial{l}^2/\partial{\Omega_{sl}}\partial{Omega(upper.tri)}, s<=l (so upper tri)
//   result was a .5*k*(k-1) by 1 matrix of above partial derivitives, useful in FI
// [[Rcpp::export]]
arma::mat sec_dev_Omegasl_Omega(const arma::mat & design, 
                                const arma::mat & Sigma,
                                const arma::mat & beta,
                                const arma::vec & mu,
                                int s, int l){
    
    arma::mat muXbeta = mu * design * beta;
    arma::mat betatXtXbeta = (arma::trans(design * beta)) * design * beta;
    arma::mat mumut = mu * mu.t();
    arma::mat Sigma_s_times_l = Sigma.col(s) * Sigma.row(l);
    arma::mat Sigma_l_times_s = Sigma.col(l) * Sigma.row(s);
    arma::mat temp;
    // FIXME: simpliy the math, it was currently (probably) 
    //   the simpliest analytical form, but not for computation
    if(s == l){
        temp = -Sigma_s_times_l - .5 * Sigma * muXbeta.t() * Sigma_s_times_l -
               .5 * Sigma_s_times_l * muXbeta.t() * Sigma -
               .5 * Sigma *  muXbeta * Sigma_s_times_l - 
               .5 * Sigma_s_times_l * muXbeta * Sigma - 
               .5 * Sigma * betatXtXbeta * Sigma_s_times_l - 
               .5 * Sigma_s_times_l * betatXtXbeta * Sigma - 
               .5 * Sigma * mumut * Sigma_s_times_l - 
               .5 * Sigma_s_times_l * mumut * Sigma;

    }
    else{
        temp = - Sigma_s_times_l - Sigma_l_times_s - 
               .5 * Sigma * muXbeta.t() * Sigma_l_times_s - 
               .5 * Sigma_l_times_s * muXbeta.t() * Sigma - 
               .5 * Sigma * muXbeta.t() * Sigma_s_times_l - 
               .5 * Sigma_s_times_l * muXbeta.t() * Sigma - 
               .5 * Sigma * muXbeta * Sigma_l_times_s - 
               .5 * Sigma_l_times_s * muXbeta * Sigma - 
               .5 * Sigma * muXbeta * Sigma_s_times_l - 
               .5 * Sigma_s_times_l * muXbeta * Sigma - 
               .5 * Sigma * betatXtXbeta * Sigma_l_times_s - 
               .5 * Sigma_l_times_s * betatXtXbeta * Sigma - 
               .5 * Sigma * betatXtXbeta * Sigma_s_times_l - 
               .5 * Sigma_s_times_l * betatXtXbeta * Sigma-
               .5 * Sigma * mumut * Sigma_l_times_s - 
               .5 * Sigma_l_times_s * mumut * Sigma - 
               .5 * Sigma * mumut * Sigma_s_times_l - 
               .5 * Sigma_s_times_l * mumut * Sigma;

    }
    //Rcout << temp << endl;
    temp += temp.t();
    temp.diag() /= 2;
    arma::mat partial_uppertri = temp(trimatu_ind(size(temp)));
    
    return(partial_uppertri); 
}


// OK now is the time to construct the FI matrix, REMEMBER TO TAKE NEGATIVE!!!
// [[Rcpp::export]]
arma::mat CAR_FI(const arma::mat & design, 
                 const arma::mat & Omega,
                 const arma::mat & beta,
                 const arma::vec & mu,
                 int k, int p){
    int dimension = k * p + k + .5 * k * (k+1);
    arma::mat FI_mat(dimension,dimension,fill::zeros);
    arma::mat Sigma = inv(Omega);

    // positions for each parameter set:
    arma::uvec beta_ind = arma::linspace<arma::uvec>(0, k * p - 1 , k * p);
    arma::uvec mu_ind = k * p + arma::linspace<arma::uvec>(0, k - 1 , k );
    arma::uvec Omega_uptri_ind = (p + 1) * k + 
                                arma::linspace<arma::uvec>(0, .5 * k * (k+1) - 1 , 
                                                           .5 * k * (k+1) );
    arma::uvec ind_temp(1,fill::zeros);
    // These are not the most efficient loops, but good for now during developing 
    // First get beta parts, by column (so same with vectorise()): 
    // We will fill from 0 to k * p - 1 column of the FI
    
    for(int j = 0; j < k ; ++j){
        for(int i = 0 ; i < p ; ++i){
            ind_temp(0) = j * p + i; // we will work on this column 
            FI_mat(beta_ind,ind_temp) = sec_dev_betaij_beta(design, Sigma, i, j);
            FI_mat(mu_ind,ind_temp) = sec_dev_betaij_mu(design, Sigma, i, j);
            FI_mat(Omega_uptri_ind,ind_temp) = sec_dev_betaij_Omega(design, Sigma, beta, mu, i, j);

        }
    }
    
    // now we fill the mu part, 
    //  from k*p to (p+1) * k - 1 column
    FI_mat(mu_ind,mu_ind) = -Sigma;
    // other part for mu, basically just Omega
    for(int j = 0 ; j < k ; ++j){
        ind_temp(0) = j + k * p ; // work on this column
        FI_mat(Omega_uptri_ind,ind_temp) = 
                sec_dev_muj_Omega(design, Sigma,beta,mu,j);
    }
    //Rcout << "flag" <<endl;
    // flaged pass
    // now Omega part
    //  from (p+1) * k to [(p+1) * k + .5 * k * (k+1) - 1]
    ind_temp(0) = k * (p + 1)-1;
    for(int l = 0 ; l < k ; ++l){
        for(int s = 0 ; s <= l ; ++s){
            ind_temp(0) += 1;
            FI_mat(Omega_uptri_ind,ind_temp) = sec_dev_Omegasl_Omega(design,
                                            Sigma, beta, mu, s, l);
        }
    }
    
    // filling the blanks by blockwise symmetry:
    FI_mat(beta_ind,mu_ind) = arma::trans(FI_mat(mu_ind,beta_ind));
    FI_mat(beta_ind,Omega_uptri_ind) = arma::trans(FI_mat(Omega_uptri_ind,beta_ind));
    FI_mat(mu_ind,Omega_uptri_ind) = arma::trans(FI_mat(Omega_uptri_ind,mu_ind));
    return(-FI_mat);// FI is negative Hessian (for exponential family)
}


// [[Rcpp::export]]
arma::mat CAR_FI_graph(const arma::mat & design, 
                 const arma::mat & Omega,
                 const arma::mat & beta,
                 const arma::vec & mu,
                 int k, int p){
    int dimension = .5 * k * (k+1);
    arma::mat FI_mat(dimension,dimension,fill::zeros);
    arma::mat Sigma = inv(Omega);

    // These are not the most efficient loops, but good for now during developing 
    // now Omega part
    //  from (p+1) * k to [(p+1) * k + .5 * k * (k+1) - 1]
    int ind_temp = -1;
    for(int l = 0 ; l < k ; ++l){
        for(int s = 0 ; s <= l ; ++s){
            ind_temp += 1;
            FI_mat.col(ind_temp) = sec_dev_Omegasl_Omega(design,
                                            Sigma, beta, mu, s, l);
        }
    }
    
    return(-FI_mat);// FI is negative Hessian (for exponential family)
}

