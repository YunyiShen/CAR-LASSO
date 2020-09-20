/* This file was intend to bind FI calculation and posterior sampling 
* We will calculate the expected (or some other summary statistics) of 
* some information criterior
*/
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h> 
#include <tgmath.h>
#include <string>
#include "CAR_FI_helper.h"
using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace RcppParallel;



// [[Rcpp::export]]
arma::vec expected_G_eigen(const arma::mat & new_design,
                          const arma::mat & old_design,
                          const Rcpp::List & CAR_model,
                          int k, int p, int n_new, int n_old,
                          const string & summary ){
    const arma::mat & Omega_mcmc = CAR_model["Omega"];
    const arma::mat & beta_mcmc = CAR_model["beta"];
    const arma::mat & mu_mcmc = CAR_model["mu"];

    arma::mat Omega_temp(k,k,fill::zeros);
    arma::mat beta_temp(p,k,fill::zeros);
    arma::vec mu_temp(k,fill::zeros);
    int dimension = .5 * k * (k+1);
    arma::vec eigen_temp(dimension);
    
    arma::uvec uppertri_graph = trimatu_ind(size(Omega_temp));
    int n_mcmc = Omega_mcmc.n_rows;
    arma::vec summary_eigen(n_mcmc);
    
    arma::mat FI_mat(dimension,dimension,fill::zeros);
    for(int i = 0 ; i < n_mcmc ; ++i){
        FI_mat.zeros();
        Omega_temp.zeros();
        Omega_temp(uppertri_graph) = arma::trans(Omega_mcmc.row(i));
        Omega_temp += Omega_temp.t();
        Omega_temp.diag() /= 2;

        beta_temp = reshape(beta_mcmc.row(i),size(beta_temp));
        mu_temp = arma::trans(mu_mcmc.row(i));
        // generate the sub-FI for Graph of the old design
        // FIXME: not efficient at all
        for(int j = 0 ; j < n_old ; ++j){
            FI_mat += CAR_FI_graph(old_design.row(j), 
                 Omega_temp,beta_temp,mu_temp,
                 k, p);
        }
        // new design
        for(int j = 0 ; j < n_new ; ++j){
            FI_mat += CAR_FI_graph(new_design.row(j), 
                 Omega_temp,beta_temp,mu_temp,
                 k, p);
        }
        eig_sym( eigen_temp, FI_mat );
        if(summary == "max"){
            summary_eigen(i) = log( arma::max(eigen_temp));
        }
        else{
            summary_eigen(i) = summary == "median" ? 
                log( median(eigen_temp)) : 
                log( min(eigen_temp));
        }
    }
    return(summary_eigen);
}



// [[Rcpp::export]]
arma::vec G_det(const arma::mat & new_design,
                          const arma::mat & old_design,
                          const Rcpp::List & CAR_model,
                          int k, int p, int n_new, int n_old){
    const arma::mat & Omega_mcmc = CAR_model["Omega"];
    const arma::mat & beta_mcmc = CAR_model["beta"];
    const arma::mat & mu_mcmc = CAR_model["mu"];

    arma::mat Omega_temp(k,k,fill::zeros);
    arma::mat beta_temp(p,k,fill::zeros);
    arma::vec mu_temp(k,fill::zeros);
    double val;
    double sign;
    int dimension = .5 * k * (k+1);
    
    arma::uvec uppertri_graph = trimatu_ind(size(Omega_temp));
    int n_mcmc = Omega_mcmc.n_rows;
    arma::vec logdets(n_mcmc);
    
    arma::mat FI_mat(dimension,dimension,fill::zeros);
    for(int i = 0 ; i < n_mcmc ; ++i){
        FI_mat.zeros();
        Omega_temp.zeros();
        Omega_temp(uppertri_graph) = arma::trans(Omega_mcmc.row(i));
        Omega_temp += Omega_temp.t();
        Omega_temp.diag() /= 2;
        
        
        beta_temp = reshape(beta_mcmc.row(i),size(beta_temp));

        mu_temp = arma::trans(mu_mcmc.row(i));

        // generate the sub-FI for Graph of the old design
        // FIXME: not efficient at all
        for(int j = 0 ; j < n_old ; ++j){
            FI_mat += CAR_FI_graph(old_design.row(j), 
                 Omega_temp,beta_temp,mu_temp,
                 k, p);
        }
        
        // new design
        for(int j = 0 ; j < n_new ; ++j){
            FI_mat += CAR_FI_graph(new_design.row(j), 
                 Omega_temp,beta_temp,mu_temp,
                 k, p);
        }
        arma::log_det(val, sign, FI_mat);
        
        logdets(i) = val;
    }
    return(logdets);
}


// [[Rcpp::export]]
double expected_G_det(const arma::vec & new_design,
                          const arma::mat & old_design,
                          const Rcpp::List & CAR_model,
                          int k, int p, int n_new, int n_old){
    arma::mat new_design_mat = arma::reshape(new_design,n_new,p);

    return(mean(
        G_det(new_design_mat,old_design,CAR_model,
                           k,  p, n_new, n_old)
    ));
}



