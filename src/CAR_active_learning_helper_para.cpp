// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h> 
#include <tgmath.h>
#include "CAR_FI_helper.h"
using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace RcppParallel;


// parallel version:

struct g_det_worker : public Worker{
    // inputs
    const arma::mat & Omega_mcmc;
    const arma::mat & beta_mcmc;
    const arma::mat & mu_mcmc;
    const arma::mat & new_design;
    const arma::mat & old_design;
    const int & k; 
    const int & p; 
    const int & n_new;
    const int & n_old;
    //output vec
    arma::vec & log_dets;

    // initialize
    g_det_worker(const arma::mat & Omega_mcmc,
                const arma::mat & beta_mcmc,
                const arma::mat & mu_mcmc,
                const arma::mat & new_design,
                const arma::mat & old_design,
                const int & k, const int & p, 
                const int & n_new, const int & n_old,
                arma::vec & log_dets):
        Omega_mcmc(Omega_mcmc), beta_mcmc(beta_mcmc),
        mu_mcmc(mu_mcmc), new_design(new_design),
        old_design(old_design), k(k), p(p), n_new(n_new),
        n_old(n_old), log_dets(log_dets) {}
    
    void operator()(std::size_t begin, std::size_t end){
        arma::mat Omega_temp(k,k,fill::zeros);
        arma::mat beta_temp(p,k,fill::zeros);
        arma::vec mu_temp(k,fill::zeros);
        double val;
        double sign;
        int dimension = .5 * k * (k+1);   
        arma::uvec uppertri_graph = trimatu_ind(size(Omega_temp));
        arma::mat FI_mat(dimension,dimension,fill::zeros);
        // loop over all posterior samples
        for(int i = begin ; i < end ; ++i){
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
        
            log_dets(i) = val;
        }

    }
};

// [[Rcpp::export]]
arma::vec G_det_para(const arma::mat & new_design,
                          const arma::mat & old_design,
                          const Rcpp::List & CAR_model,
                          int k, int p, int n_new, int n_old){
    const arma::mat & Omega_mcmc = CAR_model["Omega"];
    const arma::mat & beta_mcmc = CAR_model["beta"];
    const arma::mat & mu_mcmc = CAR_model["mu"];
    int n_mcmc = Omega_mcmc.n_rows;
    arma::vec log_dets(n_mcmc);
    g_det_worker g_det(Omega_mcmc, beta_mcmc,mu_mcmc,
                new_design,old_design,k, p, 
                n_new, n_old,
                log_dets);
    
    parallelFor(0,n_mcmc,g_det);
    return(log_dets);
}

// [[Rcpp::export]]
double expected_G_det_para(const arma::vec & new_design,
                          const arma::mat & old_design,
                          const Rcpp::List & CAR_model,
                          int k, int p, int n_new, int n_old){
    arma::mat new_design_mat = arma::reshape(new_design,n_new,p);

    return(mean(
        G_det_para(new_design_mat,old_design,CAR_model,
                           k,  p, n_new, n_old)
    ));
}




// trace version:

struct g_tr_worker : public Worker{
    // inputs
    const arma::mat & Omega_mcmc;
    const arma::mat & beta_mcmc;
    const arma::mat & mu_mcmc;
    const arma::mat & new_design;
    const arma::mat & old_design;
    const int & k; 
    const int & p; 
    const int & n_new;
    const int & n_old;
    //output vec
    arma::vec & log_dets;
    
    // initialize
    g_tr_worker(const arma::mat & Omega_mcmc,
                 const arma::mat & beta_mcmc,
                 const arma::mat & mu_mcmc,
                 const arma::mat & new_design,
                 const arma::mat & old_design,
                 const int & k, const int & p, 
                 const int & n_new, const int & n_old,
                 arma::vec & log_dets):
        Omega_mcmc(Omega_mcmc), beta_mcmc(beta_mcmc),
        mu_mcmc(mu_mcmc), new_design(new_design),
        old_design(old_design), k(k), p(p), n_new(n_new),
        n_old(n_old), log_dets(log_dets) {}
    
    void operator()(std::size_t begin, std::size_t end){
        arma::mat Omega_temp(k,k,fill::zeros);
        arma::mat beta_temp(p,k,fill::zeros);
        arma::vec mu_temp(k,fill::zeros);
        double val;
        double sign;
        int dimension = .5 * k * (k+1);   
        arma::uvec uppertri_graph = trimatu_ind(size(Omega_temp));
        arma::mat FI_mat(dimension,dimension,fill::zeros);
        // loop over all posterior samples
        for(int i = begin ; i < end ; ++i){
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
            //arma::log_det(val, sign, FI_mat);
            
            log_dets(i) = log(sum(FI_mat.diag()));
        }
        
    }
};

// [[Rcpp::export]]
arma::vec G_tr_para(const arma::mat & new_design,
                     const arma::mat & old_design,
                     const Rcpp::List & CAR_model,
                     int k, int p, int n_new, int n_old){
    const arma::mat & Omega_mcmc = CAR_model["Omega"];
    const arma::mat & beta_mcmc = CAR_model["beta"];
    const arma::mat & mu_mcmc = CAR_model["mu"];
    int n_mcmc = Omega_mcmc.n_rows;
    arma::vec log_dets(n_mcmc);
    g_tr_worker g_det(Omega_mcmc, beta_mcmc,mu_mcmc,
                       new_design,old_design,k, p, 
                       n_new, n_old,
                       log_dets);
    
    parallelFor(0,n_mcmc,g_det);
    return(log_dets);
}

// [[Rcpp::export]]
double expected_G_tr_para(const arma::vec & new_design,
                           const arma::mat & old_design,
                           const Rcpp::List & CAR_model,
                           int k, int p, int n_new, int n_old){
    arma::mat new_design_mat = arma::reshape(new_design,n_new,p);
    
    return(mean(
            G_tr_para(new_design_mat,old_design,CAR_model,
                       k,  p, n_new, n_old)
    ));
}


// trace for the full FI

struct full_tr_worker : public Worker{
    // inputs
    const arma::mat & Omega_mcmc;
    const arma::mat & beta_mcmc;
    const arma::mat & mu_mcmc;
    const arma::mat & new_design;
    const arma::mat & old_design;
    const int & k; 
    const int & p; 
    const int & n_new;
    const int & n_old;
    //output vec
    arma::vec & log_dets;
    
    // initialize
    full_tr_worker(const arma::mat & Omega_mcmc,
                 const arma::mat & beta_mcmc,
                 const arma::mat & mu_mcmc,
                 const arma::mat & new_design,
                 const arma::mat & old_design,
                 const int & k, const int & p, 
                 const int & n_new, const int & n_old,
                 arma::vec & log_dets):
        Omega_mcmc(Omega_mcmc), beta_mcmc(beta_mcmc),
        mu_mcmc(mu_mcmc), new_design(new_design),
        old_design(old_design), k(k), p(p), n_new(n_new),
        n_old(n_old), log_dets(log_dets) {}
    
    void operator()(std::size_t begin, std::size_t end){
        arma::mat Omega_temp(k,k,fill::zeros);
        arma::mat beta_temp(p,k,fill::zeros);
        arma::vec mu_temp(k,fill::zeros);
        double val;
        double sign;
        int dimension = k * p + k + .5 * k * (k+1);   
        arma::uvec uppertri_graph = trimatu_ind(size(Omega_temp));
        arma::mat FI_mat(dimension,dimension,fill::zeros);
        // loop over all posterior samples
        for(int i = begin ; i < end ; ++i){
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
                FI_mat += CAR_FI(old_design.row(j), 
                                       Omega_temp,beta_temp,mu_temp,
                                       k, p);
            }
            
            // new design
            for(int j = 0 ; j < n_new ; ++j){
                FI_mat += CAR_FI(new_design.row(j), 
                                       Omega_temp,beta_temp,mu_temp,
                                       k, p);
            }
            //arma::log_det(val, sign, FI_mat);
            
            log_dets(i) = log(sum(FI_mat.diag()));
        }
        
    }
};

// [[Rcpp::export]]
arma::vec full_tr_para(const arma::mat & new_design,
                     const arma::mat & old_design,
                     const Rcpp::List & CAR_model,
                     int k, int p, int n_new, int n_old){
    const arma::mat & Omega_mcmc = CAR_model["Omega"];
    const arma::mat & beta_mcmc = CAR_model["beta"];
    const arma::mat & mu_mcmc = CAR_model["mu"];
    int n_mcmc = Omega_mcmc.n_rows;
    arma::vec log_dets(n_mcmc);
    full_tr_worker full_tr(Omega_mcmc, beta_mcmc,mu_mcmc,
                       new_design,old_design,k, p, 
                       n_new, n_old,
                       log_dets);
    
    parallelFor(0,n_mcmc,full_tr);
    return(log_dets);
}

// [[Rcpp::export]]
double expected_all_tr_para(const arma::vec & new_design,
                           const arma::mat & old_design,
                           const Rcpp::List & CAR_model,
                           int k, int p, int n_new, int n_old){
    arma::mat new_design_mat = arma::reshape(new_design,n_new,p);
    
    return(mean(
            full_tr_para(new_design_mat,old_design,CAR_model,
                       k,  p, n_new, n_old)
    ));
}





// tr for only beta

// trace for the full FI

struct beta_tr_worker : public Worker{
    // inputs
    const arma::mat & Omega_mcmc;
    const arma::mat & beta_mcmc;
    const arma::mat & mu_mcmc;
    const arma::mat & new_design;
    const arma::mat & old_design;
    const int & k; 
    const int & p; 
    const int & n_new;
    const int & n_old;
    //output vec
    arma::vec & log_dets;
    
    // initialize
    beta_tr_worker(const arma::mat & Omega_mcmc,
                 const arma::mat & beta_mcmc,
                 const arma::mat & mu_mcmc,
                 const arma::mat & new_design,
                 const arma::mat & old_design,
                 const int & k, const int & p, 
                 const int & n_new, const int & n_old,
                 arma::vec & log_dets):
        Omega_mcmc(Omega_mcmc), beta_mcmc(beta_mcmc),
        mu_mcmc(mu_mcmc), new_design(new_design),
        old_design(old_design), k(k), p(p), n_new(n_new),
        n_old(n_old), log_dets(log_dets) {}
    
    void operator()(std::size_t begin, std::size_t end){
        arma::mat Omega_temp(k,k,fill::zeros);
        arma::mat beta_temp(p,k,fill::zeros);
        arma::vec mu_temp(k,fill::zeros);
        double val;
        double sign;
        int dimension = k * p + k ;   
        arma::uvec uppertri_graph = trimatu_ind(size(Omega_temp));
        arma::mat FI_mat(dimension,dimension,fill::zeros);
        // loop over all posterior samples
        for(int i = begin ; i < end ; ++i){
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
                FI_mat += CAR_FI_beta(old_design.row(j), 
                                       Omega_temp,beta_temp,mu_temp,
                                       k, p);
            }
            
            // new design
            for(int j = 0 ; j < n_new ; ++j){
                FI_mat += CAR_FI_beta(new_design.row(j), 
                                       Omega_temp,beta_temp,mu_temp,
                                       k, p);
            }
            //arma::log_det(val, sign, FI_mat);
            
            log_dets(i) = log(sum(FI_mat.diag()));
        }
        
    }
};

// [[Rcpp::export]]
arma::vec beta_tr_para(const arma::mat & new_design,
                     const arma::mat & old_design,
                     const Rcpp::List & CAR_model,
                     int k, int p, int n_new, int n_old){
    const arma::mat & Omega_mcmc = CAR_model["Omega"];
    const arma::mat & beta_mcmc = CAR_model["beta"];
    const arma::mat & mu_mcmc = CAR_model["mu"];
    int n_mcmc = Omega_mcmc.n_rows;
    arma::vec log_dets(n_mcmc);
    beta_tr_worker beta_tr(Omega_mcmc, beta_mcmc,mu_mcmc,
                       new_design,old_design,k, p, 
                       n_new, n_old,
                       log_dets);
    
    parallelFor(0,n_mcmc,beta_tr);
    return(log_dets);
}

// [[Rcpp::export]]
double expected_beta_tr_para(const arma::vec & new_design,
                           const arma::mat & old_design,
                           const Rcpp::List & CAR_model,
                           int k, int p, int n_new, int n_old){
    arma::mat new_design_mat = arma::reshape(new_design,n_new,p);
    
    return(mean(
            beta_tr_para(new_design_mat,old_design,CAR_model,
                       k,  p, n_new, n_old)
    ));
}


