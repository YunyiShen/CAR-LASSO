// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include <tgmath.h>
using namespace Rcpp;
using namespace arma;

#include <progress.hpp>
#include <progress_bar.hpp>
#include "rmultireg2_rcpp.h"
#include "CAR_LASSO_randeff_helper.h"

// [[Rcpp::export]]
List CAR_multireg_randeff_cpp(const arma::mat &data, const arma::mat &design,
                              const arma::mat &design_r, // design matrix for random effect, each ROW as a sample
                              const arma::mat &membership,
                              int n_burn_in, int n_iter, int thin_by,
                              const arma::mat &Bbar,
                              const arma::mat &A,
                              double nu, const arma::mat &V,
                              const double alpha, // prior for random effect
                              const double beta)  // prior for random effect)
{

    int k = data.n_cols;                  // number of nodes
    int p = design.n_cols;                //number of predictors
    int pr = design_r.n_cols;             // number of unique random effects catagories
    int n = data.n_rows;                  // number of samples
    int m = membership.n_cols;            // number of mumbers;
    int n_save = floor(n_iter / thin_by); //
    int i_save = 0;

    // mcmc matrices:
    arma::mat beta_mcmc(n_save, k * p); // beta mcmc
    beta_mcmc += NA_REAL;

    arma::mat Omega_mcmc(n_save, floor(k * (k + 1) / 2)); // vectorized column first, but had no diagnol
    Omega_mcmc += NA_REAL;

    arma::mat mu_mcmc(n_save, k); // mean for node 1 to k
    mu_mcmc += NA_REAL;

    arma::mat nu_mcmc(n_save, k * pr); // latent random effect
    nu_mcmc += NA_REAL;

    arma::mat xi_mcmc(n_save, k * m); // precision of random effect
    xi_mcmc += NA_REAL;

    arma::vec mu_curr = trans(mean(data)); // current value of mean
    arma::mat centered_data = data;
    centered_data.each_row() -= mu_curr.t(); // now for initialize beta, will be reused for substracting random effect in the main loop
    arma::mat Omega_curr(k, k);              // current value of Omega
    Omega_curr = inv(cov(data));
    arma::mat beta_curr = solve(design.t() * design, design.t() * (centered_data * Omega_curr)); // current value of beta

    arma::mat nu_curr(pr, k, arma::fill::zeros);
    arma::mat xi_curr(m, k, arma::fill::zeros);
    xi_curr += 1;

    arma::mat B_curr;

    arma::mat design_full(n, p + 1, arma::fill::ones);
    design_full.cols(1, p) = design;
    Progress prog(n_iter + n_burn_in, false); // progress bar

    for (int i = 0; i < (n_iter + n_burn_in); ++i)
    {

        if (Progress::check_abort())
        {
            Rcerr << "keyboard abort\n";
            return (
                Rcpp::List::create(
                    Rcpp::Named("beta") = beta_mcmc,
                    Rcpp::Named("mu") = mu_mcmc,
                    Rcpp::Named("nu") = nu_mcmc,
                    Rcpp::Named("xi") = xi_mcmc,
                    Rcpp::Named("Omega") = Omega_mcmc));
        }

        get_data_centered(centered_data, data, design_r, nu_curr, Omega_curr);

        rmultireg2 thissample(centered_data, design_full, Bbar, A, nu, V); // sample from it.
        B_curr = thissample.B * thissample.Omega;
        beta_curr = B_curr.rows(1, p);
        mu_curr = trans(B_curr.row(0));
        Omega_curr = thissample.Omega;
        // update Z

        nu_curr = update_car_nu_helper(data, design, design_r, membership, beta_curr,
                                       mu_curr, xi_curr, Omega_curr,
                                       k, pr, n);

        update_xi_helper(xi_curr, nu_curr, membership, alpha, beta, k, pr, m);

        if ((i - n_burn_in) >= 0 && (i + 1 - n_burn_in) % thin_by == 0)
        {

            beta_mcmc.row(i_save) = trans(vectorise(beta_curr));
            Omega_mcmc.row(i_save) = trans(Omega_curr(trimatu_ind(size(Omega_curr))));
            mu_mcmc.row(i_save) = mu_curr.t();
            nu_mcmc.row(i_save) = trans(vectorise(nu_curr));
            xi_mcmc.row(i_save) = trans(vectorise(xi_curr));

            i_save++;
        }
        prog.increment();
    }

    return (
        Rcpp::List::create(
            Rcpp::Named("beta") = beta_mcmc,
            Rcpp::Named("mu") = mu_mcmc,
            Rcpp::Named("nu") = nu_mcmc,
            Rcpp::Named("xi") = xi_mcmc,
            Rcpp::Named("Omega") = Omega_mcmc));
}
