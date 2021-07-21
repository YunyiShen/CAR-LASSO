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
#include "ars_multinomial_helper.h"
#include "ars_pois_helper.h"
#include "Probit_helper.h"

// [[Rcpp::export]]
List CAR_multireg_cpp(const arma::mat &data, const arma::mat &design,
                      int n_sample, const arma::mat &Bbar,
                      const arma::mat &A,
                      double nu, const arma::mat &V)
{

    int n = design.n_rows;
    int k = data.n_cols;
    int p = design.n_cols;

    arma::mat beta_mcmc(n_sample, k * p); // beta mcmc
    beta_mcmc += NA_REAL;

    arma::mat Omega_mcmc(n_sample, floor(static_cast<double>(k * (k + 1) / 2))); // vectorized column first, but had no diagnol
    Omega_mcmc += NA_REAL;

    arma::mat mu_mcmc(n_sample, k); // mean for node 1 to k
    mu_mcmc += NA_REAL;
    arma::mat design_full(n, p + 1, fill::ones);
    design_full.cols(1, p) = design;
    Progress prog(n_sample, false); // progress bar

    for (int i = 0; i < n_sample; ++i)
    {

        if (Progress::check_abort())
        {
            Rcerr << "keyboard abort\n";
            return (
                Rcpp::List::create(
                    Rcpp::Named("beta") = beta_mcmc,
                    Rcpp::Named("mu") = mu_mcmc,
                    Rcpp::Named("Omega") = Omega_mcmc));
        }

        rmultireg2 thissample(data, design_full, Bbar, A, nu, V); // sample from it.
        arma::mat B_hat = thissample.B * thissample.Omega;
        beta_mcmc.row(i) = trans(arma::vectorise(B_hat.rows(1, p)));
        mu_mcmc.row(i) = (B_hat.row(0));
        Omega_mcmc.row(i) = trans(thissample.Omega(trimatu_ind(size(thissample.Omega))));
        prog.increment();
    }

    return (
        Rcpp::List::create(
            Rcpp::Named("beta") = beta_mcmc,
            Rcpp::Named("mu") = mu_mcmc,
            Rcpp::Named("Omega") = Omega_mcmc));
}

// [[Rcpp::export]]
List Multinomial_CAR_multireg_cpp(const arma::mat &data, const arma::mat &design,
                              int n_burn_in, int n_iter, int thin_by,
                              const arma::mat &Bbar,
                              const arma::mat &A,
                              double nu, const arma::mat &V,
                              int ns, int m, double emax)
{

    int n = design.n_rows;
    int k = data.n_cols - 1;
    int p = design.n_cols;

    int n_save = floor(static_cast<double>(n_iter / thin_by)); //
    int i_save = 0;

    // mcmc matrices:
    arma::mat beta_mcmc(n_save, k * p); // beta mcmc
    beta_mcmc += NA_REAL;

    arma::mat Omega_mcmc(n_save, floor(static_cast<double>(k * (k + 1) / 2))); // vectorized column first, but had no diagnol
    Omega_mcmc += NA_REAL;

    arma::mat mu_mcmc(n_save, k); // mean for node 1 to k
    mu_mcmc += NA_REAL;

    arma::mat Z_curr(n, k, arma::fill::randn);
    arma::mat Omega_curr = eye(k, k);
    arma::mat B_curr(p + 1, k, arma::fill::randn);

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
                    Rcpp::Named("Omega") = Omega_mcmc));
        }

        rmultireg2 thissample(Z_curr, design_full, Bbar, A, nu, V); // sample from it.
        B_curr = thissample.B * thissample.Omega;
        Omega_curr = thissample.Omega;
        // update Z

        update_Z_helper_multinomial_CAR(Z_curr, data,
                                        design, arma::vectorise(B_curr.row(0)),
                                        B_curr.rows(1, p), Omega_curr, k, p, n, ns, m, emax);
        //Rcout << Z_curr(0) << endl;
        if ((i - n_burn_in) >= 0 && (i + 1 - n_burn_in) % thin_by == 0)
        {

            beta_mcmc.row(i_save) = trans(vectorise(B_curr.rows(1, p)));
            Omega_mcmc.row(i_save) = trans(Omega_curr(trimatu_ind(size(Omega_curr))));
            mu_mcmc.row(i_save) = B_curr.row(0);


            i_save++;
        }
        prog.increment();
    }

    return (
        Rcpp::List::create(
            Rcpp::Named("beta") = beta_mcmc,
            Rcpp::Named("mu") = mu_mcmc,
            Rcpp::Named("Omega") = Omega_mcmc));
}

// [[Rcpp::export]]
List Pois_CAR_multireg_cpp(const arma::mat &data, const arma::mat &design,
                              int n_burn_in, int n_iter, int thin_by,
                              const arma::mat &Bbar,
                              const arma::mat &A,
                              double nu, const arma::mat &V,
                              int ns, int m, double emax)
{

    int n = design.n_rows;
    int k = data.n_cols;
    int p = design.n_cols;

    int n_save = floor(static_cast<double>(n_iter / thin_by)); //
    int i_save = 0;

    // mcmc matrices:
    arma::mat beta_mcmc(n_save, k * p); // beta mcmc
    beta_mcmc += NA_REAL;

    arma::mat Omega_mcmc(n_save, floor(static_cast<double>(k * (k + 1) / 2))); // vectorized column first, but had no diagnol
    Omega_mcmc += NA_REAL;

    arma::mat mu_mcmc(n_save, k); // mean for node 1 to k
    mu_mcmc += NA_REAL;

    arma::mat Z_curr(n, k, arma::fill::randn);
    arma::mat Omega_curr = eye(k, k);
    arma::mat B_curr(p + 1, k, arma::fill::randn);

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
                    Rcpp::Named("Omega") = Omega_mcmc));
        }

        rmultireg2 thissample(Z_curr, design_full, Bbar, A, nu, V); // sample from it.
        B_curr = thissample.B * thissample.Omega;
        Omega_curr = thissample.Omega;
        // update Z

        update_Z_helper_Pois_CAR(Z_curr, data,
                                        design, arma::vectorise(B_curr.row(0)),
                                        B_curr.rows(1, p), Omega_curr, k, p, n, ns, m, emax);
        //Rcout << Z_curr(0) << endl;
        if ((i - n_burn_in) >= 0 && (i + 1 - n_burn_in) % thin_by == 0)
        {

            beta_mcmc.row(i_save) = trans(vectorise(B_curr.rows(1, p)));
            Omega_mcmc.row(i_save) = trans(Omega_curr(trimatu_ind(size(Omega_curr))));
            mu_mcmc.row(i_save) = B_curr.row(0);


            i_save++;
        }
        prog.increment();
    }

    return (
        Rcpp::List::create(
            Rcpp::Named("beta") = beta_mcmc,
            Rcpp::Named("mu") = mu_mcmc,
            Rcpp::Named("Omega") = Omega_mcmc));
}


// [[Rcpp::export]]
List Probit_CAR_multireg_cpp(const arma::mat &data, const arma::mat &design,
                              int n_burn_in, int n_iter, int thin_by,
                              const arma::mat &Bbar,
                              const arma::mat &A,
                              double nu, const arma::mat &V)
{

    int n = design.n_rows;
    int k = data.n_cols;
    int p = design.n_cols;

    int n_save = floor(static_cast<double>(n_iter / thin_by)); //
    int i_save = 0;

    // mcmc matrices:
    arma::mat beta_mcmc(n_save, k * p); // beta mcmc
    beta_mcmc += NA_REAL;

    arma::mat Omega_mcmc(n_save, floor(static_cast<double>(k * (k + 1) / 2))); // vectorized column first, but had no diagnol
    Omega_mcmc += NA_REAL;

    arma::mat mu_mcmc(n_save, k); // mean for node 1 to k
    mu_mcmc += NA_REAL;

    arma::mat Z_curr(n, k, arma::fill::randn);
    arma::mat Omega_curr = eye(k, k);
    arma::mat B_curr(p + 1, k, arma::fill::randn);

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
                    Rcpp::Named("Omega") = Omega_mcmc));
        }

        rmultireg2 thissample(Z_curr, design_full, Bbar, A, nu, V); // sample from it.
        B_curr = thissample.B * thissample.Omega;
        Omega_curr = thissample.Omega;
        // update Z

        update_Z_helper_CAR(Z_curr, data,
                                        design, arma::vectorise(B_curr.row(0)),
                                        B_curr.rows(1, p), Omega_curr, k, p, n);
        //Rcout << Z_curr(0) << endl;
        if ((i - n_burn_in) >= 0 && (i + 1 - n_burn_in) % thin_by == 0)
        {

            beta_mcmc.row(i_save) = trans(vectorise(B_curr.rows(1, p)));
            Omega_mcmc.row(i_save) = trans(Omega_curr(trimatu_ind(size(Omega_curr))));
            mu_mcmc.row(i_save) = B_curr.row(0);


            i_save++;
        }
        prog.increment();
    }

    return (
        Rcpp::List::create(
            Rcpp::Named("beta") = beta_mcmc,
            Rcpp::Named("mu") = mu_mcmc,
            Rcpp::Named("Omega") = Omega_mcmc));
}