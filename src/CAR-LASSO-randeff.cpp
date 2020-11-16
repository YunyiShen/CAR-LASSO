// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include <tgmath.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
#include "CAR_LASSO_helper.h"
#include "CAR_LASSO_randeff_helper.h"

/*
 * Main sampling functions for Conditional Auto Regression LASSO, 
 * Basic idea was to embed Graphical LASSO into a normal LASSO using the hirechical structure
 *   described by Wang (2012) and Park and Casella 2008
 * In this model average structure offer some extra information of conditional correlation
 * 
 * This one was not adaptive, i.e. lambda_Omega was a fixed for all entries of Omega and 
 *   similar for betas
 * 
 * A CAR w/ random effect can be reparameterize into a model s.t.: 
 * Y_i~N(Sigma (X_ibeta+mu+nu_{i,.}),Sigma)
 *     Here we have nu matrix similar to beta matrix, each column (random effect for a certain node) 
 *        was iid while rows were independent, that's for different nodes
 * nu_{i,j}~N(0,tt_{i,j}) tt_{ij} is the random effect's precision, we also have
 * tt = membership * xi, membership is the design matrix for variance component, must be orthagonal and 0,1
 * xi_{i,j}~Gamma(alpha,beta) iid // we have some prior on the precision
 */

/* 
Input:
  @ data: a matrix with column as nodes and row as samples
  @ design: a design matrix of common input to the network, should have same # of rows as data
  @ design_r: a design matrix for the random effect
  @ membership: "design" matrix for variance component of random effect
  @ n_iter: number of saved sampling iterations in the Gibbs sampler
  @ n_burn_in: number of burn in before sampling
  @ thin_by: subsampling steps, integer
  @ r_beta, r_Omega: shape parameter for shrinkage parameter lambda of beta and Omega
  @ delta_beta, delta_Omega: RATE parameter for lambda prior
  @ progress: whether to show a progress bar from C++

Output:
  A list with component:
  @ beta: a matrix with each row as an MCMC sample, 
    columns are the vectorization of beta, 
    while beta matrix has p row and k columns
  @ mu: a matrix with each row as an MCMC sample, columns are intercept vectors
  @ nu: a matrix with each row as an MCMC sample, columns as latent random effect
  @ Omega: a matrix with each row as an MCMC sample, 
    columns are the upper diagnol entries of precision matrix Omega
  @ xi: a matrix with each row as an MCMC sample, 
    columns are the precision of random effect
  @ lambda: a matrix with only two columns, first was for beta, second was for Omega
    each row was an MCMC sample of shrinkage parameter lambda


*/

// [[Rcpp::export]]
List CAR_LASSO_randeff_Cpp(const arma::mat & data,     // col composition data, ROW as a sample
                           const arma::mat & design,   // design matrix, each ROW as a sample
                           const arma::mat & design_r, // design matrix for random effect, each ROW as a sample
                           const arma::mat & membership, //design matrix for variance component, should have pr rows and m columns
                           const int n_iter,          // how many iterations?
                           const int n_burn_in,       // burn in
                           const int thin_by,         // thinning?
                           const double r_beta,       // prior on lambda of beta
                           const double delta_beta,
                           const double r_Omega, // prior on lambda of Omega
                           const double delta_Omega,
                           const double alpha, // prior for random effect
                           const double beta,  // prior for random effect
                           bool progress)
{                             // whether to report progress
    int k = data.n_cols;      // number of nodes
    int p = design.n_cols;    //number of predictors
    int pr = design_r.n_cols; // number of unique random effects catagories
    int n = data.n_rows;      // number of samples
    int m = membership.n_cols; // number of mumbers;
    int n_save = floor(n_iter / thin_by); //
    int i_save = 0;

    // mcmc matrices:
    arma::mat beta_mcmc(n_save, k * p); // beta mcmc
    beta_mcmc += NA_REAL;

    arma::mat Omega_mcmc(n_save, floor(k * (k + 1) / 2)); // vectorized column first, but had no diagnol
    Omega_mcmc += NA_REAL;

    arma::mat mu_mcmc(n_save, k); // mean for node 1 to k
    mu_mcmc += NA_REAL;

    arma::mat nu_mcmc(n_save, k * pr); // LASSO parameter for beta and B
    nu_mcmc += NA_REAL;

    arma::mat xi_mcmc(n_save, k*m); // LASSO parameter for beta and B
    xi_mcmc += NA_REAL;

    arma::mat lambda_mcmc(n_save, 2); // LASSO parameter for beta and B
    lambda_mcmc += NA_REAL;

    //arma::vec tau2_curr = randg<arma::vec> (k*p,distr_param(r_beta,delta_beta)); // current latent variable tau^2, for prior of beta
    arma::vec tau2_curr = randg<arma::vec>(k * p, distr_param(1.0, .01)); // current latent variable tau^2, for prior of beta

    //Rcout << tau2_curr <<endl;
    arma::vec mu_curr = trans(mean(data)); // current value of mean
    arma::mat centered_data = data;
    centered_data.each_row() -= mu_curr.t(); // now for initialize beta, will be reused for substracting random effect in the main loop
    arma::mat Omega_curr(k, k);              // current value of Omega
    Omega_curr = inv(cov(data));
    arma::mat beta_curr = solve(design.t() * design, design.t() * (centered_data * Omega_curr)); // current value of beta

    arma::mat nu_curr(pr, k, arma::fill::zeros);
    arma::mat xi_curr(m,k, arma::fill::zeros);
    xi_curr += 1;

    double lambda2_beta = randg<double>(distr_param(r_beta, 1 / delta_beta));   // current value of squared LASSO parameter of \beta
    double lambda_Omega = randg<double>(distr_param(r_Omega, 1 / delta_Omega)); // current value of squared LASSO parameter of B

    double Omega_r_post = (r_Omega + (k * (k + 1) / 2));
    double Omega_delta_post;

    Progress prog((n_iter + n_burn_in), progress); // progress bar

    // main loop
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
                    Rcpp::Named("Omega") = Omega_mcmc,
                    Rcpp::Named("lambda") = lambda_mcmc));
        }
        // block update start:
        // calculate centered data
        get_data_centered(centered_data, data, design_r, nu_curr, Omega_curr);

        // Update lambda_Omega
        Omega_delta_post = (delta_Omega + sum(sum(abs(Omega_curr))) / 2);
        lambda_Omega = R::rgamma(Omega_r_post, 1 / Omega_delta_post);

        //Update betas:
        beta_curr = update_car_beta_helper(centered_data, design, mu_curr,
                                           tau2_curr, Omega_curr,
                                           k, p, n);

        // update Omega

        update_car_Omega_helper(Omega_curr, centered_data, design,
                                mu_curr, beta_curr,
                                lambda_Omega,
                                k, p, n);

        get_data_centered(centered_data, data, design_r, nu_curr, Omega_curr);

        // Update mu

        mu_curr = update_car_mu_helper(centered_data, design, beta_curr,
                                       Omega_curr,
                                       k, p, n);

        // Update tau2 for beta
        tau2_curr = update_car_tau2_helper(beta_curr, lambda2_beta,
                                           Omega_curr, k, p, n);

        // Update lambda_beta

        lambda2_beta = R::rgamma(r_beta + k * p, 1 / (delta_beta + sum(tau2_curr) / 2));

        // update random effect
        nu_curr = update_car_nu_helper(data, design, design_r, membership,beta_curr,
                                       mu_curr, xi_curr, Omega_curr,
                                       k, pr, n);

        update_xi_helper(xi_curr, nu_curr, membership,alpha, beta, k, pr, m);

        // saving the state
        if ((i - n_burn_in) >= 0 && (i + 1 - n_burn_in) % thin_by == 0)
        {

            beta_mcmc.row(i_save) = trans(vectorise(beta_curr));
            Omega_mcmc.row(i_save) = trans(Omega_curr(trimatu_ind(size(Omega_curr))));
            mu_mcmc.row(i_save) = mu_curr.t();
            nu_mcmc.row(i_save) = trans(vectorise(nu_curr));
            xi_mcmc.row(i_save) = trans(vectorise(xi_curr));
            lambda_mcmc(i_save, 0) = sqrt(lambda2_beta);
            lambda_mcmc(i_save, 1) = lambda_Omega;

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
            Rcpp::Named("Omega") = Omega_mcmc,
            Rcpp::Named("lambda") = lambda_mcmc));
}
