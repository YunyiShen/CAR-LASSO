// ((Rcpp::depends(RcppArmadillo)))
#include <RcppArmadillo.h> // to use sparse matrix
#include <tgmath.h>
using namespace Rcpp;
static const double pi = 3.141592653589793238462643383280;

// ((Rcpp::depends(RcppProgress)))
#include <progress.hpp>
#include <progress_bar.hpp>
#include "helper.h"

List resize( const List& x, int n ){
  int oldsize = x.size() ;
  List y(n) ;
  for( int i=0; i<oldsize; i++) y(i) = x(i) ;
  return y ;
}


// Inner function to simulate random uniforms in a matrix:
arma::mat RandMat(int nrow, int ncol)
{
  arma::mat Res = arma::randu(nrow,ncol);
  return(Res);
}

// Computes maximal and minimal probability of node flipping:
arma::vec PplusMinMax(int i, 
                      const arma::mat & J, 
                      const arma::vec & s, 
                      const arma::vec & h, 
                      const arma::vec & responses){
  // The function computes the probability that node i is in Response 1 instead of 0, given all other nodes, which might be missing.
  // Output: minimal and maximal probablity
  
  arma::vec H0(2,fill::ones);
  H0 *= h(i) * responses(0);
  arma::vec H1(2,fill::ones);
  H1 *= h(i) * responses(1);
  
  arma::vec Res(2,fill::zeros);
  
  int N = J.n_rows;
  arma::vec TwoOpts(2,fill::zeros);
  
  for (int it = 0 ; it < N ; ++ it)
  {
    if (i != it)
    {
      if (!R_IsNA(s(it)))
      {
        H0(0) += J(it,i) * responses(0) * s(it);
        H0(1) += J(it,i) * responses(0) * s(it);
        H1(0) += J(it,i) * responses(1) * s(it);
        H1(1) += J(it,i) * responses(1) * s(it); 
      } else 
      {
        
        TwoOpts(0) = J(it,i) * responses(1) * responses(0);
        TwoOpts(1) = J(it,i) * responses(1) * responses(0);
        
        if (TwoOpts(1) > TwoOpts(0))
        {
          H1(0) += TwoOpts(0);
          H1(1) += TwoOpts(1);
          
          H0(0) += J(it,i) * responses(0) * responses(0);
          H0(1) += J(it,i) * responses(0) * responses(1);
        } else 
        {
          H1(0) += TwoOpts(1);
          H1(1) += TwoOpts(0);          
          
          H0(0) += J(it,i) * responses(0) * responses(1);
          H0(1) += J(it,i) * responses(0) * responses(0);
        }
      }
    }
  }
  
  Res(0) = exp(H1(0)) / ( exp(H0(0)) + exp(H1(0)) );
  Res(1) = exp(H1(1)) / ( exp(H0(1)) + exp(H1(1)) );
  
  
  return(Res);
}


// Inner function:
arma::vec IsingEx(const arma::mat& graph, 
                  const arma::mat& thresholds, 
                  double beta, 
                  int nIter, 
                  const arma::vec& responses, 
                  bool exact){
  // Parameters and results vector:
  int N = graph.n_rows;
  arma::vec state(N, fill::zeros);
  state += NA_real;
  double u;
  arma::vec P(2,fill::zeros);
  int maxChain = 100;
  List U(1);
  int minT = 0;
  bool anyNA = true;
  
  do{ 
    // Resize U if needed:
    if (minT > 0)
    {
      U = resize(U, minT+1);
    }
    
    // Generate new random numbers:
    U[minT] = RandMat(nIter, N); // for coupling
    
    // Initialize states:
    for (int i=0; i<N; ++i)
    {
      if (exact)
      {
        state(i) = NA_REAL;
      } else 
      {
        state(i) = R::runif(0,1) < 0.5 ? responses(1) : responses(0);
      }
    }    
    
    // START ALGORITHM CFTP
    for (int t = minT ; t > -1;  t--){
      
      for (int it = 0 ; it < nIter ; it++){
        
        arma::mat Ucur = U[t];
        for (int node = 0 ; node < N ; node++)
        {
          u = Ucur(it, node);
          P = PplusMinMax(node, graph, state, thresholds, responses);
          if (u < P(0)){
            state(node) = responses(1);
          } 
          else if (u >= P(1)){
            state(node) = responses(0);
          } 
          else{
            state(node) = NA_REAL;
          }
        }
      }
    }
    
    anyNA = false;
    if (exact)
    {
      if (minT < maxChain)
      {
        for (int i=0; i<N; i++)
        {
          if (R_IsNA( state(i)))
          {
            anyNA = true;
            break; // no need to check for all states to determin whether there is NA
          }
        } 
      } 
    }    
    minT++;
    
  } while (anyNA);
  
  // Rf_PrintValue(wrap(minT));
  return(state);
}


// FUNCTIONS FOR METROPOLIS SAMPLER //
double Pplus(int i, 
             const arma::mat& J, 
             const arma::vec& s, 
             const arma::vec& h, 
             const arma::vec& responses){
  // The function computes the probability that node i is in Response 1 instead of 0, given all other nodes, which might be missing.
  // Output: minimal and maximal probablity
  int N = J.n_rows;
  double H0 = h(i) * responses(0); // relevant part of the Hamiltonian for state = 0
  double H1 = h(i) * responses(1); // relevant part of the Hamiltonian for state = 1
  
  
  for (int it = 0; it<N; ++it){
    if (i != it)
    {
      H0 += J(it,i) * responses(0) * s(it);
      H1 += J(it,i) * responses(1) * s(it);
    }
  }
  
  return(exp(beta * H1) / ( exp(beta * H0) + exp(beta * H1) ));// MH ratio here, need changing
}


arma::vec IsingMet(const arma::mat& graph, 
                   const arma::vec& thresholds, 
                   int nIter, 
                   const arma::vec& responses){
  // Parameters and results vector:
  int N = graph.n_rows;
  arma::vec state_temp(N,fill::randu);
  arma::vec state(N);
  state(find(state_temp<0.5)) = response(0);
  state(find(state_temp>=0.5)) = response(1);
  
  double u;
  double P;
  
  // START ALGORITHM
  for (int it = 0 ; it < nIter ; ++it){
    for (int node=0 ; node < N ; ++node){
      
      u = R::runif(0,1);
      P = Pplus(node, graph, state, thresholds, responses);
      if (u < P)
      {
        state(node) = responses(1);
      } else 
      {
        state(node) = responses(0);
      } 
      
    }
  }
  
  return(state);
}


///ISING PROCESS SAMPLER: return column vectors
// [[Rcpp::export]]
arma::mat IsingProcess(int nSample, 
                       const arma::mat& graph, 
                       const arma::vec& thresholds, 
                       double beta, 
                       arma::vec responses)
{
  // Parameters and results vector:
  int N = graph.n_rows;
  arma::vec state_temp(N,fill::randu);
  arma::vec state(N);
  state(find(state_temp<0.5)) = response(0);
  state(find(state_temp>=0.5)) = response(1);
  
  double u;
  double P;
  
  
  arma::mat Res(N,nSample);
  int node;
  
  // START ALGORITHM
  for (int it = 0 ; it < nSample ; ++it){
    node = floor(R::runif(0,N));
    u = R::runif(0,1);
    P = Pplus(node, graph, state, thresholds, beta, responses);
    if (u < P)
    {
      state(node) = responses(1);
    } else 
    {
      state(node) = responses(0);
    }
    
    Res.col(it) = state;
  }
  
  return(Res);
}

// OVERAL FUNCTION //
// [[Rcpp::export]]
arma::mat IsingSamplerCpp(int n, 
                          const arma::mat& graph, 
                          const arma::vec& thresholds, 
                          int nIter, 
                          const arma::vec responses, 
                          bool exact){
  int Ni = graph.n_rows; // number of nodes
  arma::mat Res(Ni,n,fill::zeros);
  Res += NA_REAL;
  arma::vec state(Ni);
  if (exact){
    for (int s = 0 ; s < n ; ++s){
      state = IsingEx(graph, thresholds, nIter, responses, exact);
      Res.col(s) = state;
    }
  } 
  else{
    for (int s = 0 ; s < n ; ++s){
      state = IsingMet(graph, thresholds, nIter, responses);
      Res.col(s) = state;
    }
  }
  
  return(Res); //give column vectors
}

// Hamiltonian:
// [[Rcpp::export]]
double H(const arma::mat& J, 
         const arma::vec& s, 
         const arma::vec& h)
{
  double Res = 0;
  Res = -s.t() * h - .5 * s.t() * J * s;
  return(Res);
}



// log MH ratio 

double logMH_Ising(const arma::mat & data,
                   const arma::mat & thr_prop,
                   const arma::mat & thr_curr,
                   const arma::mat & J_prop,
                   const arma::mat & J_curr,
                   const arma::vec & responses,
                   double log_prior_beta_curr,
                   double log_prior_J_curr,
                   double log_prior_beta_prop,
                   double log_prior_J_prop,             
                   bool exact,
                   int nIter,
                   int n,int k){
  
  arma::vec Z_prime_i;
  double log_q_theta_curr_data = 0;
  double log_q_theta_prop_data = 0;
  double log_q_theta_curr_Z_prime = 0;
  double log_q_theta_prop_Z_prime = 0;
  
  for(int i = 0 ; i < n ; ++i){
    Z_prime_i = IsingSamplerCpp(1, J_prop, thr_prop.col(i), 
                                nIter, responses, exact)
    
    log_q_theta_curr_data -= H(J_curr,data.col(i),thr_curr.col(i));
    log_q_theta_prop_data -= H(J_prop,data.col(i),thr_prop.col(i));
    
    log_q_theta_curr_Z_prime -= H(J_curr,Z_prime_i,thr_curr.col(i));
    log_q_prop_curr_Z_prime -= H(J_prop,Z_prime_i,thr_prop.col(i));
  }
  
  double log_MH = log_prior_beta_prop + 
    log_prior_J_prop + 
    log_q_theta_prop_data + 
    log_q_theta_curr_Z_prime -
    log_prior_beta_curr - 
    log_prior_J_curr - 
    log_q_theta_curr_data - 
    log_q_theta_prop_Z_prime;
  
  return(log_MH);
  
}



Rcpp::List Ising_LASSO_Cpp(const arma::mat & data,
                           const arma::mat & design,
                           const int n_iter, // how many iteractions?
                           const int n_burn_in, // burn in
                           const int thin_by, // thinning?
                           const double r_beta, // prior on lambda of beta
                           const double delta_beta,
                           const double r_B,
                           const double delta_B,
                           const arma::vec & propsd_mu,
                           const arma::vec & propsd_beta,
                           const arma::vec & propsd_J,
                           const arma::vec & propsd_lambda,// beta first then J
                           bool exact, // whether use CFTP for the sampler in sampling Z*
                           int nIter, // still, interactions in Ising sampling 
                           bool progress,
                           bool verbos,int reportby){
  arma::vec responses(2);
  responses(0) = min(min(data));
  responses(1) = max(max(data));
  int k = data.n_rows; // number of nodes
  int p = design.n_cols; //number of predictors
  int n = data.n_cols; // number of samples
  
  int accept_Ising = 0;
  double accept_lambda = 0;
  
  int n_save = floor(n_iter/thin_by); //
  int i_save = 0;  
  
  // mcmc matrices:
  arma::mat beta_mcmc(n_save,k * p); // beta mcmc
  beta_mcmc += NA_REAL; // initial with NAs
  
  arma::mat J_mcmc(n_save,floor(.5*(k-1)*k));// adjacency matrix
  J_mcmc += NA_REAL;
  
  arma::mat mu_mcmc(n_save , k); // mean for node 1 to k
  mu += NA_REAL;
  
  arma::mat lambda(n_save , 2); // LASSO parameter for beta and J
  lambda += NA_REAL;
  
  // some initial values
  arma::mat J_mat_curr(k,k,fill::zeros);
  arma::uvec upperdiag = trimatu_ind(size(J_curr),1);
  
  arma::vec beta_curr(k*p,fill::zeros);
  arma::mat betamat_curr = reshape(beta_curr,p,k);
  arma::vec mu_curr(k,fill::zeros);
  arma::vec J_curr(upperdiag.n_elem,fill::zeros);
  
  
  arma::vec beta_prop;
  arma::mat beta_mat_prop;
  arma::vec J_prop;
  arma::mat J_mat_prop;
  
  double lambda_beta_curr = 1;
  double lambda_J_curr = 1;
  
  double lambda_beta_prop;
  double lambda_J_prop;
  
  double log_prior_beta_curr = dLaplace_Cpp(beta_curr,0,lambda_beta_curr,true);
  double log_prior_J_curr = dLaplace_Cpp(J_curr,0,lambda_beta_curr,true);
  
  double log_prior_beta_prop;
  double log_prior_J_prop;
  
  double logMH;
  
  
  double logPostlambda_beta_curr = sum(dLaplace_Cpp(beta_curr),0,lambda_beta_curr,true) + 
    R::dgamma(lambda_beta_curr,r_beta,1/delta_beta,true);
  double logPostlambda_J_curr = sum(dLaplace_Cpp(J_curr),0,lambda_J_curr,true) + 
    R::dgamma(lambda_J_curr,r_J,1/delta_J,true);
  
  double logPostlambda_beta_prop;
  double logPostlambda_J_curr;
  
  
  arma::mat thr_curr(size(data),fill::zeros);
  arma::mat thr_prop = thr_curr;
  
  Progress p((n_iter+n_burn_in), progress); 
  
  // main loop
  for(int i = 0 ; i < (n_iter+n_burn_in) ; ++i){
    
    if (Progress::check_abort()){
      Rcerr << "keyboard abort\n";
      return(Rcpp::List::create(
          Rcpp::Named("beta") = beta_mcmc,
          Rcpp::Named("mu") = mu_mcmc,
          Rcpp::Named("J") = J_mcmc,
          Rcpp::Named("lambda") = lambda_mcmc
      );
    }
    
    // update Ising parameters
    beta_prop = propsd_beta % randn(size(beta_curr)) + beta_curr;
    beta_mat_prop = reshape(beta_prop,p,k);
    J_prop = propsd_J % randn(size(J_curr)) + J_curr;
    J_mat_prop = J_mat_prop.zeros();
    J_mat_prop(upperdiag) = J_prop;
    J_mat_prop += J_mat_prop.t();
    mu_prop = propsd_mu % randn(size(mu_curr)) + mu_curr;
    
    thr_prop = trans(design * beta_mat_prop);
    thr_prop.each_col() += mu_prop;
    
    // MH ratio
    
    logMH = logMH_Ising(data,
                        thr_prop,thr_curr,
                        J_mat_prop,J_mat_curr,
                        responses,
                        log_prior_beta_curr,log_prior_J_curr,
                        log_prior_beta_prop,log_prior_J_prop,             
                        exact,nIter,n,k);
    
    if(log(R::runif(0,1)) <= logMH){ // accept
      beta_curr = beta_prop;
      beta_mat_curr = beta_mat_prop;
      
      J_curr = J_prop;
      J_mat_curr = J_mat_prop;
      
      mu_curr = mu_prop;
      
      thr_curr = thr_prop;
      
      log_prior_beta_curr = log_prior_beta_prop;
      
      log_prior_J_curr = log_prior_beta_prop;
      
      accept_Ising++;
      
    }
    
    // updating lambdas
    lambda_beta_prop = lambda_beta_curr + propsd_lambda(0) * R::rnorm(0,1);
    lambda_J_prop = lambda_J_curr + propsd_lambda(1) * R::rnorm(0,1);
    
    if(lambda_beta_prop>0){
      logPostlambda_beta_prop = sum(dLaplace_Cpp(beta_curr),0,lambda_beta_prop,true) + 
        R::dgamma(lambda_beta_prop,r_beta,1/delta_beta,true);
      if(log(R::runif(0,1)) <= logPostlambda_beta_prop-logPostlambda_beta_curr){
        lambda_beta_curr = lambda_beta_prop;
        logPostlambda_beta_curr = logPostlambda_beta_prop;
        accept_lambda += .5 ;
      }
    }
    
    if(lambda_J_prop>0){
      logPostlambda_J_prop = sum(dLaplace_Cpp(J_curr),0,lambda_J_prop,true) + 
        R::dgamma(lambda_J_prop,r_J,1/delta_J,true);
      if(log(R::runif(0,1)) <= logPostlambda_J_prop-logPostlambda_J_curr){
        lambda_J_curr = lambda_J_prop;
        logPostlambda_J_curr = logPostlambda_J_prop;
        accept_lambda += .5 ;
      }
    }
    
    if(verbos && (i+1) % reportby ==0){
      Rcout << "Grand acceptance rate for Ising:" << 100 * accept_Ising/(i+1) << "%\n" <<endl;
      Rcout << "Grand acceptance rate for LASSO:" << 100 * accept_lambda/(i+1) << "%\n" <<endl;
    }
    
    // 
    
    if( (i-burn_in)>=0 && (i+1-burn_in)%thin_by == 0 ){
      lambda_mcmc(i_save,0) = lambda_beta_curr;
      lambda_mcmc(i_save,1) = lambda_J_curr;
      
      beta_mcmc.row(i) = beta_curr.t();
      J_mcmc.row(i) = J_curr.t();

      i_save++ ;
    }
    
    
    
    p.increment(); // progress bar
  }
  
  return(Rcpp::List::create(
      Rcpp::Named("beta") = beta_mcmc,
      Rcpp::Named("mu") = mu_mcmc,
      Rcpp::Named("J") = J_mcmc,
      Rcpp::Named("lambda") = lambda_mcmc
  );
  
  
}
