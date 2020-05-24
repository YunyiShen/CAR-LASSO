// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h> // to use sparse matrix
#include <tgmath.h>
using namespace Rcpp;
using namespace arma;
#define pi = 3.141592653589793238462643383280;

/* Convention:
 *    single sample will be saved as col vectors
 */

// Auto-Poisson

arma::vec Auto_Poisson_Gibbs_Single_Cpp(const arma::sp_mat &graph, // graph
                                        const arma::vec &thresholds, // mean 
                                        const int & Winsorized, // right censor
                                        const int & nIter){  // number of Gibbs iters
  int N = graph.n_rows;
  int temp;
  arma::vec Res(N);
  
  for (int i = 0 ; i < N ; ++i){
    Res(i) = R::rpois(thresholds(i)); // initializing 
  }
  
  for(int i = 0 ; i < nIter ; ++i){ // outer loop, for number of samples
    for(int j = 0 ; j < N ; ++j){
      temp = R::rpois(exp(thresholds(j) + 0.5 * (Res.t() * graph.col(j)))); // Gibbs scan
      Res(i) = temp > Winsorized ? Winsorized : temp; // right censored 
    }
  }
  return(Res);
} 


// [[Rcpp::export]]
arma::mat Auto_Poisson_Gibbs_Batch_Cpp(const int & Nsample, // number of samples
                                       const arma::sp_mat &graph, // graph
                                       const arma::vec &thresholds, // mean 
                                       const int & Winsorized, // right censor
                                       const int & nIter){  
  int N = graph.n_rows;
  arma::mat Res(N,Nsample)
  for(int i = 0 ; i < Nsample ; ++i ){
    Res.col(i) = Auto_Poisson_Gibbs_Single_Cpp(graph,thresholds,Winsorized,nIter);
  }
  return(Res);
}

// [[Rcpp::export]]
double Auto_Poisson_Pseudo_likelihood_Cpp(const arma::mat Sample,
                                          const arma::sp_mat & graph,
                                          const arma::vec & thresholds){
  int Nsample = Sample.n_cols; //samples
  int N = Sample.n_rows; // nodes
  double Res = 0;
  arma::mat Interactions = graph * Sample ; // interaction terms
  for(int i = 0 ; i < Nsample ; ++i){
    arma::vec Lambda = exp(thresholds + Interactions.col(i)); // lambda of Poisson of that repeat
    for(int j = 0 ; j < N ; ++j){
      Res += R::dpois(Sample(j,i),Lambda(j),true); // likelihood, use Poisson approximation 
    } 
  }
  
  return(Res);
}

// Q-function for auto poisson 
double Auto_Poisson_Q_Cpp(const arma::mat Sample,
                          const arma::sp_mat & graph,
                          const arma::vec & thresholds){
  // not sure if this line work sonce sum(...) will give a row vec
  arma::rowvec thr_contribution = thresholds.t() * Sample - sum( lgamma(Sample + 1),0); // non-interaction part of Besag's Q-function, see Augustin et al. 2006
  double Res = 0;
  int Nsample = Sample.n_cols;
  
  for(int i = 0 ; i < Nsample ; ++i){
    Res += thr_contribution(i) + .5 * (Sample.col(i).t() * graph * Sample.col(i)); // Besag's Q function for (spatial) auto Poisson model 
  }
  return(Res);
}


// logLik of a normal response using composition as predictor
double logLik_Normal_Response(const arma::vec & data,
                              const arma::mat & Design,
                              const arma::vec & beta_fix,
                              const arma::mat & composition,
                              const arma::vec & beta_composition,
                              const double & sigma){
  int N = Design.n_rows;
  arma::vec mu = Design * beta_fix + composition.t() * beta_composition;
  double Res = 0;
  for(int i = 0 ; i < N ; ++i){
    Res += R::dnorm(data(i) , mu(i) , true);
  }
  return(Res);
}


// Laplace distribution density, for L1 regularization
double dLaplace_Cpp(const double &x,const double &mu, const double &b, bool log){
  double logd = - abs(x-mu)/b - log(2*b);
  return(log ? logd : exp(logd));
}



/* TODO:
 * 1. Coupling from the past for one perfect sample from the auto poisson model
 * 2. Auto-NegBin model sampling, 
 * 3. derive the Q-function 
 */



// Auto-normal

/* Starting here is auto-normal, the main idea was fitting a SAR and convert it to CAR if needed, see Ver Hoef et al. 2018
 * The Reason for using SAR rather than CAR is that SAR asks for a ralative easier to check condition 
 *   condition on B (I-B) non-singular, while CAR asks for C to have positive eigen value
 *   We can use independent Laplace prior in SAR to get similar result as LASSO.
 * This is all because Sigma has to be positive semi-definite.
 * Also, matrix B can be asymmertic, 
 *   however this will cause some explaning problem, 
 *    as B is no longer an adjacency matrix of the faithfuil graph related with the joint distribution 
 */



double Loglik_SAR_Cpp(const arma::vec & Z,
                      const arma::vec & mu,
                      const arma::vec & sigma,
                      const arma::sp_mat B){
  int N = Z.n_elem;
  arma::sp_mat I_minus_B = - B;
  I_minus_B.diag() += 1; // I - B
  arma::sp_mat Sigma(N,N);
  Simga.diag() = sigma; // Var for noise
 
 // calculating normalizing constant
 // The cov matrix of SAR is : S = (I-B)^{-1} \Sigma (I^T-B^T)^{-1}
 // Thus log(det(S)) = log( det(\Sigma) ) - 2 log (det(I-B))
  double logC = - (N/2) * log (2 * pi) - // 2pi
    .5 * (sum(log(sigma))) + // Sigma
    log(det(I_minus_B)); // I-B part
  
  arma::vec try_solve;
  spsolve(try_solve, I_minus_B, Z);
  if(!try_solve(0)){
    return(Rcpp::R_NegInf); // if I-B not inversible, return -Inf 
  }
  
  double logH = - 0.5 * X.t() * spsolve(I_minus_B , Sigma * try_solve); // Hamiltonian
  return(logC+logH);
}

Rcpp::List Sigma_to_CAR_Cpp(const arma::mat & Sigma){
  
  arma::sp_mat Q = solve(Sigma); // percision matrix
  arma::mat R = -Sigma; // R matrix, see Ver Hoef et al. 2018
  diag(R) = 0;
  arma::vec M = 1/diag(Q);
  arma::sp_mat invD(size(Q));
  diag(invD) = M ;
  return(Rcpp::List::create(
    Rcpp::Named("C") = invD * R, // C matrix per Ver Hoef et al. 2018
    Rcpp::Named("M") = M // diag of M matrix per Ver Hoef et al. 2018
  ));
}

arma::mat Sample_SAR_Cpp(const int Nsample,
                         const arma::vec & mu,
                         const arma::vec & sigma,
                         const arma::sp_mat B){
  int N = mu.n_elem;
  arma::sp_mat I_minus_B = - B;
  I_minus_B.diag() += 1; // I - B
  arma::sp_mat Sigma(N,N);
  Simga.diag() = sigma; // Var for noise
  arma::mat try_solve;
  spsolve(try_solve, I_minus_B);// This will be all false if I_minus_B is singular
  if(!try_solve(0,0)){
    return(try_solve); // if I-B not inversible, return 0s
  }
  
  arma::mat full_Sigma = try_solve * Sigma * try_solve.t(); //cov matrix of this Sar model
  arma::mat Res;
  mvnrnd(Res, mu, full_Sigma, Nsample )); // this will give column vectors
  return(Res);
}



