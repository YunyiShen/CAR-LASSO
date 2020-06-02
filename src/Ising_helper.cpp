// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h> // to use sparse matrix
#include "helper.h"
using namespace Rcpp;
using namespace arma;
// FUNCTIONS FOR EXACT SAMPLING //

// Ising Sampling tested 20200601


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
                  int nIter, 
                  const arma::vec& responses, 
                  bool exact){
  // Parameters and results vector:
  int N = graph.n_rows;
  arma::vec state(N, fill::zeros);
  state += NA_REAL;
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
  
  return(exp( H1) / ( exp( H0) + exp(H1) ));// MH ratio here, need changing
}


arma::vec IsingMet(const arma::mat& graph, 
                   const arma::vec& thresholds, 
                   int nIter, 
                   const arma::vec& responses){
  // Parameters and results vector:
  int N = graph.n_rows;
  arma::vec state_temp(N,fill::randu);
  arma::vec state(N,fill::ones);
  state(find(state_temp<0.5)) *= responses(0);
  state(find(state_temp>=0.5)) *= responses(1);
  
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
arma::mat IsingProcess(int nSample, 
                       const arma::mat& graph, 
                       const arma::vec& thresholds, 
                       arma::vec responses)
{
  // Parameters and results vector:
  int N = graph.n_rows;
  arma::vec state_temp(N,fill::randu);
  arma::vec state(N,fill::ones);
  state(find(state_temp<0.5)) *= responses(0);
  state(find(state_temp>=0.5)) *= responses(1);
  
  double u;
  double P;
  
  
  arma::mat Res(N,nSample);
  int node;
  
  // START ALGORITHM
  for (int it = 0 ; it < nSample ; ++it){
    node = floor(R::runif(0,N));
    u = R::runif(0,1);
    P = Pplus(node, graph, state, thresholds, responses);
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
double H(const arma::mat& J, 
         const arma::vec& s, 
         const arma::vec& h)
{
  arma::mat Res;
  Res = (-s.t() * h - .5 * s.t() * J * s);
  return(Res(0,0));
}

double PartitionCpp(
    const arma::mat & graph,
    const arma::mat & thr,
    const arma::vec &responses){
  
  int N = graph.n_rows;
  int n_possible = pow(2,N);
  int t=0;
  double Z=0;
  
  arma::vec temp(N);
  
  for(int i = 0; i<n_possible ; ++i){
    t = i;
    for(int j = 0; j<N ; ++j){
      temp(j) = responses(t % 2);//use binary number coding
      t = t>>1;
    }
    Z+=exp(- H(graph,temp,thr));
  }
  
  return(Z);
  
  
}



// [[Rcpp::export]]
double Ising_PseudoLikelihood_Cpp(const arma::mat & x, 
                                  const arma::mat & graph, 
                                  const arma::vec & thresholds, 
                                  const arma::vec & responses, 
                                  bool logis){
  int k = x.n_rows;
  int n = x.n_cols;
  
  double logPS = 0;
  arma::mat e;
  
  // for every response, compute log PL and sum:
  for (int i = 0 ; i < n ; ++i){
    // for every node, compute log PL and sum:
    for (int j = 0 ; j < k ; ++j){
      e = thresholds(j) + 0.5  * graph.row(j) * x.col(i);
      logPS += x(i,j) * e(0,0)  -  log(exp(responses(0)*e(0,0)) + exp(responses(1)*e(0,0)));
    }
  }
  
  
  return(logis? logPS : exp(logPS));
}





