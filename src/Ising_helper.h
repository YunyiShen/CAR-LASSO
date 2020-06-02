#ifndef ISING_HELPER_H
#define ISING_HELPER_H

arma::vec PplusMinMax(int i, 
                      const arma::mat & J, 
                      const arma::vec & s, 
                      const arma::vec & h, 
                      const arma::vec & responses);

arma::vec IsingEx(const arma::mat& graph, 
                  const arma::mat& thresholds, 
                  int nIter, 
                  const arma::vec& responses, 
                  bool exact);

double Pplus(int i, 
             const arma::mat& J, 
             const arma::vec& s, 
             const arma::vec& h, 
             const arma::vec& responses);

arma::vec IsingMet(const arma::mat& graph, 
                   const arma::vec& thresholds, 
                   int nIter, 
                   const arma::vec& responses);


arma::mat IsingProcess(int nSample, 
                       const arma::mat& graph, 
                       const arma::vec& thresholds, 
                       arma::vec responses);

arma::mat IsingSamplerCpp(int n, 
                          const arma::mat& graph, 
                          const arma::vec& thresholds, 
                          int nIter, 
                          const arma::vec responses, 
                          bool exact);

double H(const arma::mat& J, 
         const arma::vec& s, 
         const arma::vec& h);

double PartitionCpp(
    const arma::mat & graph,
    const arma::mat & thr,
    const arma::vec &responses);




double Ising_PseudoLikelihood_Cpp(const arma::mat & x, 
                                  const arma::mat & graph, 
                                  const arma::vec & thresholds, 
                                  const arma::vec & responses, 
                                  bool logis);

#endif
