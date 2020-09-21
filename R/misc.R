stein_loss <- function(Omega, Omega_hat) {
    stein_loss_cpp(Omega, Omega_hat)
}

get_graph <- function(CAR_sample,k, summary = "mean"){
  Omega <- matrix(0,k,k)
  Omega[upper.tri(Omega,T)] = apply(CAR_sample$Omega,2,summary)
  Omega <- Omega+t(Omega)
  diag(Omega) <- 0.5 * diag(Omega)
  return(Omega)
}
