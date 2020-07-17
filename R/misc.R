stein_loss <- function(Omega_hat, Omega) {
    p <- nrow(Omega)
    Sigma_hat <- solve(Omega_hat)
    loss_prod <- Sigma_hat %*% Omega
    sum(diag(loss_prod)) - log(det(loss_prod)) - p
}