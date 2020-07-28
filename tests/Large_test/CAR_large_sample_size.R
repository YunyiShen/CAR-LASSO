library(Matrix)
library(RcppArmadillo)
library(Rcpp)
library(RcppProgress)

rm(list = ls())

sourceCpp("./src/CAR-LASSO.cpp")
sourceCpp("./src/graphical-LASSO.cpp")

#ks <- c(10, 20, 30)
#ns <- 100 * 2^(2:7) 
#p <- 2
#sp <- c(0.1,0.2,0.3)

#rep <- 200

ks <- c(10)
ns <- 100 * 2^(2) 
p <- 2
sp <- c(0.1)
rep <- 2

stein_loss_car <- data.frame(matrix(NA, rep, 9))
colnames(stein_loss_car) <- c("loss", "algo", "k", 
                              "n", "sparse", "beta_mean", "beta_sd",
                              "mu_mean", "mu_sd")
stein_loss_glasso <- stein_loss_car

set.seed(42)

for (k in ks) {
    I <- diag(rep(1,k))
    for (n in ns) {
        for (sparse in sp) {
            for (r in 1:rep) {
                cat("k =",k," n =",n, " sp =", sparse, " rep #", r, "\n")
                # generate random graph
                B <- rsparsematrix(k,k,0.2)
                omega <- diag(rgamma(k,5,.1))
                Omega <- t(I-B) %*% omega %*% (I-B)
                Omega <- as.matrix(Omega)
                Sigma <- solve(Omega)

                # generate data
                Design <- 1.0* (matrix(rnorm(n*p,0,1),n,p))
                beta <- matrix(rnorm(p*k,1,3),p,k)


                mu <- rnorm(k,1,1)


                Xbeta <- Design %*% beta

                Z <- matrix(NA,n,k)


                for (i in 1:n) {
                    Z[i,] <- MASS::mvrnorm(1,Sigma %*% (Xbeta[i,]+mu),Sigma)
                }

                cat("Running CAR:\n")
                CAR_test <- CAR_LASSO_Cpp(Z,  Design, n_iter = 25000, 
                          n_burn_in = 5000, thin_by = 10, 
                          r_beta = 1, delta_beta = .01,
                          r_Omega = 1,delta_Omega = .01,
                          progress = T)


                CAR_Graph <- 0 * Omega
                CAR_Graph[upper.tri(CAR_Graph, T)] <- 
                        apply(CAR_test$Omega, 2, mean)
                CAR_Graph <- CAR_Graph + t(CAR_Graph)
                diag(CAR_Graph) <- 0.5 * diag(CAR_Graph)
                CAR_stein_loss <- stein_loss_cpp(CAR_Graph, Omega)

                stein_loss_car[r,] <- c(CAR_stein_loss, "CAR", k ,n , sparse, 
                                        1, 3, 1, 1)
                cat("Stein's loss CAR:" , CAR_stein_loss, "\n\n")
                cat("Running Glasso:\n")
                Glasso <- Graphical_LASSO_Cpp(Z, 25000, 5000, 10, 1, .01, T)

                Glasso_Graph <- 0 * Omega
                Glasso_Graph[upper.tri(Glasso_Graph, T)] <- 
                        apply(Glasso$Omega, 2, mean)
                Glasso_Graph <- Glasso_Graph + t(Glasso_Graph)
                diag(Glasso_Graph) <- 0.5 * diag(Glasso_Graph)

                Glasso_stein_loss <- stein_loss_cpp(Glasso_Graph,Omega)
                stein_loss_glasso[r,] <- c(Glasso_stein_loss, 
                                            "Glasso", k ,n , sparse, 
                                            1, 3, 1, 1)
                cat("Stein's loss Glasso:" , Glasso_stein_loss, "\n\n")
                write.csv(stein_loss_glasso, 
                          "./tests/Large_test/Res/Glasso_sample_size_betamu1311.csv")
                write.csv(stein_loss_car, 
                          "./tests/Large_test/Res/car_sample_size_betamu1311.csv")
            }
        }
    }
}
