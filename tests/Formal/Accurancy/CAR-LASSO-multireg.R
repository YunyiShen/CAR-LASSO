library(Matrix)
library(RcppArmadillo)
library(Rcpp)
library(RcppProgress)
library(bayesm)

sourceCpp("./src/SRG-LASSO.cpp")
sourceCpp("./src/CAR-LASSO.cpp")
sourceCpp("./src/graphical-LASSO.cpp")
sourceCpp("./src/CAR-ALASSO.cpp")
source("./tests/Formal/Accurancy/misc.R")
source("./tests/Formal/Accurancy/Graph_generator.R")
settings_link <- "./tests/Formal/Accurancy/settings/"
res_loss_file <- "./tests/Formal/Accurancy/results/loss_CAR.csv"
res_graph_Omega_file <- "./tests/Formal/Accurancy/results/graph_Omega_CAR.csv"
res_graph_beta_file <- "./tests/Formal/Accurancy/results/graph_beta_CAR.csv"


ks <- c(100,30)
ss <- c(.2,.5)
ps <- c(10,5)

nrep <- 50
sigma_models <- 1:6

thr <- c( 1e-2, 1e-3 )

n_res_loss <- length(ks) * length(ss) * 
              length(ps) * length(sigma_models) * 
              nrep * 4
res_loss <- data.frame(matrix(NA,nrow = n_res_loss,ncol = 9))
colnames(res_loss) <- c("k","p","n","s","mod",
                        "irep","algo","logL2beta",
                        "steinOmega")

res_graph_Omega <- data.frame(matrix(NA,nrow = n_res_loss,ncol = 19))
colnames(res_graph_Omega) <- c("k","p","n","s","mod",
                        "irep","algo","TP_bayes","TN_bayes",
                        "FP_bayes","FN_bayes","TP_thrh","TN_thrh",
                        "FP_thrh","FN_thrh","TP_thrl","TN_thrl",
                        "FP_thrl","FN_thrl")

res_graph_beta <- data.frame(matrix(NA,nrow = n_res_loss,ncol = 19))
colnames(res_graph_beta) <- c("k","p","n","s","mod",
                        "irep","algo","TP_bayes","TN_bayes",
                        "FP_bayes","FN_bayes","TP_thrh","TN_thrh",
                        "FP_thrh","FN_thrh","TP_thrl","TN_thrl",
                        "FP_thrl","FN_thrl")

i_res_loss <- 1
i_res_graph_Omega <- 1
i_res_graph_beta <- 1

for(k in ks) {
    n <- 50 * (k==30) + 200 * (k==100)
    for(p in ps){
        for(s in ss){

            beta <- read.csv(paste0(
                settings_link, "beta_k", k, "_p", p , "_s",s, ".csv"
            ))
            beta <- as.matrix(beta)
            Design <- read.csv(paste0(
                settings_link, "design_p",p,"_n",n,".csv"
            ))
            Design <- as.matrix(Design)

            for(mod in sigma_models){
                graph_tmp <- do.call(paste0("g_model",mod),list(k=k))
                Omega <- graph_tmp$Omega
                Sigma <- graph_tmp$Sigma
                
                for(i in 1:nrep){
                    cat("k =",k,",p =",p,",s =",s,",mod =",mod,",rep =",i,"\n")
                    Z <- makedata(Sigma,Design,beta,n)

                    # CAR
                    cat("CAR:\n")
                    sample_CAR <- CAR_LASSO_Cpp(Z,  Design, n_iter = 10000, 
                          n_burn_in = 5000, thin_by = 10, 
                          r_beta = 1, delta_beta = .01,
                          r_Omega = 1,delta_Omega = .01,
                          progress = T)
                    
                    Omega_CAR <- get_graph(sample_CAR,k)
                    beta_CAR <- get_beta(sample_CAR,p)

                    res_loss[i_res_loss,1:6] <- c(k,p,n,s,mod,i)
                    res_loss$algo[i_res_loss] <- "CAR-LASSO"
                    res_loss$logL2beta[i_res_loss] <- log_l2(beta,beta_CAR)
                    res_loss$steinOmega[i_res_loss] <- stein_loss_cpp(Omega_CAR, Omega)
                    write.csv(res_loss,res_loss_file)
                    i_res_loss <- i_res_loss + 1

                    # A-CAR
                    cat("CAR-A:\n")
                    sample_CAR_A <- CAR_ALASSO_Cpp(Z,  Design, n_iter = 10000, 
                          n_burn_in = 5000, thin_by = 10, 
                          r_beta = 1+0*beta, delta_beta = .01 + 0 * beta,
                          r_Omega = rep(1,.5*(k+1)*k),
                          delta_Omega = rep(.01,.5*(k+1)*k),
                          progress = T)
                    
                    Omega_CAR_A <- get_graph(sample_CAR_A,k)
                    beta_CAR_A <- get_beta(sample_CAR_A,p)

                    res_loss[i_res_loss,1:6] <- c(k,p,n,s,mod,i)
                    res_loss$algo[i_res_loss] <- "CAR-ALASSO"
                    res_loss$logL2beta[i_res_loss] <- log_l2(beta,beta_CAR_A)
                    res_loss$steinOmega[i_res_loss] <- stein_loss_cpp(Omega_CAR_A, Omega)
                    write.csv(res_loss,res_loss_file)
                    i_res_loss <- i_res_loss + 1

                    # SRG
                    cat("SRG:\n")
                    sample_SRG <- SRG_LASSO_Cpp(Z,  Design, n_iter = 10000, 
                          n_burn_in = 5000, thin_by = 10, 
                          r_beta = 1, delta_beta = .01,
                          r_Omega = 1,delta_Omega = .01,
                          progress = T)
                    
                    Omega_SRG <- get_graph(sample_SRG,k)
                    beta_SRG <- get_beta(sample_SRG,p)

                    res_loss[i_res_loss,1:6] <- c(k,p,n,s,mod,i)
                    res_loss$algo[i_res_loss] <- "SRG-LASSO"
                    res_loss$logL2beta[i_res_loss] <- log_l2(beta %*% Omega,beta_SRG)
                    res_loss$steinOmega[i_res_loss] <- stein_loss_cpp(Omega_SRG, Omega)
                    write.csv(res_loss,res_loss_file)
                    i_res_loss <- i_res_loss + 1

                    # multireg
                    cat("multireg:\n")
                    sample_multireg <- multireg_Sample(Z,Design,k,p)
                    res_loss[i_res_loss,1:6] <- c(k,p,n,s,mod,i)
                    res_loss$algo[i_res_loss] <- "multireg"
                    res_loss$logL2beta[i_res_loss] <- log_l2(beta,sample_multireg$beta)
                    res_loss$steinOmega[i_res_loss] <- stein_loss_cpp(sample_multireg$Omega, Omega)
                    write.csv(res_loss,res_loss_file)
                    i_res_loss <- i_res_loss + 1

                    # now evaluate the graph_Omega
                    ## CAR-A
                    res_graph_Omega[i_res_graph_Omega,1:6] <- c(k,p,n,s,mod,i)
                    res_graph_Omega$algo[i_res_graph_Omega] <- "CAR-ALASSO"
                    res_graph_Omega[i_res_graph_Omega,8:11] <- 
                        get_counts_Omega(abs(Omega_CAR_A/sample_multireg$Omega)>.5,Omega!=0)
                    res_graph_Omega[i_res_graph_Omega,1:4+11] <- 
                        get_counts_Omega(abs(Omega_CAR_A)>thr[1],Omega!=0)
                    res_graph_Omega[i_res_graph_Omega,1:4+15] <- 
                        get_counts_Omega(abs(Omega_CAR_A)>thr[2],Omega!=0)
                    write.csv(res_graph_Omega,res_graph_Omega_file)
                    i_res_graph_Omega <- i_res_graph_Omega + 1

                    ## CAR
                    res_graph_Omega[i_res_graph_Omega,1:6] <- c(k,p,n,s,mod,i)
                    res_graph_Omega$algo[i_res_graph_Omega] <- "CAR-LASSO"
                    res_graph_Omega[i_res_graph_Omega,8:11] <- 
                        get_counts_Omega(abs(Omega_CAR/sample_multireg$Omega)>.5,Omega!=0)
                    res_graph_Omega[i_res_graph_Omega,1:4+11] <- 
                        get_counts_Omega(abs(Omega_CAR)>thr[1],Omega!=0)
                    res_graph_Omega[i_res_graph_Omega,1:4+15] <- 
                        get_counts_Omega(abs(Omega_CAR)>thr[2],Omega!=0)
                    write.csv(res_graph_Omega,res_graph_Omega_file)
                    i_res_graph_Omega <- i_res_graph_Omega + 1

                    ## SRG
                    res_graph_Omega[i_res_graph_Omega,1:6] <- c(k,p,n,s,mod,i)
                    res_graph_Omega$algo[i_res_graph_Omega] <- "SRG-LASSO"
                    res_graph_Omega[i_res_graph_Omega,8:11] <- 
                        get_counts_Omega(abs(Omega_SRG/sample_multireg$Omega)>.5,Omega!=0)
                    res_graph_Omega[i_res_graph_Omega,1:4+11] <- 
                        get_counts_Omega(abs(Omega_SRG)>thr[1],Omega!=0)
                    res_graph_Omega[i_res_graph_Omega,1:4+15] <- 
                        get_counts_Omega(abs(Omega_SRG)>thr[2],Omega!=0)
                    write.csv(res_graph_Omega,res_graph_Omega_file)
                    i_res_graph_Omega <- i_res_graph_Omega + 1

                    ## multireg
                    res_graph_Omega[i_res_graph_Omega,1:6] <- c(k,p,n,s,mod,i)
                    res_graph_Omega$algo[i_res_graph_Omega] <- "multireg"
                    res_graph_Omega[i_res_graph_Omega,1:4+11] <- 
                        get_counts_Omega(abs(sample_multireg$Omega)>thr[1],Omega!=0)
                    res_graph_Omega[i_res_graph_Omega,1:4+15] <- 
                        get_counts_Omega(abs(sample_multireg$Omega)>thr[2],Omega!=0)
                    write.csv(res_graph_Omega,res_graph_Omega_file)
                    i_res_graph_Omega <- i_res_graph_Omega + 1


                    # graph beta
                    ## CAR-A
                    res_graph_beta[i_res_graph_beta,1:6] <- c(k,p,n,s,mod,i)
                    res_graph_beta$algo[i_res_graph_beta] <- "CAR-ALASSO"
                    res_graph_beta[i_res_graph_beta,8:11] <- 
                        get_counts_beta(abs(beta_CAR_A/sample_multireg$beta)>.5,beta!=0)
                    res_graph_beta[i_res_graph_beta,1:4+11] <- 
                        get_counts_beta(abs(beta_CAR_A)>thr[1],beta!=0)
                    res_graph_beta[i_res_graph_beta,1:4+15] <- 
                        get_counts_beta(abs(beta_CAR_A)>thr[2],beta!=0)
                    write.csv(res_graph_beta,res_graph_beta_file)
                    i_res_graph_beta <- i_res_graph_beta + 1

                    ## CAR
                    res_graph_beta[i_res_graph_beta,1:6] <- c(k,p,n,s,mod,i)
                    res_graph_beta$algo[i_res_graph_beta] <- "CAR-LASSO"
                    res_graph_beta[i_res_graph_beta,8:11] <- 
                        get_counts_beta(abs(beta_CAR/sample_multireg$beta)>.5,beta!=0)
                    res_graph_beta[i_res_graph_beta,1:4+11] <- 
                        get_counts_beta(abs(beta_CAR)>thr[1],beta!=0)
                    res_graph_beta[i_res_graph_beta,1:4+15] <- 
                        get_counts_beta(abs(beta_CAR)>thr[2],beta!=0)
                    write.csv(res_graph_beta,res_graph_beta_file)
                    i_res_graph_beta <- i_res_graph_beta + 1

                    ## SRG
                    res_graph_beta[i_res_graph_beta,1:6] <- c(k,p,n,s,mod,i)
                    res_graph_beta$algo[i_res_graph_beta] <- "SRG-LASSO"
                    res_graph_beta[i_res_graph_beta,8:11] <- 
                        get_counts_beta(abs((beta_SRG%*%Omega_SRG)/sample_multireg$beta)>.5,beta!=0)
                    res_graph_beta[i_res_graph_beta,1:4+11] <- 
                        get_counts_beta(abs(beta_SRG%*%Omega_SRG)>thr[1],beta!=0)
                    res_graph_beta[i_res_graph_beta,1:4+15] <- 
                        get_counts_beta(abs(beta_SRG%*%Omega_SRG)>thr[2],beta!=0)
                    write.csv(res_graph_beta,res_graph_beta_file)
                    i_res_graph_beta <- i_res_graph_beta + 1

                    ## multireg
                    res_graph_beta[i_res_graph_beta,1:6] <- c(k,p,n,s,mod,i)
                    res_graph_beta$algo[i_res_graph_beta] <- "multireg"
                    res_graph_beta[i_res_graph_beta,1:4+11] <- 
                        get_counts_beta(abs(sample_multireg$beta)>thr[1],beta!=0)
                    res_graph_beta[i_res_graph_beta,1:4+15] <- 
                        get_counts_beta(abs(sample_multireg$beta)>thr[2],beta!=0)
                    write.csv(res_graph_beta,res_graph_beta_file)
                    i_res_graph_beta <- i_res_graph_beta + 1

                }
            }   
        }
    }
}
