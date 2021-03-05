library(Matrix)
library(RcppArmadillo)
library(Rcpp)
library(RcppProgress)
library(bayesm)


sourceCpp("./src/graphical-LASSO.cpp")
sourceCpp("./src/graphical-ALASSO.cpp")
source("./tests/Formal/Accurancy/misc.R")
source("./tests/Formal/Accurancy/Graph_generator.R")
settings_link <- "./tests/Formal/Accurancy/settings/"
res_loss_file <- "./tests/Formal/Accurancy/results/loss_Glasso.csv"
res_graph_Omega_file <- "./tests/Formal/Accurancy/results/graph_Omega_Glasso.csv"


ks <- c(100,30)
#ks <- c(30)
ss <- c(.2,.5)
#ss <- c(.2)
ps <- c(10,5)
#ps <- c(5)

nrep <- 50
sigma_models <- 1:6

retry <- 5

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

                    # GLASSO
                    cat("GLASSO:\n")
                    for(jj in 1:retry){
                        sample_GL <- tryCatch( Graphical_LASSO_Cpp(Z, 
                                        10000, 5000, 10, 1, .01, T),
                                error = function(e){return(list())} )
                        if(length(sample_GL)>0) break
                        else cat("retry ",jj,"\n")
                    }
                    Omega_GL <- get_graph(sample_GL,k)

                    res_loss[i_res_loss,1:6] <- c(k,p,n,s,mod,i)
                    res_loss$algo[i_res_loss] <- "GLASSO"
                    
                    res_loss$steinOmega[i_res_loss] <- stein_loss_cpp(Omega_GL, Omega)
                    write.csv(res_loss,res_loss_file)
                    i_res_loss <- i_res_loss + 1

                    # A-GLASSO
                    cat("GLASSO-A:\n")
                    for(jj in 1:retry){
                        sample_GL_A <- tryCatch( Graphical_ALASSO_Cpp(Z,  n_iter = 10000, 
                          n_burn_in = 5000, thin_by = 10, 
                          lambda_a = rep(1,.5*(k-1)*k),
                          lambda_b = rep(.01,.5*(k-1)*k),
                          progress = T),
                          error = function(e){return(list())} )
                        if(length(sample_GL_A)>0) break
                        else cat("retry ",jj,"\n")
                    }
                    
                    Omega_GL_A <- get_graph(sample_GL_A,k)

                    res_loss[i_res_loss,1:6] <- c(k,p,n,s,mod,i)
                    res_loss$algo[i_res_loss] <- "GALASSO"
                    
                    res_loss$steinOmega[i_res_loss] <- stein_loss_cpp(Omega_GL_A, Omega)
                    write.csv(res_loss,res_loss_file)
                    i_res_loss <- i_res_loss + 1

                    # ad hoc
                    cat("ad-hoc:\n")
                    for(jj in 1:retry){
                        sample_adhoc <- tryCatch( MASS::ginv(cov(Z)),
                            error = function(e){return(list())})
                        if(length(sample_adhoc)>0) break
                        else cat("retry ",jj,"\n")
                    }
                    if(length(sample_adhoc)>0){
                        Omega_adhoc <- as.matrix(sample_adhoc)
                        

                        res_loss[i_res_loss,1:6] <- c(k,p,n,s,mod,i)
                        res_loss$algo[i_res_loss] <- "ad-hoc"
                        
                        res_loss$steinOmega[i_res_loss] <- stein_loss_cpp(Omega_adhoc, Omega)
                    }
                    write.csv(res_loss,res_loss_file)
                    i_res_loss <- i_res_loss + 1

                    # multireg
                    cat("multireg:\n")
                    for(jj in 1:retry){
                        sample_multireg <- tryCatch( multireg_mu0_Sample(Z,Design,k,p),
                                            error = function(e){return(list())} )
                        if(length(sample_multireg)>0) break
                        else cat("retry ",jj,"\n")
                    }
                    res_loss[i_res_loss,1:6] <- c(k,p,n,s,mod,i)
                    res_loss$algo[i_res_loss] <- "multireg_mu0"
                    
                    res_loss$steinOmega[i_res_loss] <- stein_loss_cpp(sample_multireg$Omega, Omega)
                    write.csv(res_loss,res_loss_file)
                    i_res_loss <- i_res_loss + 1

                    # now evaluate the graph_Omega
                    ## GLASSO-A
                    res_graph_Omega[i_res_graph_Omega,1:6] <- c(k,p,n,s,mod,i)
                    res_graph_Omega$algo[i_res_graph_Omega] <- "GALASSO"
                    res_graph_Omega[i_res_graph_Omega,8:11] <- 
                        get_counts_Omega(abs(Omega_GL_A/sample_multireg$Omega)>.5,Omega!=0)
                    res_graph_Omega[i_res_graph_Omega,1:4+11] <- 
                        get_counts_Omega(abs(Omega_GL_A)>thr[1],Omega!=0)
                    res_graph_Omega[i_res_graph_Omega,1:4+15] <- 
                        get_counts_Omega(abs(Omega_GL_A)>thr[2],Omega!=0)
                    write.csv(res_graph_Omega,res_graph_Omega_file)
                    i_res_graph_Omega <- i_res_graph_Omega + 1

                    ## GLASSO
                    res_graph_Omega[i_res_graph_Omega,1:6] <- c(k,p,n,s,mod,i)
                    res_graph_Omega$algo[i_res_graph_Omega] <- "GLASSO"
                    res_graph_Omega[i_res_graph_Omega,8:11] <- 
                        get_counts_Omega(abs(Omega_GL/sample_multireg$Omega)>.5,Omega!=0)
                    res_graph_Omega[i_res_graph_Omega,1:4+11] <- 
                        get_counts_Omega(abs(Omega_GL)>thr[1],Omega!=0)
                    res_graph_Omega[i_res_graph_Omega,1:4+15] <- 
                        get_counts_Omega(abs(Omega_GL)>thr[2],Omega!=0)
                    write.csv(res_graph_Omega,res_graph_Omega_file)
                    i_res_graph_Omega <- i_res_graph_Omega + 1

                    ## ad hoc
                    res_graph_Omega[i_res_graph_Omega,1:6] <- c(k,p,n,s,mod,i)
                    res_graph_Omega$algo[i_res_graph_Omega] <- "adhoc"
                    
                    res_graph_Omega[i_res_graph_Omega,1:4+11] <- 
                        get_counts_Omega(abs(Omega_adhoc)>thr[1],Omega!=0)
                    res_graph_Omega[i_res_graph_Omega,1:4+15] <- 
                        get_counts_Omega(abs(Omega_adhoc)>thr[2],Omega!=0)
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

                }
            }   
        }
    }
}
