library(Matrix)
library(RcppArmadillo)
library(Rcpp)
library(RcppProgress)

rm(list = ls())


Rcodes <- list.files("./R",".R$", full.names = T)
lapply(Rcodes,source)
sourceCpp("./src/CAR-ALASSO.cpp")
sourceCpp("./src/CAR-LASSO.cpp")
sourceCpp("./src/graphical-intercept-LASSO.cpp")
sourceCpp("./src/CAR-ALASSO-hir.cpp")
sourceCpp("./src/CAR-LASSO-hir.cpp")
sourceCpp("./src/multireg_wrap.cpp")

set.seed(42)
dt <- simu_AR1()
car_res <- CARlasso(y1+y2+y3+y4+y5~x1+x2+x3+x4+x5, data = dt, adaptive = TRUE)
plot(car_res,tol = 0.05)
#' # with horseshoe inference
car_res <- horseshoe(car_res)
plot(car_res)

glassores <- bGlasso(data = dt[,1:5])
plot(glassores, tol = 0.05) 
