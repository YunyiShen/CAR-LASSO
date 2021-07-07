## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(CARlasso)

## ----ar1data------------------------------------------------------------------
set.seed(42)
dt <- simu_AR1(n=100,k=5, rho=0.7)
head(dt)


## ----ar1example_first---------------------------------------------------------
car_res <- CARlasso(y1+y2+y3+y4+y5~x1+x2+x3+x4+x5, data = dt, adaptive = TRUE)
plot(car_res,tol = 0.05)

## ----horseshoe_1--------------------------------------------------------------
# with horseshoe inference
car_res <- horseshoe(car_res)
plot(car_res)


## ----comp_data----------------------------------------------------------------
mgp154[1:5,1:7]

## ----compositional1-----------------------------------------------------------
gut_res <- CARlasso(Alistipes+Bacteroides+
                        Eubacterium+Parabacteroides+all_others~
                        BMI+Age+Gender+Stratum,
                    data = mgp154,link = "logit", 
                    adaptive = TRUE, n_iter = 2000, 
                    n_burn_in = 1000, thin_by = 2)

## ----horseshoe_comp-----------------------------------------------------------
# horseshoe will take a while, as it needs to sample the latent normal too
gut_res <- horseshoe(gut_res)
plot(gut_res)

## ----counting-----------------------------------------------------------------
gut_res <- CARlasso(Alistipes+Bacteroides+
                        Eubacterium+Parabacteroides+all_others~
                        BMI+Age+Gender+Stratum,
                    data = mgp154,link = "log", 
                    adaptive = TRUE, 
                    r_beta = 0.1, # default sometimes cause singularity in Poisson model due to exponential transformation, slightly change can fix it.
                    n_iter = 2000, 
                    n_burn_in = 1000, thin_by = 2)
# horseshoe will take a while, as it's currently implemented in R rather than C++
gut_res <- horseshoe(gut_res)
plot(gut_res)

