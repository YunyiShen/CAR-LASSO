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
dt <- dt[,1:5]
head(dt)


## ----ar1example_first---------------------------------------------------------
glassores <- bGlasso(data = dt)
plot(glassores)

## ----comp_data----------------------------------------------------------------
dt <- mgp154[,c("Alistipes","Bacteroides",
                        "Eubacterium","Parabacteroides","all_others")]

## ----compositional1-----------------------------------------------------------
gut_res <- bGlasso( data = dt,link = "logit", 
                    n_iter = 2000, 
                    n_burn_in = 1000, thin_by = 2)
plot(gut_res)

## ----counting-----------------------------------------------------------------
gut_res <- gut_res <- bGlasso( data = dt[,1:4],link = "log", 
                    n_iter = 2000, 
                    n_burn_in = 1000, thin_by = 2)
plot(gut_res)

