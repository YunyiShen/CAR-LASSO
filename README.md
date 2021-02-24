![](https://github.com/YunyiShen/CAR-LASSO/workflows/R-CMD-check/badge.svg)
# CAR-LASSO
Conditional Auto-Regressive LASSO in R

R implementation of 

Shen, Yunyi, and Claudia Solis-Lemus. "Bayesian Conditional Auto-Regressive LASSO Models to Learn Sparse Networks with Predictors." arXiv preprint arXiv:2012.08397 (2020).

The package is not yet on CRAN, to install it, use:
```r
devtools::install_github("YunyiShen/CAR-LASSO")
```

Example run, we generated data from a 5-node AR1 model with each node has a specific treatment then use CAR-ALASSO to reconstruct such network and plot the result. 

```r
set.seed(42)
dt <- simu_AR1(n=100,k=5, rho=0.7)
car_res <- CARlasso(y1+y2+y3+y4+y5~x1+x2+x3+x4+x5, data = dt, adaptive = TRUE)
plot(car_res,tol = 0.05)
# with horseshoe inference
car_res <- horseshoe(car_res)
plot(car_res)
```
