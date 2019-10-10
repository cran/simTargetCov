
[![Build Status](https://travis-ci.org/AnthonyChristidis/simTargetCov.svg?branch=master)](https://travis-ci.com/AnthonyChristidis/simTargetCov) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/simTargetCov)](https://cran.r-project.org/package=simTargetCov) [![Downloads](http://cranlogs.r-pkg.org/badges/simTargetCov)](https://cran.r-project.org/package=simTargetCov)

simTargetCov
============

This package transforms or simulates data with a target empirical covariance matrix supplied by the user.

------------------------------------------------------------------------

### Installation

You can install the **stable** version on [R CRAN](https://cran.r-project.org/package=simTargetCov).

``` r
install.packages("simTargetCov", dependencies = TRUE)
```

You can install the **development** version from [GitHub](https://github.com/AnthonyChristidis/simTargetCov).

``` r
library(devtools)
devtools::install_github("AnthonyChristidis/simTargetCov")
```

### Usage

``` r
# Function to create target covariance matrix with kernel set to r
target_cor <- function(r, p){
 Gamma <- diag(p)
 for(i in 1:(p-1)){
   for(j in (i+1):p){
     Gamma[i,j] <- Gamma[j,i] <- r^(abs(i-j))
   }
 }
 return(Gamma)
}

# Transformation of data to target empirical covariance
dat.target.cov <- simTargetCov(X = MASS::mvrnorm(30, mu = rep(0,6),
                               Sigma = target_cor(0.5,6)),
                               target = target_cor(0.5,6))
round(cov(dat.target.cov), 2)

# Simulation of data with target empirical covariance
sim.target.cov <- simTargetCov(n = 30, p = 6, target = target_cor(0.5,6))
round(cov(sim.target.cov), 2)
```

### License

This package is free and open source software, licensed under GPL (&gt;= 2).
