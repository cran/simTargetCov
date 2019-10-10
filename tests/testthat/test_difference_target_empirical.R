# --------------------------------------------------------
# Test Script - Error for Resulting Empirical Correlation
# --------------------------------------------------------

# Required libraries
library(simTargetCov)

# Context of test script
context("Verify matching between target and resulting empirical covariance matrix.")

# There should be an error if we want to compute the IF TS, and no returns are provided
test_that("Difference between target and resulting empirical correlation should be near 0.", {

  # Target correlation
  target_cor <- function(r, p){
    Gamma <- diag(p)
    for(i in 1:(p-1)){
      for(j in (i+1):p){
        Gamma[i,j] <- Gamma[j,i] <- r^(abs(i-j))
      }
    }
    return(Gamma)
  }

  # Case where data is provided by user
  dat.target.cov <- simTargetCov(X = MASS::mvrnorm(30, mu = rep(0,6), Sigma=target_cor(0.5,6)), target = target_cor(0.5,6))
  expect_false(any(abs(cov(dat.target.cov) - target_cor(0.5,6)) > 1e-10), FALSE)

  # Case where data is simulated internally by the function
  sim.target.cov <- simTargetCov(n = 30, p = 6, target = target_cor(0.5,6))
  expect_false(any(abs(cov(sim.target.cov) - target_cor(0.5,6)) > 1e-10), FALSE)

  }
)

