#' @import stats
#'
#' @title Data Transformation or Simulation with Target Empirical Covariance Matrix
#'
#' @description \code{simTargetCov} transforms or simulates data with a target empirical covariance matrix supplied by the user.
#'
#' @param n Number of observations for data matrix output.
#' @param p Number of variables for data matrix output.
#' @param target Target empirical covariance for data matrix output.
#' @param X Data matrix for transformation.
#'
#' @export
#'
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#'
#' @examples
#' # Function to create target covariance matrix with kernel set to r
#' target_cor <- function(r, p){
#'   Gamma <- diag(p)
#'   for(i in 1:(p-1)){
#'     for(j in (i+1):p){
#'       Gamma[i,j] <- Gamma[j,i] <- r^(abs(i-j))
#'     }
#'   }
#'   return(Gamma)
#' }
#'
#' # Transformation of data to target empirical covariance
#' dat.target.cov <- simTargetCov(X = MASS::mvrnorm(30, mu = rep(0,6),
#'                                Sigma = target_cor(0.5,6)),
#'                                target = target_cor(0.5,6))
#' round(cov(dat.target.cov), 2)
#'
#' # Simulation of data with target empirical covariance
#' sim.target.cov <- simTargetCov(n = 30, p = 6, target = target_cor(0.5,6))
#' round(cov(sim.target.cov), 2)
#'
simTargetCov <- function(n, p, target, X=NULL){

  # Input check for target
  if(!is.matrix(target))
    stop("The argument target must be a matrix.") else{
      if(nrow(target)!=ncol(target))
        stop("The argument target must be a square matrix.")
      if(ncol(target) <=1)
        stop("The number of rows and columns of target must be greater than one.")
    }

  # Input check for X, n and p
  if(!is.null(X)){
    if(!is.matrix(X) && !is.data.frame(X))
      stop("X must be a matrix or a data frame.")
    if(any(ncol(X)<=1, nrow(X)<=1))
      stop("The number of rows and columns of X must be greater than 1.")
    if(any(nrow(target)!=ncol(X), ncol(target)!=ncol(X)))
      stop("The dimension of the target matrix does not match the number of variables in the matrix X.")
  } else{
    if(any(nrow(X)<=1, ncol(X)<=1))
      stop("The arguments n and p must be greater than 1.")
    if(any(nrow(target)!=p, ncol(target)!=p))
      stop("The dimension of the target matrix does not match the number of variables p.")
  }

  # Case where the original data is provided by the user
  if(is.null(X))
    X <- MASS::mvrnorm(n, mu=rep(0,p), Sigma=target)

  # Storing the column means of the data
  X.colMeans <- colMeans(X)

  # Covariance matrix of scaled data
  X.tilde <- scale(X)
  S <- cov(X.tilde)
  # Spectral decomposition
  res <- eigen(S, symmetric=TRUE)
  P <- res$vectors;
  # Random orthogonal columns matrix
  Z <- X.tilde%*%P
  Z.tilde <- scale(Z)
  # Transformation by rotation of data matrix
  res.target.cor <- eigen(target, symmetric=TRUE)
  Y <- t(res.target.cor$vectors%*%sqrt(diag(res.target.cor$values))%*%t(Z.tilde))

  # Adding the column means back to original data
  Y <- sapply(1:ncol(Y), function(k, Y, col.mean){Y[,k]+col.mean[k]}, Y=Y, col.mean=X.colMeans)

  # Return data with target empirical covariance matrix
  return(Y)
}










