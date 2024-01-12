#' Base pglmm
#'
#' Runs phylogenetic generalized linerar model
#' @param Y
#' a column vector of 1s and 0s.
#' @param vcv
#' a variance covariance matrix for the pglmm estimation.
#' @param maxit
#' a control parameter dictating the maximum number of iterations in the optimization.
#' @param ss
#' initial estimates of the parameter for the random effect that scales the variance (it's kind of unnecessary).
#' @param tolpql
#' a control parameter dictating the tolerance for convergence in the PQL estimates of the mean components of the binomial GLMM (it's kind of unnecessary).
#' @param maxitpql
#' a control parameter dictating the maximum number of iterations in the PQL estimates of the mean components of the binomial GLMM (it's kind of unnecessary).
#' @return The probability matrix
#' @export
pglmm <- function(Y, vcv, maxit=40, ss=0.1, tolpql=10^-6, maxitpql=200){

  spp <- nrow(vcv)
  X <- rep(1,length(Y))
  Zi <- diag(spp)
  XX <- diag(spp)
  b <- rep(1,spp)

  Zi <- chol(vcv) %*% Zi
  B <- mean(Y)
  B <- log(B/(1-B))
  mu <- exp(X * B)/(1 + exp(X * B))
  Z <- X * B

  estss <- ss
  estB <- B
  oldestss <- 10^6
  oldestB <- 10^6
  it <- 0

  while ((((estss - oldestss) * (estss - oldestss) > tolpql*tolpql) |
          ((estB - oldestB) * (estB - oldestB) > tolpql*tolpql)) & (it <= maxitpql)){

          it <- it + 1

          oldestss <- estss
          oldestB <- estB

          B <- B_est(Y, X, Zi, Z, mu, B, b, XX, ss, tolpql, maxitpql)

          H <- Z - X * B

          opt <- try(optim(fn = plmm_binary_LL, par = ss, H=H, X=X, Zi=Zi, mu=mu, method = "L-BFGS-B", control = list(maxit = maxit)), TRUE)

          if (class(opt) == "try-error"){
            stop("Estimation of ss failed. You could try with a smaller s2.init, but this might not help.")
          }

          ss <- abs(opt$par)
          LL <- opt$value
          estss <- ss
          estB <- B
  }

  return(mu)

}

#' Probability estimation
#'
#' Estimates a probability matrix from an incidence matrix and a variance covariance matrix.
#' @param mat
#' an incidence matrix.
#' @param vcv
#' a variance covariance matrix for the pglmm estimation.
#' @param perspective
#' this parameter can be set to either "rows" or "columns" . It specifies the perspective taken for the estimation of the probability matrix if 'pmat' is NULL and 'cormed' is set (see details).
#' @param degree
#' this parameter determines if the randomizations preserve the degree of "columns", "rows", "both" or "sample".
#' @param maxit
#' a control parameter dictating the maximum number of iterations in the optimization.
#' @param ss
#' initial estimates of the parameter for the random effect that scales the variance (it's kind of unnecessary).
#' @param tolpql
#' a control parameter dictating the tolerance for convergence in the PQL estimates of the mean components of the binomial GLMM (it's kind of unnecessary).
#' @param maxitpql
#' a control parameter dictating the maximum number of iterations in the PQL estimates of the mean components of the binomial GLMM (it's kind of unnecessary).
#' @return The probability matrix
#' @export
probability_estimation <- function(mat, vcv, degree="both", perspective = "rows", maxit=40, ss=0.01, tolpql=10^-6, maxitpql=200){

  sstry <- unique(c(ss,0.1,0.01,0.5,10^-5,10^-10))

  if (perspective=="rows"){
    inter <- mat
  } else{
    inter <- t(mat)
  }

  prob <- c()
  for (i in 1:ncol(inter)){
    n <- sum(inter[,i])

    if(n!=0){
      for(k in 1:length(sstry)){
        result <- try(pglmm(Y = inter[,i], vcv = vcv, maxit = maxit, ss = sstry[k], tolpql = tolpql, maxitpql = maxitpql), TRUE)
        if (class(result) != "try-error"){
          if (degree!="sample"){result <- result*(1./n)}
          prob <- cbind(prob, matrix(result, nrow(inter),1))
          break
        }
      }
      if (class(result) == "try-error"){
        # warning(paste("Estimation of B and ss failed. Check for lack of variation in Y. You could try with a smaller s2.init, but this might not help."))
        prob <- cbind(prob, matrix(rep(1.0/nrow(inter), nrow(inter)), nrow(inter),1))
      }
    } else {
      result <- rep(0,nrow(inter))
      prob <- cbind(prob, matrix(result, nrow(inter),1))
    }
  }

  if (perspective=="columns"){
    prob <- t(prob)
  }

  return(prob)
}
