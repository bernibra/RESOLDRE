#' Base pglmm
#'
#' Probability estimation
#' @return The probability matrix
#' @export
pglmm <- function(Y, vcv){
  iter <- 0
  max_iter <- 40
  maxit <- 40
  spp <- nrow(vcv)
  ss <- 0.002
  tolpql <- 10^-6
  maxitpql <- 200;
  X <- rep(1,length(Y))
  Zi <- diag(spp)
  XX <- diag(spp)
  b <- rep(1,spp)
  Zi <- chol(vcv) %*% Zi
  B <- mean(Y)
  B <- log(B/(1-B))
  mu <- exp(X * B)/(1 + exp(X * B))
  estss <- ss
  estB <- B
  oldestss <- 10^6
  oldestB <- 10^6
  it <- 0
  exitflag <- 0
  rcondflag <- 0
  Z <- X * B
  nested <- list(NULL)

  while ((((estss - oldestss) * (estss - oldestss) > tolpql*tolpql) |
          ((estB - oldestB) * (estB - oldestB) > tolpql*tolpql)) & (it <= maxitpql)){

          it <- it + 1

          oldestss <- estss
          oldestB <- estB

          B <- B_est(Y, X, Zi, Z, mu, B, b, XX, ss, tolpql, maxitpql)

          H <- Z - X * B

          opt <- optim(fn = plmm_binary_LL, par = ss, H=H, X=X, Zi=Zi, mu=mu, method = "L-BFGS-B",lower = 0, upper = 100, control = list(maxit = maxit))

          ss <- abs(opt$par)
          LL <- opt$value
          estss <- ss
          estB <- B
  }

  return(mu)

}
