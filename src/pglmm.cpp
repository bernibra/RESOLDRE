//' @useDynLib RESOLDRE
//' @importFrom Rcpp sourceCpp

#include "RcppArmadillo.h"


using namespace Rcpp;


// [[Rcpp::depends(RcppArmadillo)]]


arma::mat plmm_binary_iV(double par, const arma::mat& Zi, const arma::colvec& St, const arma::vec& mu) {

  arma::mat U;
  arma::mat Ut;
  arma::mat iV;
  arma::mat iA;

  Ut = diagmat(par * St) * Zi;
  U = trans(Ut);

  iA = diagmat(mu % (1 - mu));
  iV = Ut * iA * U;

  iV = iA - iA * U * inv(arma::eye<arma::mat>(U.n_rows,U.n_rows) + iV) * Ut * iA;

  return iV;
}


arma::mat plmm_binary_V(double par, const arma::mat& Zi, const arma::colvec& St, const arma::mat& C) {

  arma::mat U;
  arma::mat Ut;

  Ut = diagmat(par * St) * Zi;
  U = trans(Ut);

  return C + U * Ut;

}


//'@export
// [[Rcpp::export]]
double plmm_binary_LL(double par, const arma::colvec& H, const arma::colvec& X, const arma::mat& Zi, const arma::colvec& mu) {

  par = std::abs(par);
  arma::mat U;
  arma::mat Ut;
  arma::mat iV;
  arma::mat iA;
  double LL;

  Ut = diagmat(par * X) * Zi;
  U = trans(Ut);

  iA = diagmat(mu % (1 - mu));
  iV = Ut * iA * U;
  iV = arma::eye<arma::mat>(U.n_rows,U.n_rows) + iV;

  LL=log(std::abs(arma::det(iV))) - log(std::abs(arma::det(iA)));
  if (isinf(LL)==1){
    LL = 2 * accu(log(diagmat(arma::chol(iV)))) - log(std::abs(arma::det(iA)));
  }

  iV = iA - iA * U * inv(iV) * Ut * iA;

  LL = 0.5 * (LL + sum(H % (iV * H)) + log(std::abs(sum(X % (iV * X)))));
  return(LL);
}


//'PGLMM
//'
//'@param x A matrix
//'@param vcv A matrix
//'@export
// [[Rcpp::export]]
long double B_est(const arma::colvec& Y, const arma::colvec& X, const arma::mat& Zi, Rcpp::NumericVector Zr, Rcpp::NumericVector mur, long double B, Rcpp::NumericVector br, const arma::mat& XX, double ss, double tolpql, int maxitpql){
  arma::colvec Z(Zr.begin(), Zr.size(), false);
  arma::colvec b(br.begin(), br.size(), false);
  arma::colvec mu(mur.begin(), mur.size(), false);
  long double estBm=B;
  double oldestBm = 10^6;
  arma::mat C;
  arma::mat iV;
  int itm = 0;

  while (((estBm - oldestBm) * (estBm - oldestBm) > tolpql*tolpql) && (itm <= maxitpql)){
    itm = itm + 1;
    oldestBm = estBm;
    iV = plmm_binary_iV(ss, Zi, X, mu);
    Z = X * B + b + (Y - mu)/(mu % (1 - mu));
    B = sum(X % ( iV * Z))/(long double)sum(X % ( iV * X));
    C = diagmat(1.0/(mu % (1 - mu)));
    C = plmm_binary_V(ss, Zi, X, C) - C;
    b = C * iV * (Z - X * B);
    mu = exp(XX * b + X * B)/(1 + exp(XX * b + X * B));
    estBm = B;
    if (isinf(B)==1){
      Rcpp:stop("Estimation of B failed. Check for lack of variation in Y. You could try with a smaller s2.init, but this might not help.");
    }
  }
  Z = X * B + b + (Y - mu)/(mu % (1 - mu));
  return B;
}

