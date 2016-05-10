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


//'B_estim
//'
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
    if (isinf(B) | isnan(B)){
      Rcpp:stop("Estimation of B failed. Check for lack of variation in Y. You could try with a smaller s2.init, but this might not help.");
    }
  }
  Z = X * B + b + (Y - mu)/(mu % (1 - mu));
  return B;
}








//
// #include <RcppGSL.h>
// #include <gsl/gsl_min.h>
//
// // [[Rcpp::depends(RcppGSL)]]
//
// typedef struct { arma::vec H; arma::vec X; arma::mat Zt; arma::vec mu; }arguments;
//
// double plmm_binary_LL2(double par, void* p) {
//
//   arguments * params  = (arguments *)p;
//   arma::vec H = (params->H);
//   arma::vec X = (params->X);
//   arma::mat Zt = (params->Zt);
//   arma::vec mu = (params->mu);
//   par = std::abs(par);
//   arma::mat U;
//   arma::mat Ut;
//   arma::mat iV;
//   arma::mat iA;
//   double LL;
//
//   Ut = diagmat(par * X) * Zt;
//   U = trans(Ut);
//
//   iA = diagmat(mu % (1 - mu));
//   iV = Ut * iA * U;
//   iV = arma::eye<arma::mat>(U.n_rows,U.n_rows) + iV;
//
//   LL=log(std::abs(arma::det(iV))) - log(std::abs(arma::det(iA)));
//   if (isinf(LL)==1){
//     LL = 2 * accu(log(diagmat(arma::chol(iV)))) - log(std::abs(arma::det(iA)));
//   }
//
//   iV = iA - iA * U * inv(iV) * Ut * iA;
//
//   LL = 0.5 * (LL + sum(H % (iV * H)) + log(std::abs(sum(X % (iV * X)))));
//   return(LL);
// }
//
// //'PGLMM2
// //'
// //'@param x A matrix
// //'@param vcv A matrix
// //'@export
// // [[Rcpp::export]]
// List pglmm2(arma::vec Y, arma::mat vcv){
//   int status;
//   int iter = 0, max_iter = 1000;
//   double a1 = 0.0, a2 = 100.0;
//   const gsl_min_fminimizer_type *T;
//   gsl_min_fminimizer *s;
//   gsl_function F;
//
//   int spp = vcv.n_rows;
//   double B, ss = 0.5, tolpql = 10^-6, estBm, oldestBm, oldestB, oldestss, estB, estss;
//   double m = ss;
//   double exitflag, rcondflag, maxitpql = 200;
//   int it, itm;
//   arma::vec X = arma::ones<arma::vec>(spp);
//   arma::mat Zi;
//   arma::mat XX = arma::eye<arma::mat>(spp,spp);
//   arma::vec b = arma::zeros<arma::vec>(spp);
//   arma::mat C;
//   arma::vec Z;
//   arma::vec mu;
//   arma::vec H;
//   arma::mat iV;
//
//   Zi = arma::chol(vcv) * XX;
//   B = mean(Y);
//   B = log(B/(1-B));
//   mu = exp(X * B)/(1 + exp(X * B));
//   estss = ss;
//   estB = B;
//   oldestss = 10^6;
//   oldestB = 10^6;
//   it = 0;
//   exitflag = 0;
//   rcondflag = 0;
//
//   while ((((estss - oldestss) * (estss - oldestss) > tolpql*tolpql) ||
//          ((estB - oldestB) * (estB - oldestB) > tolpql*tolpql)) && (it <= maxitpql)){
//     it = it + 1;
//     oldestss = estss;
//     oldestB = estB;
//     estBm = B;
//     oldestBm = 10^6;
//     itm = 0;
//
//     while (((estBm - oldestBm) * (estBm - oldestBm) > tolpql*tolpql) && (itm <= maxitpql)){
//       itm = itm + 1;
//       oldestBm = estBm;
//       iV = plmm_binary_iV(ss, Zi, X, mu);
//       Z = X * B + b + (Y - mu)/(mu % (1 - mu));
//       B = sum(X % ( iV * Z))/(long double)sum(X % ( iV * X));
//       C = diagmat(1.0/(mu % (1 - mu)));
//       C = plmm_binary_V(ss, Zi, X, C) - C;
//       b = C * iV * (Z - X * B);
//       mu = exp(XX * b + X * B)/(1 + exp(XX * b + X * B));
//       estBm = B;
//       if (isinf(B)==1){
//         Rcpp:stop("Estimation of B failed. Check for lack of variation in Y. You could try with a smaller s2.init, but this might not help.");
//       }
//     }
//     Z = X * B + b + (Y - mu)/(mu % (1 - mu));
//     H = Z - X * B;
//
//     iter = 0;
//     a1 = 0.0;
//     a2 = 100.0;
//
//     F.function = &plmm_binary_LL2;
//
//     arguments params = {H, X, Zi, mu};
//
//
//     F.params = &params;
//
//
//     T = gsl_min_fminimizer_brent;
//     s = gsl_min_fminimizer_alloc (T);
//     gsl_min_fminimizer_set (s, &F, m, a1, a2);
//
//     do
//     {
//       iter++;
//       status = gsl_min_fminimizer_iterate (s);
//
//       m = gsl_min_fminimizer_x_minimum (s);
//       a1 = gsl_min_fminimizer_x_lower (s);
//       a2 = gsl_min_fminimizer_x_upper (s);
//
//       status = gsl_min_test_interval (a1, a2, 0.001, 0.0);
//
//     }
//     while (status == GSL_CONTINUE && iter < max_iter);
//
//     ss = m;
//     estss = ss;
//     estB = B;
//
//   }
//
//   Rcpp::checkUserInterrupt();
//   gsl_min_fminimizer_free (s);
//
//   return List::create(Named("m") = m, Named("mu") = mu, Named("ss") = estss, Named("B") = estB);
// }
