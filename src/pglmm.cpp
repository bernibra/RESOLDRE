//' @useDynLib RESOLDRE
//' @importFrom Rcpp sourceCpp
#include <Rcpp.h>
using namespace Rcpp;

//'PGLMM
//'
//'@param x A matrix
//'@param vcv A matrix
//'@export
// [[Rcpp::export]]
NumericMatrix pglmm(NumericMatrix x, NumericMatrix vcv){
  int nrow = x.nrow(), ncol = x.ncol();
  NumericMatrix y(nrow, ncol);
  for (int i=0; i<nrow; i++){
    for (int j=0; j<ncol; j++){
      y(i,j)=1;
    }
  }
  return y;
}
