//' @useDynLib RESOLDRE
//' @importFrom Rcpp sourceCpp
#include <Rcpp.h>
using namespace Rcpp;

//'Fill the matrix with ones
//'
//'@param x A matrix
//'@export
// [[Rcpp::export]]
NumericMatrix randomize(NumericMatrix x){
  int nrow = x.nrow(), ncol = x.ncol();
  NumericMatrix y(nrow, ncol);
  for (int i=0; i<nrow; i++){
    for (int j=0; j<ncol; j++){
      y(i,j)=1;
    }
  }
  return y;
}
