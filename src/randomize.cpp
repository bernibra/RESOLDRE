//' @useDynLib RESOLDRE
//' @importFrom Rcpp sourceCpp


#include <Rcpp.h>
using namespace Rcpp;

// Checking a proposed single link swap
int check_uniswap(NumericMatrix x, NumericMatrix pmat, double (*dfunc)(double,double,double,double), NumericMatrix links, int rn1, int rn2){
  double prob=dfunc( pmat(links(rn1,0),links(rn2,1)) , pmat(links(rn1,0),links(rn1,1)) , pmat(links(rn2,0),links(rn1,1)) , pmat(links(rn2,0),links(rn2,1)));
  if (rn1==rn2){
    //checking that we haven't got the same random number
    return 0;
  } else if ( unif_rand()>=prob ){
    return 0;
  } else if ( links(rn1,0)==links(rn2,0) || links(rn1,1)==links(rn2,1) ){
    //checking whether this change don't modify the interaction matrix
    return 1;
  } else if ( links(rn1,0)==links(rn2,1) || links(rn2,0)==links(rn1,1) ){
    //checking that we don't create a cannibal link
    return 0;
  } else if ( x(links(rn1,0),links(rn2,1))==1 ||
    x(links(rn2,0),links(rn1,1))==1 ||
    x(links(rn2,1),links(rn1,0))==1 ||
    x(links(rn1,1),links(rn2,0))==1 ){
    //Checking that the new elements or their permutations aren't already present in the links list
    return 0;
  } else {
    x(links(rn1,0),links(rn2,1))=1;
    x(links(rn1,0),links(rn1,1))=0;
    x(links(rn2,0),links(rn1,1))=1;
    x(links(rn2,0),links(rn2,1))=0;

    links(rn1,1)=links(rn1,1)+links(rn2,1);
    links(rn2,1)=links(rn1,1)-links(rn2,1);
    links(rn1,1)=links(rn1,1)-links(rn2,1);

    return 1;
  }
}

// Checking a proposed double link swap
int check_biswap(NumericMatrix x, NumericMatrix pmat, double (*dfunc)(double,double,double,double), NumericMatrix links, int rn1, int rn2, int rn3, int rn4){
  double prob=dfunc( pmat(links(rn1,rn3),links(rn2,1-rn4)) , pmat(links(rn1,rn3),links(rn1,1-rn3)) , pmat(links(rn2,rn4),links(rn1,1-rn3)) , pmat(links(rn2,rn4),links(rn2,1-rn4)));
  if ( rn1==rn2 ){
    //Checking that we are not selecting the same interaction
    return 0;
  } else if ( unif_rand()>=prob ){
    Rcout << "Buyaaa! " << std::endl;
    return 0;
  } else if ( links(rn1,rn3)==links(rn2,rn4) ){
    //Case in which nothing changes
    return 1;
  } else if ( x(links(rn2,rn4),links(rn1,1-rn3))==1 ||
              x(links(rn1,rn3),links(rn2,1-rn4))==1 ||
              x(links(rn1,1-rn3),links(rn2,rn4))==1 ||
              x(links(rn2,1-rn4),links(rn1,rn3))==1 ){
    //Checking that we are not creating interactions that are in conflict with others
    return 0;
  } else if ( links(rn1,rn3)==links(rn2,1-rn4) ||
              links(rn2,rn4)==links(rn1,1-rn3) ){
    //checking that we don't create a cannibal link
    return 0;
  } else {
    x(links(rn2,rn4),links(rn1,1-rn3))=1;
    x(links(rn1,rn3),links(rn2,1-rn4))=1;
    x(links(rn1,1-rn3),links(rn2,rn4))=1;
    x(links(rn2,1-rn4),links(rn1,rn3))=1;

    x(links(rn1,rn3),links(rn1,1-rn3))=0;
    x(links(rn1,1-rn3),links(rn1,rn3))=0;
    x(links(rn2,rn4),links(rn2,1-rn4))=0;
    x(links(rn2,1-rn4),links(rn2,rn4))=0;

    links(rn1,rn3)=links(rn2,rn4)+links(rn1,rn3);
    links(rn2,rn4)=links(rn1,rn3)-links(rn2,rn4);
    links(rn1,rn3)=links(rn1,rn3)-links(rn2,rn4);
    return 1;
  }
}

// Random integer
int randint(int size){
  return unif_rand()*size;
}

double probability_new(double p11, double p12, double p21, double p22){
  // Given: A <- B; C <- D
  //p11=p(A <- D); p12=p(A <- B); p21=p(C <- B); p22=p(C <- D);
  return p11*p21;
}

double probability_old(double p11, double p12, double p21, double p22){
  // Given: A <- B; C <- D
  //p11=p(A <- D); p12=p(A <- B); p21=p(C <- B); p22=p(C <- D);
  return (1-p12)*(1-p22);
}

double probability_all(double p11, double p12, double p21, double p22){
  // Given: A <- B; C <- D
  //p11=p(A <- D); p12=p(A <- B); p21=p(C <- B); p22=p(C <- D);
  return 0.5*(p11*(1-p12)+p21*(1-p22));
}


//'Randomize
//'
//'Fill the matrix with ones
//'@param x A matrix
//'@export
// [[Rcpp::export]]
List randomize(NumericMatrix x, NumericMatrix pmat, Nullable<NumericMatrix> unilinks_null, Nullable<NumericMatrix> bilinks_null, int N, int fprob){

  GetRNGstate();
  double puni, pbi;
  double (*dfunc)(double,double,double,double);
  int swaps=0;
  NumericMatrix unilinks, bilinks, mat=Rcpp::clone(x);

  if (unilinks_null.isNotNull()) {
    unilinks = Rcpp::clone(unilinks_null.get());
    puni = unilinks.nrow();
  } else {
    puni=0;
  }
  if (bilinks_null.isNotNull()) {
    bilinks = Rcpp::clone(bilinks_null.get());
    pbi = 2*bilinks.nrow();
  } else {
    pbi = 0;
  }

  if (puni+pbi == 0){
    return List::create(Named("matrix") = x, Named("swaps") = 1 );
  } else {
    puni = puni/(double)(puni+pbi);
  }

  if (fprob == 0){
    dfunc=&probability_new;
  } else if (fprob == 1){
    dfunc=&probability_old;
  } else {
    dfunc=&probability_all;
  }

  for (int j=0; j<N; j++) {
    if(unif_rand()<puni){
      swaps=swaps+check_uniswap(mat, pmat, dfunc, unilinks,randint(unilinks.nrow()),randint(unilinks.nrow()));
    } else{
      swaps=swaps+check_biswap(mat, pmat, dfunc, bilinks,randint(bilinks.nrow()),randint(bilinks.nrow()),randint(2),randint(2));
    }
  }

  PutRNGstate();
  return List::create(Named("matrix") = mat, Named("swaps") = swaps );
}

