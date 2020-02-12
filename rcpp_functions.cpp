// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double log_sum_expS(NumericVector foo1, NumericVector foo2){
  
  double msf;
  NumericVector sumfoos = foo1 + foo2;
  double exps;
  double lse;
  int n = foo1.size();
  
  msf = max(sumfoos);
  exps = 0.0; 
  for(int i=0; i< n; i++){
    exps+=exp(sumfoos[i] -msf);
  }
  
  lse = msf + log(exps);
  
  return(lse);
}