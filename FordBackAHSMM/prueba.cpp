#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
DataFrame prueba(NumericVector d, NumericVector st,NumericVector pdf, double pi, int jj, NumericVector N){
NumericVector D(d);
NumericVector PDF(pdf);
NumericVector ST(st);
double P(pi);
NumericVector NN(N);
int J = jj;
double FF;
double Observ;


FF = 0;
Observ = 1;
for (int u = 1; u <= 1 + 1; u++) 
    {							
      Observ *= PDF[1 - u + 1] / N[1 - u + 1];
      
      if (u < 1 + 1) 
      {
        FF += Observ * D[u] * ST[1 - u + 1];
      }
      else 
      {
       FF += Observ * D[2] * P;
      }
    }


return DataFrame::create(Named("F")= FF);
}