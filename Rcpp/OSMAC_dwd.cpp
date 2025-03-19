#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericMatrix OSMAC_dwd(NumericMatrix& x, double r, double psi, NumericVector theta) {
  int nrow = x.nrow(), ncol = x.ncol();
  NumericMatrix out(2*r,ncol+1);
  NumericVector tmp(ncol);
  
  int count = 0;
  double randomnum=0;
  for (int i = 0; i < nrow; i++) {
    NumericVector tmp1= x(i,_);
    tmp1.erase(0);
    double u = sum(tmp1*theta)*x(i,0);
    NumericVector hbar = x(i,0)*(u>0.5 ? -1.0/(4.0*u*u) : -1.0)*tmp1;
    double prob = r*sqrt(sum(hbar*hbar))/psi;
    randomnum = runif(1)[0];
    if(randomnum<prob){
      tmp = x.row(i);
      tmp.push_back(prob);
      out(count,_) = tmp;
      count++;
    }
  }
  NumericMatrix Out = out(Range(0,count-1),_);
  return Out;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//


