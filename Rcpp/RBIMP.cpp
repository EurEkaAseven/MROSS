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
NumericMatrix RBIMP(NumericMatrix& x, double r, double psi, NumericVector theta) {
  int nrow = x.nrow(), ncol = x.ncol(), npos = 0, nneg = 0, tdataN=0, sumY=0;
  NumericMatrix out(3*r,ncol+1), Z2(2*r,ncol+1);
  NumericVector tmp(ncol), xsumpos(ncol), xsumneg(ncol), gradient(ncol), hbar(ncol), hbar_tmp(ncol);
  
  int count = 0;
  double randomnum=0;
  for (int i = 0; i < nrow; i++) {
    NumericVector tmp1= x(i,_);
    tmp1.erase(0);
    double pr = 1.0/(1.0+exp(-sum(tmp1*theta)));
    if(abs(x(i,0)-pr)<0.001){
      if(x(i,0)==1){
        xsumpos +=x(i,_);
        npos++;
      }else{
        xsumneg +=x(i,_);
        nneg++;
      }
    }else{
      hbar = (x(i,0)-pr)*tmp1;
      tdataN++;
      sumY +=x(i,0);
      gradient += hbar;
      double prob = r*sqrt(sum(hbar*hbar))/psi;
      randomnum = runif(1)[0];
      if(randomnum<prob){
        tmp = x.row(i);
        tmp.push_back(prob);
        out(count,_) = tmp;
        hbar_tmp = hbar;
        hbar_tmp.push_front(x(i,0));
        hbar_tmp.push_front(1.0);
        Z2(count,_) = hbar_tmp;
        count++;
      }
    }
  }
  NumericMatrix Out = out(Range(0,count-1),_);
  hbar_tmp = (xsumpos/npos);
  hbar_tmp.push_back(npos);
  out(count,_) = hbar_tmp;
  hbar_tmp = (xsumneg/nneg);
  hbar_tmp.push_back(nneg);
  out(count+1,_) = hbar_tmp;
  gradient.push_front(sumY);
  gradient.push_front(tdataN);
  out(count+2,_) = gradient;
for(int j=count+3; j< 2*count+3; j++ ){
    out(j,_) = Z2(j-count-3,_);
  }
  
  return out(Range(0,2*count+2),_);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//


