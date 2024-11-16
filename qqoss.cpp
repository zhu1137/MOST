
#include <RcppArmadillo.h>
#include <cmath>
#include <armadillo>
#include <queue>
#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
arma::vec lev(const arma::mat& X) {
  arma::mat U, V;
  arma::vec s;
  arma::svd_econ(U, s, V, X, "left");
  arma::vec row_crossprods(U.n_rows);
  
  for (arma::uword i = 0; i < U.n_rows; ++i) {
    arma::rowvec row = U.row(i);
    row_crossprods(i) = arma::dot(row, row);
  }
  
  arma::vec result = row_crossprods / X.n_cols;
  
  return result;
}

// [[Rcpp::export]]
arma::mat dum(arma::mat X1, arma::vec q) {
  int rows=X1.n_rows;
  int cols=X1.n_cols;
  int sum_q = 0;
  for (int i = 0; i < q.size(); i++) {
    sum_q += q[i];
  }
  arma::mat X2(rows, sum_q - cols,  arma::fill::zeros);
  
  for (int i = 0; i < rows; i++) {
    std::vector<arma::vec> e(cols);
    for (int j = 0; j < cols; j++) {
      e[j] = vector<double>(q[j] - 1, 0);
      if (X1(i,j) != q[j]) {
        e[j](X1(i,j) - 1) = 1;
      }
    }
    arma::vec row_X2;
    for (int j = 0; j < e.size(); j++) {
      row_X2.insert_rows(row_X2.n_rows, e[j]);
    }
    X2.row(i) = row_X2.t();
  }
  return X2;
}

// [[Rcpp::export]]
arma::mat DRN(arma::vec x, arma::vec p, int n) {
  int N = x.n_elem;
  arma::vec ps = arma::cumsum(p);
  arma::vec r = arma::randu(n);
  arma::mat RN(n, 1);
  for (int i = 0; i < n; i++) {
    int j = 0;
    while (r(i) > ps(j)) j++;
    RN(i) = x(j);
  }
  return RN;
}



// [[Rcpp::export]]
arma::vec bottom_k(arma::vec x, unsigned int k) {
  arma::vec x2 = x; // save a copy of x
  arma::vec ind(k); // save the indexes of the smallest k numbers
  std::nth_element(x.begin(), x.begin() + k - 1, x.end()); // std::greater<double>());
for(int ii=0, i=0; i<int(x.n_elem) && ii<int(k); i++){
  if(x2[i] <= x[k-1])  ind[ii++] = i;  // +1 for R
}
return ind;
}

// [[Rcpp::export]]
arma::vec top_k(arma::vec x, unsigned int k) {
  return bottom_k(-x,k);
}


// [[Rcpp::export]]
double delta_fun(arma::rowvec A, arma::rowvec b, arma::rowvec q) {
  int n = arma::size(b, 1);
  arma::rowvec D(n, fill::zeros); 
  for (int j = 0; j < n; j++) {
    if (A(j) == b(j)) {
      D(j) = q(j);
    }
  }
  double result = arma::sum(D);
  return result;
}

// [[Rcpp::export]]
arma::vec Dscr_oss(arma::mat X, arma::vec xa, arma::mat y, double ya, double tPow) {
  // set tPow=2 for D2 and tPow=4 for D4
  // y is a row vector
  int n=X.n_rows;
  int p=X.n_cols;
  arma::vec B = zeros<vec>(n);
  for(int i=0; i<n; i++){
    B(i) = pow(accu(X.row(i)==y)+p-xa(i)/2-ya/2,tPow); // current used
    //B(i) = pow(p-accu(abs(X.row(i)-y))/2,tPow); // not work
  }
  return B;
}

// [[Rcpp::export]]
arma::uvec OSS_cpp(arma::mat x, int k, double tPow=2){
  int n=x.n_rows;
  arma::uvec candi=linspace<uvec>(1,n,n);
  arma::uvec ind=linspace<uvec>(1,k,k);
  arma::vec L=sum(pow(x,2),1);
  arma::vec xa=L;
  uword mm=L.index_max();
  ind(0)=candi(mm);
  candi.shed_row(mm);
  L.shed_row(mm);
  
  arma::mat sx=sign(x);
  double r=log(n/k)/log(k);
  for(int i=1; i<k; i++){
    if(i==1)
      L=Dscr_oss(sx.rows(candi-1),xa.elem(candi-1),sx.row(ind(i-1)-1),xa(ind(i-1)-1),tPow);
    else
      L=L+Dscr_oss(sx.rows(candi-1),xa.elem(candi-1),sx.row(ind(i-1)-1),xa(ind(i-1)-1),tPow);
    
    
    mm=L.index_min();
    ind(i)=candi(mm);
    candi.shed_row(mm);
    L.shed_row(mm);
    
    int nc=floor(n/pow(i,r));
    //Rcout << ind(i) << std::endl;
    //double nc=n/pow(i,r)/L.n_elem;
    if((i>1) & (L.n_elem>double(nc))){
      //arma::uvec tt=as<arma::uvec>(bottom_k(L,nc));
      //Rcout << bottom_k(L,nc) << std::endl;
      arma::uvec tt=arma::conv_to<arma::uvec>::from(bottom_k(L,nc));
      L=L.elem(tt);
      candi=candi.elem(tt);
    }
  }
  return ind;
}

// [[Rcpp::export]]
arma::vec Dscr_woss(arma::mat X, arma::vec xa, arma::mat y, double ya, arma::rowvec qv, double tPow) {
  int n=X.n_rows;
  int p=X.n_cols;
  arma::vec B = zeros<vec>(n);
  for(int i=0; i<n; i++){
    B(i) = pow(delta_fun(X.row(i), y, qv) + (p -xa(i)/2-ya/2),tPow); // current used
  }
  return B;
}


// [[Rcpp::export]]
arma::uvec WOSS_cpp(arma::mat x, int k, arma::rowvec qv, double tPow=2){
  int n=x.n_rows;
  arma::uvec candi=linspace<uvec>(1,n,n);
  arma::uvec ind=linspace<uvec>(1,k,k);
  // arma::rowvec qnv=qv.subvec(x.n_cols, x.n_cols+z.n_cols-1);
  // arma::vec L=sum(pow(z*diagmat(sqrt(qnv)),2),1);
  arma::vec L=sum(pow(x,2),1);
  arma::vec xa=L;
  uword mm=L.index_max();
  ind(0)=candi(mm);
  candi.shed_row(mm);
  L.shed_row(mm);
  
  arma::mat sx=sign(x);
  double r=log(n/k)/log(k);
  for(int i=1; i<k; i++){
    if(i==1)
      L=Dscr_woss(sx.rows(candi-1),xa.elem(candi-1),sx.row(ind(i-1)-1),xa(ind(i-1)-1), qv, tPow);
    else
      L=L+Dscr_woss(sx.rows(candi-1),xa.elem(candi-1),sx.row(ind(i-1)-1),xa(ind(i-1)-1), qv, tPow);
    
    
    mm=L.index_min();
    ind(i)=candi(mm);
    candi.shed_row(mm);
    L.shed_row(mm);
    
    int nc=floor(n/pow(i,r));
    if((i>1) & (L.n_elem>double(nc))){
      arma::uvec tt=arma::conv_to<arma::uvec>::from(bottom_k(L,nc));
      L=L.elem(tt);
      candi=candi.elem(tt);
    }
  }
  return ind;
}


// [[Rcpp::export]]
arma::vec Dscr_qqs(arma::mat OX, arma::mat OZS, arma::vec OZE, arma::rowvec b, double bE, arma::rowvec qv, double h, int tPow) {
  int n = OX.n_rows;
  int p = OZS.n_cols;
  arma::mat W = join_rows(OX, OZS);
  arma::vec B = zeros<vec>(n);
  for (int i = 0; i < n; i++) {
    B(i) = pow(h * (OZS.n_cols - OZE(i)/2 - bE/2) + delta_fun(W.row(i), b, qv), tPow);
  }
  return B;
}


// [[Rcpp::export]]
arma::uvec OAJ2_qqs(arma::mat x, arma::mat z, int k, arma::rowvec qv, double h, int tPow=2){
  int n=x.n_rows;
  arma::uvec candi=linspace<uvec>(1,n,n);
  arma::uvec ind=linspace<uvec>(1,k,k);
  arma::vec L=sum(pow(z,2),1);
  arma::vec za=L;
  uword mm=L.index_max();
  ind(0)=candi(mm);
  candi.shed_row(mm);
  L.shed_row(mm);
  
  arma::mat sz = sign(z);
  arma::mat sw = join_rows(x, sz);
  double r=log(n/k)/log(k);
  for(int i=1; i<k; i++){
    if(i==1)
      L = Dscr_qqs(x.rows(candi-1), sz.rows(candi-1), za.elem(candi-1), sw.row(ind(i-1)-1), za(ind(i-1)-1), qv, h, tPow);
    else
      L = L + Dscr_qqs(x.rows(candi-1), sz.rows(candi-1), za.elem(candi-1), sw.row(ind(i-1)-1), za(ind(i-1)-1), qv, h, tPow);
    
    mm=L.index_min();
    ind(i)=candi(mm);
    candi.shed_row(mm);
    L.shed_row(mm);
    
    int nc=floor(n/pow(i,r));
    if((i>1) & (L.n_elem>double(nc))){
      arma::uvec tt=arma::conv_to<arma::uvec>::from(bottom_k(L,nc));
      L=L.elem(tt);
      candi=candi.elem(tt);
    }
  }
  return ind;
}



// [[Rcpp::export]]
arma::vec Dscr_ds(arma::mat OX, arma::mat OZS, arma::vec OZE, arma::rowvec ba, arma::rowvec bb, double bE, arma::rowvec qv, double h, int tPow) {
  int n = OX.n_rows;
  int p = OZS.n_cols;
  arma::mat W = join_rows(OX, OZS);
  arma::vec B = zeros<vec>(n);
  for (int i = 0; i < n; i++) {
    // B(i) = pow(h * (OZS.n_cols - OZE(i)/2 - bE/2) + delta_fun(W.row(i), b, qv), tPow);
    B(i) = pow(delta_fun(OX.row(i), ba, qv), tPow) + pow(accu(OZS.row(i)==bb) + OZS.n_cols - OZE(i)/2 - bE/2,tPow);
  }
  return B;
}


// [[Rcpp::export]]
arma::uvec OAJ2_ds(arma::mat x, arma::mat z, int k, arma::rowvec qv, double h, int tPow=2){
  int n=x.n_rows;
  arma::uvec candi=linspace<uvec>(1,n,n);
  arma::uvec ind=linspace<uvec>(1,k,k);
  arma::vec L=sum(pow(z,2),1);
  arma::vec za=L;
  uword mm=L.index_max();
  ind(0)=candi(mm);
  candi.shed_row(mm);
  L.shed_row(mm);
  
  arma::mat sz = sign(z);
  //arma::mat sw = join_rows(x, sz);
  double r=log(n/k)/log(k);
  for(int i=1; i<k; i++){
    if(i==1)
      L = Dscr_ds(x.rows(candi-1), sz.rows(candi-1), za.elem(candi-1), x.row(ind(i-1)-1), sz.row(ind(i-1)-1), za(ind(i-1)-1), qv, h, tPow);
    else
      L = L + Dscr_ds(x.rows(candi-1), sz.rows(candi-1), za.elem(candi-1), x.row(ind(i-1)-1), sz.row(ind(i-1)-1), za(ind(i-1)-1), qv, h, tPow);
    
    mm=L.index_min();
    ind(i)=candi(mm);
    candi.shed_row(mm);
    L.shed_row(mm);
    
    int nc=floor(n/pow(i,r));
    if((i>1) & (L.n_elem>double(nc))){
      arma::uvec tt=arma::conv_to<arma::uvec>::from(bottom_k(L,nc));
      L=L.elem(tt);
      candi=candi.elem(tt);
    }
  }
  return ind;
}





