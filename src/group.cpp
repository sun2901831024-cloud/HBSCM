# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp ;



// [[Rcpp::export]]
arma::mat genXYg(arma::mat x,
                 arma::mat y,
                 arma::mat idBreak,
                 arma::mat G,
                 int g) {
  
  
  int k = x.n_cols ;
  int nobs = x.n_rows ;
  int T = idBreak.n_rows ;
  int k2 = y.n_cols ;
  
  arma::mat XY(k*T, k2, arma::fill::zeros) ;
  arma::mat XYg(k*g, k2, arma::fill::zeros) ;
  
  int n_start = 0 ;
  int n_end = 0 ;
  
  for (int i = 0; i < T; i++) {
    int j = G(i, 0) - 1;
    n_start = idBreak(i, 0) ;
    if (i < (T - 1)) {
      n_end = idBreak(i+1, 0) - 1 ;
    } 
    else {
      n_end = nobs - 1 ;
    }
    XY.rows(i*k, (i+1)*k - 1) = x.rows(n_start, n_end).t() * y.rows(n_start, n_end) ;
    XYg.rows(j*k, (j+1)*k - 1)+=XY.rows(i*k, (i+1)*k - 1);
  }  
  return(XYg) ;
}




// iterated alpha sampling for each unit
// [[Rcpp::export]]

arma::mat iterGenAlphag(arma::mat XX,
                       arma::mat XY,
                       arma::mat A0,
                       int N,
                       int k,
                       int r,
                       double sigma2) {
  // int k1 = k - r ;
  
  arma::mat invA0 = arma::inv(A0) ;
  
  arma::mat bara(k, 1, arma::fill::zeros) ;
  arma::mat A1(k, k, arma::fill::zeros) ;
  
  // store all random effects
  arma::mat Alpha(N, k, arma::fill::zeros) ;
  // store factor loadings
  arma::mat L(N, r, arma::fill::zeros) ;
  
  arma::mat a(k, 1, arma::fill::zeros) ;
  
  for (int i = 0; i < N; i++) {
    A1 = arma::inv(XX.rows(i*k, (i+1)*k - 1) / sigma2 + invA0) ;
    bara = A1 * XY.rows(i*k, (i+1)*k - 1) / sigma2 ;
    // we need to add restrictions on the factor loadings 
    //a = sampleSubAlpha(bara, A1, r, i+1) ; // i+1 unit indicator
    Alpha.row(i) = bara.t() ;
  }
  
  // if (r > 0) {
  //	L = Alpha.cols(k1, k - 1) ;
  // 	Alpha.shed_cols(k1, k - 1) ;
  // }
  // List result ;
  // if (k1 > 0) {
  // 	result["Alpha"] = Alpha ;
  // }
  // if (r > 0) {
  // 	result["L"] = L ;
  // }
  return(Alpha) ;
}


// iterated alpha sampling for each unit: a N*k matrix 
// [[Rcpp::export]]
arma::mat sampleAlphag(arma::mat X,
                      arma::mat y,
                      arma::mat A0,
                      arma::mat idBreak, 
                      arma::mat G,     // indicator vector
                      int N,
                      int r,
                      int g,
                      double sigma2) {
  
  int k = X.n_cols ;
  
  arma::mat XY = genXYg(X, y, idBreak,G,g) ;
  arma::mat XX = genXYg(X, X, idBreak,G,g) ;
  
  arma::mat Alpha = iterGenAlphag(XX, XY, A0, g, k, r, sigma2) ;
  
  return(Alpha) ;
}





