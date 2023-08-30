#define ARMA_NO_DEBUG

#include <armadillo>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo, BH)]]

#include <Rcpp.h>
#include <vector>

using namespace Rcpp;
using namespace arma;
using namespace RcppArmadillo;


// [[Rcpp::plugins(cpp11)]]

////////////////////////////////////////////////////////////////////////
//some functions for preparation
////////////////////////////////////////////////////////////////////////

// linear regression in Rcpp
// [[Rcpp::export]]
Rcpp::List fastLmV(const arma::mat& X, const arma::colvec& y, const int k){ 
  int n = X.n_rows; 
  
  arma::colvec coef = arma::solve(X, y);    // fit model y ~ X
  arma::colvec res  = y - X*coef;           // residuals
  
  // std.errors of coefficients
  double s2 = std::inner_product(res.begin(), res.end(), res.begin(), 0.0)/(n - k); // k is the number of covariates including intercept
  
  arma::colvec std_err = arma::sqrt(s2 * arma::diagvec(arma::pinv(arma::trans(X)*X)));
  
  arma::colvec t_value = coef/std_err; 
  return Rcpp::List::create(Named("coefficients") = coef,
                            Named("stderr")       = std_err,
                            //Named("df.residual")  = n - k,
                            Named("residuals")    = res,
                            Named("t.value")      = t_value);
}

//compute test statistic without indicators (for computing theta hat star)
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
arma::vec computeTest_C(arma::mat mydata, int n, int num_other_alpha,  int num_other_beta){
  arma::vec Y_n = mydata.col(0);
  arma::vec S_n = mydata.col(1);
  arma::vec M_n_Y = mydata.col(2);
  arma::vec M_n_S = mydata.col(3);
  
  Rcpp::List lm_beta = fastLmV( arma::conv_to<arma::mat>::from(M_n_Y), arma::conv_to<arma::vec>::from(Y_n), (num_other_beta+2) );
  Rcpp::List lm_alpha = fastLmV( arma::conv_to<arma::mat>::from(S_n),  arma::conv_to<arma::vec>::from(M_n_S), (num_other_alpha+1) );
  
  double Z_beta = as<arma::vec>(lm_beta["t.value"])(0);
  double Z_alpha = as<arma::vec>(lm_alpha["t.value"])(0);
      
  arma::vec z_stat = {Z_beta, Z_alpha};
  return z_stat;
}



//compute test statistic without indicators (for computing theta hat star)
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
arma::vec computeMed_C(arma::mat mydata, int n, int num_other_alpha,  int num_other_beta){
  arma::vec Y_n = mydata.col(0);
  arma::vec S_n = mydata.col(1);
  arma::vec M_n_Y = mydata.col(2);
  arma::vec M_n_S = mydata.col(3);
  
  Rcpp::List lm_beta = fastLmV( arma::conv_to<arma::mat>::from(M_n_Y), arma::conv_to<arma::vec>::from(Y_n), (num_other_beta+2) );
  Rcpp::List lm_alpha = fastLmV( arma::conv_to<arma::mat>::from(S_n),  arma::conv_to<arma::vec>::from(M_n_S), (num_other_alpha+1) );
  
  double Z_beta = as<arma::vec>(lm_beta["coefficients"])(0);
  double Z_alpha = as<arma::vec>(lm_alpha["coefficients"])(0);
      
  arma::vec z_stat = {Z_beta, Z_alpha};
  return z_stat;
}


