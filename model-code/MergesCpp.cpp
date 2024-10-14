#include <RcppArmadillo.h> // new 'lighter' header
#include <boost/math/special_functions/digamma.hpp>

// Finding merge candidates - Cpp functions

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double kl_dirichlet(arma::vec alpha, arma::vec beta) {
  int n = alpha.size();
  double sum_alpha = sum(alpha);
  double term1 = lgamma(sum_alpha) - lgamma(sum(beta));
  double term2 = 0;
  double term3 = 0;
  
  for (int i = 0; i < n; i++) {
    term2 += lgamma(beta(i)) - lgamma(alpha(i));
    term3 += (alpha(i) - beta(i)) * (boost::math::digamma(alpha(i)) - boost::math::digamma(sum_alpha));
  }
  
  return term1 + term2 + term3;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat kl_divsCalc(double n_clusts, arma::vec used_clusts, arma::cube eps, double D, double maxNCat){
  arma::mat v(n_clusts, n_clusts);
  for (int i = 0; i < (n_clusts - 1); i++){
    for (int j = (i + 1); j < n_clusts; j++){
      double sum = 0; // Sum value
      for(int d = 0; d < D; d++){
        
        arma::vec clust1(maxNCat);
        arma::vec clust2(maxNCat);
        
        for (int l = 0; l < maxNCat; l++) {
          clust1(l) = eps(used_clusts(i) - 1, l, d);
          clust2(l) = eps(used_clusts(j) - 1, l, d);
        }
        
        sum += kl_dirichlet(clust1, clust2);
      }
      v(i, j) = sum;
    }
  }
  return v;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double bha_dirichlet(arma::vec alpha, arma::vec beta) {
  double term1 = 0.5 * (lgamma(arma::sum(alpha)) + lgamma(arma::sum(beta)));
  
  double term2 = 0.5 * arma::sum(lgamma(alpha) + lgamma(beta));
  
  double term3 = lgamma(0.5 * arma::sum(alpha + beta));
  
  double term4 = arma::sum(lgamma(0.5 * (alpha + beta)));
  
  return term2 - term1 + term3 - term4;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat bha_divsCalc(double n_clusts, arma::vec used_clusts, arma::cube eps, double D, double maxNCat){
  arma::mat v(n_clusts, n_clusts);
  for (int i = 0; i < (n_clusts - 1); i++){
    for (int j = (i + 1); j < n_clusts; j++){
      double sum = 0; // Sum value
      for(int d = 0; d < D; d++){
        
        arma::vec clust1(maxNCat);
        arma::vec clust2(maxNCat);
        
        for (int l = 0; l < maxNCat; l++) {
          clust1(l) = eps(used_clusts(i) - 1, l, d);
          clust2(l) = eps(used_clusts(j) - 1, l, d);
        }
        
        sum += bha_dirichlet(clust1, clust2);
      }
      v(i, j) = sum;
    }
  }
  return v;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat corrBin_Calc(double n_clusts, arma::vec used_clusts, arma::cube eps, double D, double maxNCat){
  arma::mat v(n_clusts, n_clusts);
  for (int i = 0; i < (n_clusts - 1); i++){
    for (int j = (i + 1); j < n_clusts; j++){
      arma::vec clust1(D);
      arma::vec clust2(D);
        
      for (int d = 0; d < D; d++) {
        clust1(d) = eps(used_clusts(i) - 1, 0, d);
        clust2(d) = eps(used_clusts(j) - 1, 0, d);
      }
      v(i,j) = arma::as_scalar(arma::cor(clust1, clust2));
    }
  }
  return v;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec globalcorr_Calc(double k_shard, arma::vec proposal, arma::cube eps, double D){
  arma::vec v(k_shard);
  for (int i = 0; i < k_shard; i++){
      arma::vec clust1(D);
      
      for (int d = 0; d < D; d++) {
        clust1(d) = eps(i, 0, d);
      }
      
      v(i) = arma::as_scalar(arma::cor(proposal, clust1));
  }
  return v;
}



