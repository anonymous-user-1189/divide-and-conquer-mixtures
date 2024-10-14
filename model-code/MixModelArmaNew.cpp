#include <RcppArmadillo.h> // new 'lighter' header
#include <boost/math/special_functions/digamma.hpp>

// Functions in RcppArmadillo for MerDel/sMerDel

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube ElogphiCalc(arma::cube eps, double K, double D, double N, double maxNCat, arma::mat X){
  arma::cube v(K, D, N);
  for (int k = 0; k < K; k++){
    for (int d = 0; d < D; d++){
      double sum = 0; // Sum value
      for(int i = 0; i < maxNCat; i++){
        sum += eps(k, i, d);
      }
      for (int n = 0; n < N; n++){
        double j = X(n, d);
        v(k, d, n) = boost::math::digamma(eps(k, j-1, d)) - boost::math::digamma(sum);
      }
    }
  }
  return v;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube ElogphiLCalc(arma::cube eps, double K, double D, double maxNCat){
  arma::cube v(K, maxNCat, D);
  for (int d = 0; d < D; d++){
    for (int k = 0; k < K; k++){
      double sum = 0; // Sum value
      for(int i = 0; i < maxNCat; i++){
        sum += eps(k, i, d);
      }
      for (int l = 0; l < maxNCat; l++){
        if (eps(k, l, d) != 0){
          v(k, l, d) = boost::math::digamma(eps(k, l, d)) - boost::math::digamma(sum);
        } else{
          v(k, l, d) = 0;
        }
      }
    }
  }
  return v;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat logrhonkCalc(arma::vec Elogpi, arma::cube Elogphi, double K, double D, double N){
  arma::mat v(N, K);
  for (int n = 0; n < N; n++){
    for (int k = 0; k < K; k++){
      double sum = 0; // Sum value
      for(int d = 0; d < D; d++){
        sum += Elogphi(k, d, n);
      }
      v(n, k) = Elogpi(k) + sum;
    }
  }
  return v;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube epsCalc(double K, double maxNCat, double D, double N, arma::mat prioreps, arma::mat X, arma::mat rnk){
  arma::cube v(K, maxNCat, D);
  for (int k = 0; k < K; k++){
    for (int d = 0; d < D; d++){
      for (int l = 0; l < maxNCat; l++){
        double sum = 0; // Sum value
        for(int n = 0; n < N; n++){
          if(X(n, d) == l+1){
            sum += rnk(n, k);
          } 
        }
        v(k, l, d) = prioreps(d, l) + sum;
      }
    }
  }
  return v;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube firstepsCalc(double K, double maxNCat, double D, double N, arma::mat prioreps, arma::mat X, arma::vec clusterInit){
  arma::cube v(K, maxNCat, D);
  for (int k = 0; k < K; k++){
    for (int d = 0; d < D; d++){
      for (int l = 0; l < maxNCat; l++){
        double sum = 0; // Sum value
        for(int n = 0; n < N; n++){
          if(clusterInit(n) == k + 1 && X(n,d) == l+1){
            sum += 1;
          } 
        }
        v(k, l, d) = prioreps(d, l) + sum;
      }
    }
  }
  return v;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat CpriorepsCalc(arma::mat prioreps, double K, double D, arma::vec nCat){
  arma::mat v(K, D);
  for (int d = 0; d < D; d++){
    double varCat = nCat(d);
    for (int k = 0; k < K; k++){
      double sum1 = 0; // Sum value
      for(int j = 0; j < varCat; j++){
        sum1 += prioreps(j, d);
      }
      double sum2 = 0;
      for (int j = 0; j < varCat; j++){
        sum2 += lgamma(prioreps(j, d));
      }
      v(k, d) = lgamma(sum1) - sum2;
    }
  }
  return v;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat CpostepsCalc(arma::cube eps, double K, double D, double maxNCat){
  arma::mat v(K, D);
  for (int d = 0; d < D; d++){
    for (int k = 0; k < K; k++){
      double sum1 = 0; // Sum value
      for(int j = 0; j < maxNCat; j++){
        sum1 += eps(k, j, d);
      }
      double sum2 = 0;
      for (int j = 0; j < maxNCat; j++){
        if (eps(k, j, d) != 0){
          sum2 += lgamma(eps(k, j, d));
        } else{
          sum2 += 0;
        }
      }
      v(k, d) = lgamma(sum1) - sum2;
    }
  }
  return v;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat sumDElogphiCalc(arma::cube Elogphi, double K, double D, double N){
  arma::mat v(N, K);
  for (int n = 0; n < N; n++){
    for (int k = 0; k < K; k++){
      double sum = 0; // Sum value
      for(int d = 0; d < D; d++){
        sum += Elogphi(k, d, n);
      }
      v(n, k) = sum;
    }
  }
  return v;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube priorepsminusoneCalc(arma::mat prioreps, double K, double D, double maxNCat){
  arma::cube v(K, maxNCat, D);
  for (int l = 0; l < maxNCat; l++){
    for (int k = 0; k < K; k++){
      for(int d = 0; d < D; d++){
        if (prioreps(l, d) != 0){
          v(k, l, d) = prioreps(l, d) - 1;
        } else{
          v(k, l, d) = 0;
        }
      }
    }
  }
  return v;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube epsminusoneCalc(arma::cube eps, double K, double D, double maxNCat){
  arma::cube v(K, maxNCat, D);
  for (int l = 0; l < maxNCat; l++){
    for (int k = 0; k < K; k++){
      for(int d = 0; d < D; d++){
        if (eps(k, l, d) != 0){
          v(k, l, d) = eps(k, l, d) - 1;
        } else{
          v(k, l, d) = 0;
        }
      }
    }
  }
  return v;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube epsminuspriorepsCalc(arma::cube eps, arma::mat prioreps, double K, double D, double maxNCat){
  arma::cube v(K, maxNCat, D);
  for (int l = 0; l < maxNCat; l++){
    for (int k = 0; k < K; k++){
      for(int d = 0; d < D; d++){
        if (eps(k, l, d) != 0){
          v(k, l, d) = eps(k, l, d) - prioreps(l, d);
        } else{
          v(k, l, d) = 0;
        }
      }
    }
  }
  return v;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat priorepsminusoneCalcGlob(arma::mat prioreps, double D, double maxNCat){
  arma::mat v(maxNCat, D);
  for (int l = 0; l < maxNCat; l++){
    for(int d = 0; d < D; d++){
      if (prioreps(l, d) != 0){
        v(l, d) = prioreps(l, d) - 1;
      } else{
        v(l, d) = 0;
      }
    }
  }
  return v;
}