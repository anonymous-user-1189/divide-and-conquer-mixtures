#####
#Variational E Step
#####

library(matrixStats)
expectStep <- function(X, model){

  alpha <- model$alpha
  eps <- model$eps
  
  N = dim(X)[1]
  D = dim(X)[2]
  K = length(model$alpha)
  nCat <- as.vector(apply(X, 2, max)) #number of categories in each variable
  maxNCat <- max(nCat)
  
  Elogpi <- digamma(alpha) - digamma(sum(alpha)) #k dim vector, expectation of log pi k
  
  Elogphi <- ElogphiCalc(eps, K, D, N, maxNCat, X)
  
  logrhonk <- logrhonkCalc(Elogpi, Elogphi, K, D, N) #calculate rho_nk
  lse <- rowLogSumExps(logrhonk)
  rnk <- rnkCalc(logrhonk, lse, N, K)
  
  labels <- apply(rnk, 1, which.max) #k with the highest responsibility is the cluster that zn is assigned to
  
  model$rnk <- rnk #update responsibilities in model
  model$labels <- labels #update labels in model

  return(model)

}
