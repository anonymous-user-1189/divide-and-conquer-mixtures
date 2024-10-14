#####
#ELBO calculation function - after M Step for efficiency
#####

ELBOCalcMStep <- function(X, model, prior){
  N = dim(X)[1]
  D = dim(X)[2]
  K = length(model$alpha)
  
  prior2 = list(alpha = rep(prior$alpha, K),
                eps = t(prior$eps))
  prioralpha <- prior2$alpha
  prioreps <- prior2$eps
  
  alpha <- model$alpha
  eps <- model$eps
  rnk <- model$rnk
  
  nCat <- as.vector(apply(X, 2, max)) #number of categories in each variable
  maxNCat <- max(nCat)
  
  Elogpi <- digamma(alpha) - digamma(sum(alpha)) #Taken from E step
  ElogphiL <- ElogphiLCalc(eps, K, D, maxNCat)
  
  Tk <- alpha - prioralpha
  
  #(log) normalising constants of Dirichlet
  Cprioralpha <- lgamma(sum(prioralpha)) - sum(lgamma(prioralpha))
  Cpostalpha <- lgamma(sum(alpha)) - sum(lgamma(alpha))
  Cprioreps <- CpriorepsCalc(prioreps, K, D, nCat)
  Cposteps <- CpostepsCalc(eps, K, D, maxNCat)
  
  #Matrix of epsilon parameters -1, where all 0's remain 0 (as these are unused parameters)
  priorepsminusone <- priorepsminusoneCalc(prioreps, K, D, maxNCat)
  epsminusone <- epsminusoneCalc(eps, K, D, maxNCat)
  epsminusprioreps <- epsminuspriorepsCalc(eps, prioreps, K, D, maxNCat)
  
  Exp1 <- sum(epsminusprioreps * ElogphiL) #E(logp(X|Z,phi))
  
  Exp2 <- sum(Tk * Elogpi) #E(logp(Z|pi))
  
  Exp3 <- sum((prioralpha - 1)*Elogpi) + Cprioralpha #E(logp(pi))
  
  Exp4 <- sum((priorepsminusone)*ElogphiL) + sum(Cprioreps) #E(logp(phi)) 
  
  logrnk <- log(rnk)
  logrnk[logrnk == -Inf] <- 0
  Exp5 <- sum(rnk * logrnk) #E(q(Z)) 
  
  Exp6 <- sum((alpha - 1)*Elogpi) + Cpostalpha #E(logq(pi))
  
  Exp7 <- sum((epsminusone)*ElogphiL) + sum(Cposteps) #Elogq(phi)
  
  ELBO <- Exp1 + Exp2 + Exp3 + Exp4 - Exp5 - Exp6 - Exp7
  
}
