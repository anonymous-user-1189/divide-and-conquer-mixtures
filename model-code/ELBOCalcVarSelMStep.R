#####
#ELBO calculation function - variable selection after M Step 
#####

ELBOCalcVarSelMStep <- function(X, model, prior){
  N = dim(X)[1]
  D = dim(X)[2]
  K = length(model$alpha)
  
  prior2 = list(alpha = rep(prior$alpha, K),
                eps = t(prior$eps))
  prioralpha <- prior2$alpha
  prioreps <- prior2$eps
  a <- prior$a
  
  alpha <- model$alpha
  eps <- model$eps
  rnk <- model$rnk
  c <- model$c
  nullphi <- model$nullphi
  
  nCat <- as.vector(apply(X, 2, max)) #number of categories in each variable
  maxNCat <- max(nCat)
  
  Elogpi <- digamma(alpha) - digamma(sum(alpha)) #Taken from E step
  Elogphi <- ElogphiCalc(eps, K, D, N, maxNCat, X)
  ElogphiL <- ElogphiLCalc(eps, K, D, maxNCat)
  
  Elogdelta <- digamma(c + a) - digamma(2*a + 1)
  Elogminusdelta <- digamma(1 - c + a) - digamma(2*a + 1)
  
  Tk <- alpha - prioralpha
  
  #(log) normalising constants of Dirichlet - STAYS THE SAME COMPARED TO NO VARIABLE SELECTION
  Cprioralpha <- lgamma(sum(prioralpha)) - sum(lgamma(prioralpha))
  Cpostalpha <- lgamma(sum(alpha)) - sum(lgamma(alpha))
  Cprioreps <- CpriorepsCalc(prioreps, K, D, nCat)
  Cposteps <- CpostepsCalc(eps, K, D, maxNCat)
  Cpriordelta <- lgamma(a + a) - 2 * lgamma(a)
  Cpostdelta <- as.vector(CpostdeltaCalc(c, a, D))
  
  #Calculations
  
  carray <- replicate(N, matrix(rep(c, K), ncol = D, byrow = TRUE), simplify="array") * Elogphi
  cmatrix <- cmatrixCalc(nullphi, X, c, N, D)
  sumDElogphi <- sumDElogphiCalcVarSel(carray, cmatrix, K, D, N)
  
  #Matrix of epsilon parameters -1, where all 0's remain 0 (as these are unused parameters)
  priorepsminusone <- priorepsminusoneCalc(prioreps, K, D, maxNCat)
  epsminusone <- epsminusoneCalc(eps, K, D, maxNCat)
  
  matExp1 <- rnk * sumDElogphi
  
  Exp1 <- sum(matExp1) #E(logp(X|Z,phi)) 
  
  Exp2 <- sum(Tk * Elogpi) #E(logp(Z|pi))
  
  Exp3 <- sum((prioralpha - 1)*Elogpi) + Cprioralpha #E(logp(pi)) 
  
  Exp4 <- sum((priorepsminusone)*ElogphiL) + sum(Cprioreps) #E(logp(phi)) 
  
  Exp5 <- sum((c * Elogdelta) + (1-c)*Elogminusdelta) #E(logp(gamma|delta))
  
  Exp6 <- sum((a - 1) * Elogdelta + (a - 1) * Elogminusdelta + Cpriordelta) #E(logp(delta)) 
  
  logrnk <- log(rnk)
  logrnk[logrnk == -Inf] <- 0
  Exp7 <- sum(rnk * logrnk) #E(q(Z)) 
  
  Exp8 <- sum((alpha - 1)*Elogpi) + Cpostalpha #E(logq(pi)) 
  
  Exp9 <- sum((epsminusone)*ElogphiL) + sum(Cposteps) #Elogq(phi) 
  
  matExp10 <- (c * log(c)) + ((1-c) * log(1-c))
  matExp10[is.na(matExp10)] <- 0
  
  Exp10 <- sum(matExp10) #Elogq(gamma) 
  
  Exp11 <- sum((c + a - 1) * Elogdelta + (a - c) * Elogminusdelta + Cpostdelta) #Elogq(delta) 
  
  ELBO <- Exp1 + Exp2 + Exp3 + Exp4 + Exp5 + Exp6 - Exp7 - Exp8 - Exp9 - Exp10 - Exp11
  
}
