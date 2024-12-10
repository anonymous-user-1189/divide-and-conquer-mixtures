maxStepPar <- function(X, model, prior){
  
  prioralpha <- prior$alpha
  prioreps <- prior$eps
  rnk <- model$rnk
  
  N = dim(X)[1]
  D = dim(X)[2]
  K = length(model$alpha)
  nCat <- as.vector(apply(X, 2, max)) #number of categories in each variable
  maxNCat <- max(nCat)
  
  parCalc <- foreach(k = 1:4) %dorng% {
    start_idx <- (k - 1) * 5000 + 1
    end_idx <- k * 5000
    
    rnkchunk <- rnk[start_idx:end_idx,]
    
    alphaSums <- colSums(rnkchunk)
    
    epsSums <- epsCalc(K, maxNCat, D, 5000, prioreps, X[start_idx:end_idx,], rnkchunk) - prior$eps[1,1]
    
    list(alphaSums, epsSums)
  }
  
  
  #Parameters for pi update - Dirichlet
  alpha <- prioralpha + Reduce(`+`, lapply(parCalc, '[[', 1))
  
  #Parameters for phi update - Dirichlet
  eps <- prior$eps[1,1] + Reduce(`+`, lapply(parCalc, '[[', 2))
  
  model$alpha <- alpha #update alpha* in model
  model$eps <- eps #update epsilon* in model
  return(model)
}


#library(microbenchmark)
#benchmark_results <- microbenchmark(
#  Method1 = maxStepPar(X, model, prior),             # Code block 1
#  Method2 = maxStep(X, model, prior),     # Code block 2
#  times = 10  # Number of iterations
#)

