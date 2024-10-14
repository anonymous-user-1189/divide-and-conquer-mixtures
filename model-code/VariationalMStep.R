#####
#Variational M Step
#####

maxStep <- function(X, model, prior){

  prioralpha <- prior$alpha
  prioreps <- prior$eps
  rnk <- model$rnk
  
  N = dim(X)[1]
  D = dim(X)[2]
  K = length(model$alpha)
  nCat <- as.vector(apply(X, 2, max)) #number of categories in each variable
  maxNCat <- max(nCat)
  
  #Parameters for pi update - Dirichlet
  alpha <- prioralpha + colSums(rnk)
  
  #Parameters for phi update - Dirichlet
  eps <- epsCalc(K, maxNCat, D, N, prioreps, X, rnk)
  
  model$alpha <- alpha #update alpha* in model
  model$eps <- eps #update epsilon* in model
  return(model)
}
