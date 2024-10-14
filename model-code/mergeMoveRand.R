#####
#Merge moves with random clusters selected for merge
#####

mergeMoveRand <- function(X, model, prior, currentELBO){
  
  N = dim(X)[1]
  D = dim(X)[2]
  K = length(model$alpha)
  nCat <- as.vector(apply(X, 2, max)) #number of categories in each variable
  maxNCat <- max(nCat)
  
  
  used_clusts <- unique(model$labels)
  n_clusts <- length(used_clusts)
  
  if(n_clusts > 1){
    pair <- sample(used_clusts, 2)
  }else{
    print("No merge clusters proposed")
    return(model)
  }
  
  k1 <- pair[1]
  k2 <- pair[2]
  
  mergemodel <- model
  mergemodel$rnk[,k1] <- model$rnk[,k1] + model$rnk[,k2] #add responsibilities for two clusters
  mergemodel$rnk <- mergemodel$rnk[,-k2] #remove k2 cluster
  
  mergemodel$labels <- apply(mergemodel$rnk, 1, which.max) #update current labels (keep it consistent with original E step)
  #Surely this might mean that a lot of other observations just get assigned to this new combined cluster too?
  #However I think this is how Hughes does it. If two clusters really don't explain an obs then the rnk is very close to 0 for both
  mergemodel$alpha <- model$alpha[-k2]
  mergemodel$eps <- model$eps[-k2,,] #make sure that K in the M function is now K-1
  #add c and nullphi for variable selection
  
  mergemodel <- maxStep(X, mergemodel, prior) #updates parameters for alpha, eps based on new clusters
  
  mergemodel <- expectStep(X, mergemodel) #allow obs to move in and out of new cluster
  
  mergemodel <- maxStep(X, mergemodel, prior) #update parameters for the clusters
  
  newELBO <- ELBOCalcMStep(X, mergemodel, prior)
  
  if(newELBO > currentELBO){
    mergemodel$merges <- mergemodel$merges + 1
    print("Merge accepted")
    return(mergemodel)
  }
  else{
    print("Merge rejected")
    return(model)
  } #how do we know if this has been accepted and how many moves have happened?
  
}