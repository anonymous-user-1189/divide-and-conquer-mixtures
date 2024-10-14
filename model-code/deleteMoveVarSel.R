#####
#Delete moves - variable selection
######

deleteMoveVarSel <- function(X, model, prior, currentELBO){
  
  N = dim(X)[1]
  D = dim(X)[2]
  K = length(model$alpha)
  nCat <- as.vector(apply(X, 2, max)) #number of categories in each variable
  maxNCat <- max(nCat)
  
  if (length(unique(model$labels)) < 3){
    print("Delete rejected - too few clusters")
    return(model)
  }
  
  frequencies <- table(model$labels)
  proposed_clusts <- as.numeric(names(frequencies)[frequencies <= N * 0.05]) #all possible clusters to delete
  if(length(proposed_clusts) < 1){
    proposed_clusts <- as.numeric(names(sort(frequencies))[1:3]) #take three smallest clusters
  }
  
  if(length(proposed_clusts) == 1){
    k1 <- proposed_clusts
  }else{
    k1 <- sample(proposed_clusts, 1) #cluster to be deleted, make sure this is a number
  }
  
  newX <- X[model$labels != k1,] #new data frame without observations assigned to k
  
  delmodel <- model
  delmodel$alpha <- model$alpha[-k1]
  delmodel$eps <- model$eps[-k1,,]
  #add c and nullphi for variable selection
  
  delmodel <- expectStepVarSel(newX, delmodel) #updates rnk and labels based on the newX and the alpha, eps with k1 taken out
  
  delmodel <- maxStepVarSel(newX, delmodel, prior) #updates parameters for alpha, eps based on all observations without k1
  
  delmodel <- expectStepVarSel(X, delmodel) #bring back in the deleted observations with the new parameters and K-1 clusters in model, update their cluster assignments
  
  delmodel <- maxStepVarSel(X, delmodel, prior) #update parameters for the clusters
  
  newELBO <- ELBOCalcVarSelMStep(X, delmodel, prior)
  
  if(newELBO > currentELBO){
    delmodel$deletes <- delmodel$deletes + 1
    print("Delete accepted")
    return(delmodel)
  }
  else{
    print("Delete rejected")
    return(model)
  } 
  
}
