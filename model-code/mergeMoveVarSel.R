#####
#Merge moves with variable selection - correlation (default)
#####

Rcpp::sourceCpp('MergesCpp.cpp')
mergeMoveVarSel <- function(X, model, prior, currentELBO){
  
  N = dim(X)[1]
  D = dim(X)[2]
  K = length(model$alpha)
  nCat <- as.vector(apply(X, 2, max)) #number of categories in each variable
  maxNCat <- max(nCat)
  
  eps <- model$eps
  
  
  #work out correlations, how do I save pairs? into a list? proposed_clusts is now a list
  used_clusts <- unique(model$labels)
  n_clusts <- length(used_clusts)
  all_divs <- corrBin_Calc(n_clusts, used_clusts, eps, D, maxNCat)
  all_divs[all_divs < 0.05] <- NA
  inds <- which(`dim<-`(all_divs %in% head(sort(c(all_divs), decreasing = TRUE), 3), dim(all_divs)), arr.ind = TRUE) #Take 3 smallest
  
  if(nrow(inds) == 1){
    pair <- as.vector(inds[1,])
  }else if(nrow(inds) > 1){
    pair <- as.vector(inds[sample.int(nrow(inds), 1),])
  }else{
    print("No merge clusters proposed")
    return(model)
  }
  
  k1 <- used_clusts[pair[1]]
  k2 <- used_clusts[pair[2]]
  
  mergemodel <- model
  mergemodel$rnk[,k1] <- model$rnk[,k1] + model$rnk[,k2] #add responsibilities for two clusters
  mergemodel$rnk <- mergemodel$rnk[,-k2] #remove k2 cluster
  
  mergemodel$labels <- apply(mergemodel$rnk, 1, which.max) #update current labels (keep it consistent with original E step)
  #Surely this might mean that a lot of other observations just get assigned to this new combined cluster too?
  #However I think this is how Hughes does it. If two clusters really don't explain an obs then the rnk is very close to 0 for both
  mergemodel$alpha <- model$alpha[-k2]
  mergemodel$eps <- model$eps[-k2,,] #make sure that K in the M function is now K-1
  #add c and nullphi for variable selection
  
  mergemodel <- maxStepVarSel(X, mergemodel, prior) #updates parameters for alpha, eps based on new clusters
  
  mergemodel <- expectStepVarSel(X, mergemodel) #allow obs to move in and out of new cluster
  
  mergemodel <- maxStepVarSel(X, mergemodel, prior) #update parameters for the clusters
  
  newELBO <- ELBOCalcVarSelMStep(X, mergemodel, prior)
  
  if(newELBO > currentELBO){
    mergemodel$merges <- mergemodel$merges + 1
    print("Merge accepted")
    return(mergemodel)
  }
  else{
    print("Merge rejected")
    return(model)
  } 
  
}
