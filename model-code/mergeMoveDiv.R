#####
#Merge moves with divergences
#####

Rcpp::sourceCpp('MergesCpp.cpp')

mergeMoveDiv <- function(X, model, prior, currentELBO, divergence){
  
  N = dim(X)[1]
  D = dim(X)[2]
  K = length(model$alpha)
  nCat <- as.vector(apply(X, 2, max)) #number of categories in each variable
  maxNCat <- max(nCat)
  
  used_clusts <- unique(model$labels)
  n_clusts <- length(used_clusts)

  eps <- model$eps
  
  if(divergence == "KL"){
    all_divs <- kl_divsCalc(n_clusts, used_clusts, eps, D, maxNCat)
  } else if (divergence == "bha"){
    all_divs <- bha_divsCalc(n_clusts, used_clusts, eps, D, maxNCat)
  } else {
    print("Invalid divergence, no merge performed")
    return(model)
  }
  all_divs[lower.tri(all_divs, diag = TRUE)] <- NA
  
  inds <- which(`dim<-`(all_divs %in% head(sort(c(all_divs)), 3), dim(all_divs)), arr.ind = TRUE) #Take 3 smallest
  #Works even if there are less than 3 pairs - code will not break
  
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
  mergemodel$alpha <- model$alpha[-k2]
  mergemodel$eps <- model$eps[-k2,,] #make sure that K in the M function is now K-1
  
  mergemodel <- maxStep(X, mergemodel, prior) #updates parameters for alpha, eps based on new clusters
  
  mergemodel <- expectStep(X, mergemodel) #allow obs to move in and out of new cluster
  
  mergemodel <- maxStep(X, mergemodel, prior) #update parameters for the clusters
  
  newELBO <- ELBOCalcMStep(X, mergemodel, prior)
  
  if(newELBO > currentELBO){
    mergemodel$merges <- mergemodel$merges + 1
    print("Merge accepted")
    return(mergemodel)
  } else{
    print("Merge rejected")
    return(model)
  } 
  
}
