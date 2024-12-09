#####
#Merge moves - correlation (default)
#####

Rcpp::sourceCpp('MergesCpp.cpp')

mergeMove <- function(X, model, prior, currentELBO, Kinit){
  
  N = dim(X)[1]
  D = dim(X)[2]
  K = length(model$alpha)
  maxNCat = dim(model$eps)[[2]]
  
  eps <- model$eps
  
  used_clusts <- unique(model$labels)
  n_clusts <- length(used_clusts)
  all_divs <- corrBin_Calc(n_clusts, used_clusts, eps, D, maxNCat)
  all_divs[all_divs < 0.05] <- NA
  inds <- which(`dim<-`(all_divs %in% head(sort(c(all_divs), decreasing = TRUE), 3), dim(all_divs)), arr.ind = TRUE) #Take 3 highest
  
  if(nrow(inds) == 1){
    pair <- as.vector(inds[1,])
  }else if(nrow(inds) > 1){
    pair <- as.vector(inds[sample.int(nrow(inds), 1),])
  }else{
    cat("No merge clusters proposed", "\n", sep = "")
    return(list(model = model, ELBO = currentELBO))
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
  
  mergemodel <- expectStep(X, mergemodel, prior, Kinit) #allow obs to move in and out of new cluster
  
  mergemodel <- maxStep(X, mergemodel, prior) #update parameters for the clusters
  
  newELBO <- ELBOCalcMStep(X, mergemodel, prior, Kinit)
  
  if(newELBO > currentELBO){
    mergemodel$merges <- mergemodel$merges + 1
    cat("Merge accepted, ELBO: ", newELBO, "\n", sep = "")
    return(list(model = mergemodel, ELBO = newELBO))
  }
  else{
    cat("Merge rejected, ELBO: ", currentELBO, "\n", sep = "")
    return(list(model = model, ELBO = currentELBO))
  } 
  
}
