######
#Global Merge - Random search
######

library(abind)
Rcpp::sourceCpp('MergesCpp.cpp')
Rcpp::sourceCpp('MixModelArmaNew.cpp')

globalMergeAlt <- function(shard_results, prior, Kinit){
  
  shard_results <- lapply(shard_results, function(x) {
    x[c("alpha", "eps", "labels", "rnk")]
  }) #Make sure these are the only values remaining. 
  
  num_clusters <- lapply(shard_results, '[[', 3 )
  num_clusters <- as.vector(unlist(lapply(num_clusters, function(x) length(unique(x))))) #clusters in each shard
  
  Ktotal = sum(num_clusters) #Total number of NON EMPTY clusters
  B = length(shard_results)
  
  prioralpha <- prior$alpha
  prioreps <- t(prior$eps)
  maxNCat <- dim(prioreps)[[1]]
  D <- dim(prioreps)[[2]]
  
  #New responsibilities with all clusters
  for (i in 1:B){
    new_rnk = matrix(0, nrow = dim(shard_results[[i]]$rnk)[1], ncol = Ktotal)
    
    if (i == 1){
      start_index = 1
      end_index = num_clusters[i]
    } else {
      start_index = sum(num_clusters[1:(i-1)]) + 1
      end_index = start_index + num_clusters[i] - 1
    }
    
    new_rnk[,start_index:end_index] <- shard_results[[i]]$rnk[,unique(shard_results[[i]]$labels)] #Only keep clusters with obs, rename these
    shard_results[[i]]$new_rnk <- new_rnk
    shard_results[[i]]$alpha <- shard_results[[i]]$alpha[unique(shard_results[[i]]$labels)]
    shard_results[[i]]$eps <- shard_results[[i]]$eps[unique(shard_results[[i]]$labels),,]
    
    
  }
  
  entropy <- initialentNoParCalc(shard_results) #CALCULATE ENTROPY FOR ALL BATCHES GIVEN COLUMNS. 
  currentELBO <- globalELBOCalcAlt(shard_results, entropy, prior, Kinit) #Should take into account entropy 
  mergerejects <- 0
  
  while(mergerejects < 10){ #Stop merging if we have rejected 10 merges in a row 
    
    num_clusters <- lapply(shard_results, '[[', 1 )
    num_clusters <- as.vector(unlist(lapply(num_clusters, function(x) length(x)))) #clusters in each shard
    
    fulleps <- lapply(shard_results, '[[', 2)
    fulleps <- fulleps[sapply(fulleps, length) > 0]
    fulleps <- lapply(fulleps, function(x) {
      if(length(dim(x)) == 2){
        
        dim(x) <- c(1, c(dim(x)))
      }
      return(x)
    })
    fulleps <- do.call("abind", c(fulleps, along=1))
    
    all_corrs <- corrBin_Calc(dim(fulleps)[1], c(1:dim(fulleps)[1]), fulleps, D, maxNCat)
    all_corrs[all_corrs < 0.05] <- NA
    if (mergerejects > 5){
      inds <- which(`dim<-`(all_corrs %in% head(sort(c(all_corrs), decreasing = TRUE)), dim(all_corrs)), arr.ind = TRUE) #Take all correlations > 0.05 now
    } else{
      inds <- which(`dim<-`(all_corrs %in% head(sort(c(all_corrs), decreasing = TRUE), 3), dim(all_corrs)), arr.ind = TRUE) #Take 3 highest
    }
    
    if(nrow(inds) == 1){
      pair <- as.vector(inds[1,])
    }else if(nrow(inds) > 1){
      pair <- as.vector(inds[sample.int(nrow(inds), 1),])
    }else{
      print("No merge clusters proposed")
      break
    }
    
    k1_pos <- pair[1]
    k2_pos <- pair[2]
    
    cluster_dictionary <- data.frame()
    cluster_dictionary <- data.frame(shard = rep(1:B, c(num_clusters)), og_clust = as.vector(unlist(sapply(num_clusters, function(x) if (x > 0) 1:x else numeric(0)))))
    
    i = cluster_dictionary$shard[k1_pos]
    k1 = cluster_dictionary$og_clust[k1_pos]
    m = cluster_dictionary$shard[k2_pos]
    k2 = cluster_dictionary$og_clust[k2_pos]
    
    k_shard2 <- num_clusters[m]
    
    if(k_shard2 == 1){
      mergemodel <- shard_results
      
      #Update responsibilities
      mergemodel <- lapply(mergemodel, function(lst) {
        lst$new_rnk[, k1_pos] <- lst$new_rnk[, k1_pos] + lst$new_rnk[,k2_pos]  # add k2 column to k1 in each matrix
        lst$new_rnk <- lst$new_rnk[,-k2_pos] #remove k2 column from each matrix
        return(lst)  # Return the modified list element
      })
      
      #Update alpha and eps for k1 
      mergemodel[[i]]$alpha[k1] <- mergemodel[[i]]$alpha[k1] + mergemodel[[m]]$alpha[k2] - prioralpha
      if(num_clusters[i] == 1){
        mergemodel[[i]]$eps <- mergemodel[[i]]$eps + mergemodel[[m]]$eps - prioreps
      } else {
        mergemodel[[i]]$eps[k1,,] <- mergemodel[[i]]$eps[k1,,] + mergemodel[[m]]$eps - prioreps 
      }
      
      mergemodel[[m]]$alpha <- mergemodel[[m]]$alpha[-k2] #remove k2 from alpha
      mergemodel[[m]]$eps <- numeric(0) #remove k2 from eps
      
      newentropy <- newentCalc(shard_results, mergemodel, entropy, k1_pos, k2_pos)
      newELBO <- globalELBOCalcAlt(mergemodel, newentropy, prior, Kinit)#FILL OUT
      
      if (newELBO > currentELBO){
        shard_results <- mergemodel
        entropy <- newentropy
        currentELBO <- newELBO
        cat("Global merge accepted, current ELBO: ", currentELBO, "\n", sep = "")
        mergerejects <- 0
      } else {
        cat("Global merge rejected, current ELBO: ", currentELBO, "\n", sep = "")
        mergerejects <- mergerejects + 1
        cat("No. of consecutive rejections: ", mergerejects, "\n", sep = "")
      }
    }
    
    if(k_shard2 > 1){
      
      mergemodel <- shard_results
      
      #Update responsibilities
      mergemodel <- lapply(mergemodel, function(lst) {
        lst$new_rnk[, k1_pos] <- lst$new_rnk[, k1_pos] + lst$new_rnk[,k2_pos]  # add k2 column to k1 in each matrix
        lst$new_rnk <- lst$new_rnk[,-k2_pos] #remove k2 column from each matrix
        return(lst)  # Return the modified list element
      })
      
      #Update alpha and eps for k1 
      mergemodel[[i]]$alpha[k1] <- mergemodel[[i]]$alpha[k1] + mergemodel[[m]]$alpha[k2] - prioralpha
      if(num_clusters[i] == 1){
        mergemodel[[i]]$eps <- mergemodel[[i]]$eps + mergemodel[[m]]$eps[k2,,] - prioreps
      } else {
        mergemodel[[i]]$eps[k1,,] <- mergemodel[[i]]$eps[k1,,] + mergemodel[[m]]$eps[k2,,] - prioreps 
      }
      
      mergemodel[[m]]$alpha <- mergemodel[[m]]$alpha[-k2] #remove k2 from alpha
      mergemodel[[m]]$eps <- mergemodel[[m]]$eps[-k2,,] #remove k2 from eps
      
      
      newentropy <- newentCalc(shard_results, mergemodel, entropy, k1_pos, k2_pos)
      newELBO <- globalELBOCalcAlt(mergemodel, newentropy, prior, Kinit)#FILL OUT
      
      if (newELBO > currentELBO){
        shard_results <- mergemodel
        entropy <- newentropy
        currentELBO <- newELBO
        cat("Global merge accepted, current ELBO: ", currentELBO, "\n", sep = "")
        mergerejects <- 0
      } else{
        cat("Global merge rejected, current ELBO: ", currentELBO, "\n", sep = "")
        mergerejects <- mergerejects + 1
        cat("No. of consecutive rejections: ", mergerejects, "\n", sep = "")
      }

    }
    
    
    
  }
  
  
  #Should have final model now
  #Find labels
  
  overalllabels <- lapply(shard_results, function(x){
    labels <- apply(x$new_rnk, 1, which.max)
    return(labels)
  })
  
  overallalpha <- unlist(lapply(shard_results, '[[', 1))
  overalleps <- lapply(shard_results, '[[', 2)
  overalleps <- overalleps[sapply(overalleps, length) > 0]
  overallrnk <- lapply(shard_results, '[[', 5)
  
  overallmodel <- list(labels = overalllabels, alpha = overallalpha, eps = overalleps, rnk = overallrnk)
  
  return(overallmodel)
}



globalELBOCalcAlt <- function(mergemodel, entropy, prior, Kinit){
  
  fullalpha <- unlist(lapply(mergemodel, '[[', 1))
  fulleps <- lapply(mergemodel, '[[', 2)
  fulleps <- fulleps[sapply(fulleps, length) > 0]
  fulleps <- lapply(fulleps, function(x) {
    if(length(dim(x)) == 2){
      
      dim(x) <- c(1, c(dim(x)))
    }
    return(x)
  })
  K <- length(fullalpha) #Total number of non-empty clusters
  
  prioralpha <- rep(prior$alpha, K) #vector
  prioreps <- t(prior$eps) #L * D matrix
  maxNCat <- dim(prioreps)[[1]]
  D <- dim(prioreps)[[2]]
  
  truealpha <- c(fullalpha, rep(prior$alpha, Kinit - K))
  trueprioralpha <- rep(prior$alpha, Kinit) #Take into account all clusters in prior
  
  Elogpi <- digamma(fullalpha) - digamma(sum(truealpha)) #Taken from E step
  ElogphiL <- lapply(fulleps, function(x) ElogphiLCalc(x, dim(x)[1], D, maxNCat)) #This calculation is independent over each K. Apply this to each eps and then sum everything
  #Don't need to worry about zombie clusters for epsilon
  
  Tk <- fullalpha - prioralpha #T_k, sufficient stat
  epsminusprioreps <- lapply(fulleps, function(x) epsminuspriorepsCalc(x, prioreps, dim(x)[1], D, maxNCat)) #s_kil, sufficient stat
  
  Cprioralpha <- lgamma(sum(trueprioralpha)) - sum(lgamma(prioralpha))
  Cpostalpha <- lgamma(sum(truealpha)) - sum(lgamma(fullalpha))
  Cprioreps <- CpriorepsCalc(prioreps, K, D, maxNCat) #makes a K x D matrix - total these for the prior constant for all non-empty clusters K 
  Cposteps <- lapply(fulleps, function(x) CpostepsCalc(x, dim(x)[1], D, maxNCat)) #makes several cluster x D matrices 
  
  priorepsminusone <- lapply(fulleps, function(x) priorepsminusoneCalc(prioreps, dim(x)[1], D, maxNCat))
  epsminusone <- lapply(fulleps, function(x) epsminusoneCalc(x, dim(x)[1], D, maxNCat))
  
  
  Exp1 <- sum(unlist(Map(function(x, y) x * y, epsminusprioreps, ElogphiL))) #E(logp(X|Z,phi)) DONE
  
  Exp2 <- sum(Tk * Elogpi) #E(logp(Z|pi)) DONE
  
  Exp3 <- sum((prioralpha - 1)*Elogpi) + Cprioralpha #E(logp(pi)), sum will sum over all k. DONE
  
  Exp4 <- sum(unlist(Map(function(x, y) x * y, priorepsminusone, ElogphiL))) + sum(Cprioreps) #E(logp(phi)) DONE 
  
  Exp5 <- entropy #DONE
  
  Exp6 <- sum((fullalpha - 1)*Elogpi) + Cpostalpha #E(logq(pi)) DONE
  
  Exp7 <- sum(unlist(Map(function(x, y) x * y, epsminusone, ElogphiL))) + sum(unlist(Cposteps)) #Elogq(phi) DONE
  
  ELBO <- Exp1 + Exp2 + Exp3 + Exp4 - Exp5 - Exp6 - Exp7
  
}

initialentNoParCalc <- function(all_shards){
  
  allents <- lapply(all_shards, '[[', 4)
  allents <- lapply(allents, function(rnk) {
    logrnk <- log(rnk)
    logrnk[logrnk == -Inf] <- 0
    return(sum(rnk * logrnk))
  })
  
  return(sum(unlist(allents)))
}


initialentParCalc <- function(all_shards){ #not sure why this is currently slower
  
  B <- length(all_shards)
  allents2 = foreach(k = 1:B) %dopar%{
    rnk <- all_shards[[k]]$rnk
    logrnk <- log(rnk)
    logrnk[logrnk == -Inf] <- 0
    return(sum(rnk * logrnk))
  }
  
  return(sum(unlist(allents2)))
  
}

newentCalc <- function(oldmodel, mergemodel, entropy, k1_pos, k2_pos){
  
  oldent <- lapply(oldmodel, '[[', 5)
  oldent <- lapply(oldent, function(rnk) {
    
    rnk <- rnk[,c(k1_pos, k2_pos)]
    logrnk <- log(rnk)
    logrnk[logrnk == -Inf] <- 0
    return(sum(rnk * logrnk))
  })
  
  newent <- lapply(mergemodel, '[[', 5)
  newent <- lapply(newent, function(rnk) {
    
    rnk <- rnk[,k1_pos]
    logrnk <- log(rnk)
    logrnk[logrnk == -Inf] <- 0
    return(sum(rnk * logrnk))
  })
  
  return( entropy - sum(unlist(oldent)) + sum(unlist(newent)))
  
}

