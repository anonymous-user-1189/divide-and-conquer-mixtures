######
#Implement MerDel with clusters selected for merging via divergence, time measured for each 'laps' model. *NO laps = 0 SIMULATION*
#'mixturemodel_merdeldiv' (make sure this isn't confused with simulation including laps = 0)
#####

library(klaR)
library(tidyverse)
library(matrixStats)
source('VariationalEStep.R')
source('VariationalMStep.R')
source('deleteMove.R')
source('mergeMoveDiv.R')
source('ELBOCalc.R')
source('ELBOCalcMStep.R')

library(Rcpp)
library(RcppArmadillo)
Rcpp::sourceCpp('MixModel.cpp')
Rcpp::sourceCpp('MixModelArmaNew.cpp')

#Variational mixtures for categorical distributions
check_convergence<- function(ELBO, iter, maxiter, tol){
  if (iter > 1 && abs(ELBO[iter * 2] - ELBO[iter * 2-1]) < tol && abs(ELBO[iter * 2-1] - ELBO[iter * 2-2] < tol) && 
      abs(ELBO[iter * 2-2] - ELBO[iter * 2-3] < tol)){
    print(paste("Stopped after iteration ",iter)) #make sure the last 3 ELBOs close to each other
    return(TRUE)
  }
  if (iter == maxiter){
    print(paste("Not converged after maximum number of iterations"))
    return(TRUE)
  }
  else{
    return(FALSE)
  }
}

#data = 'X' observed variables
#K = number of clusters
#alpha = prior parameters for clusters
#maxiter = maximum number of iterations
#tol = convergence criteria
#laps = number of E-M iterations to do before performing merge/delete. Numeric or "all" ie. no merge/delete moves

mixturemodel_merdeldiv <- function(data, K, alpha, maxiter, tol, laps, divergence){
  ################################################
  #Process dataset
  start.time <- Sys.time()
  #First make sure X (data) is a data frame and make every column into factors
  X <- as.data.frame(data)
  X[colnames(X)] <- lapply(X[colnames(X)], as.factor)
  
  #Data frame mapping factor labels to numbers
  factor_labels = data.frame()
  for (i in 1:length(colnames(X))){
    factorlist <- data.frame(factors = unique(X[,i]), value = as.numeric(unique(X[,i])))
    factordict <- cbind(data.frame(variable = colnames(X)[i]), factorlist)
    factor_labels <- rbind(factor_labels, factordict)
  }
  
  #Check for any completely empty factors as these may cause problems
  categories <- lapply(1:ncol(X), function(j)levels(X[, j])) 
  cat_lengths <- sapply(categories, length)
  if(any(cat_lengths == 0)) {
    stop("Column(s) ", paste0(colnames(X)[cat_lengths == 0], collapse = ","), " is empty")
  }
  
  X <- data.matrix(X)
  
  N = dim(X)[1]
  D = dim(X)[2]
  
  nCat <- as.vector(apply(X, 2, max)) #number of categories in each variable
  maxNCat <- max(nCat)
  
  ########################################
  #Define priors - given input alpha, all variables have a symmetric Dirichlet prior with parameter alpha
  
  prior = list(alpha = alpha) #prior for clustering pi
  prior$eps = matrix(0, D, maxNCat) #prior for clustering phi
  for(d in 1:D){
    prior$eps [d,1:nCat[d]] <- 1/nCat[d]
  }
  
  clusterInit <- klaR::kmodes(X, modes = K)$cluster #k modes: analogous to k means
  end.time <- Sys.time()
  total_time1 <- difftime(end.time, start.time, unit = "secs")
  all_models <- list()
  for (k in 1:length(laps)){
    start.time2 <- Sys.time()
    ELBO = Cl = secondELBO = rep(-Inf,maxiter*2) #record ELBO and no of clusters at each iteration
    secondELBO <- as.data.frame(cbind(secondELBO, Time = NA))
    
    model = list(alpha = rep(prior$alpha, K),
                 eps = firstepsCalc(K, maxNCat, D, N, prior$eps, X, clusterInit),
                 labels = clusterInit,
                 totalmoves = 0,
                 merges = 0,
                 deletes = 0) 
    
    
    ########################################
    ELBOcounter <- 1
    
    lap <- laps[k]
    for (iter in 1:maxiter){
      if (lap == "all"){
        print(paste("Iteration number ",iter))
        model = expectStep(X, model) #Expectation step
        .GlobalEnv$maxNCat <- maxNCat
        ELBO[iter * 2-1] = secondELBO[iter * 2-1,1] = ELBOCalc(X, model, prior) #ELBO
        secondELBO[iter * 2-1,2] <- Sys.time()
        print(ELBO[iter * 2-1])
        
        Cl[iter] = length(unique(model$labels)) #Counts number of non-empty clusters
        
        model = maxStep(X, model, prior) #Maximisation step
        ELBO[iter * 2] = secondELBO[iter * 2,1] = ELBOCalcMStep(X, model, prior) #ELBO
        secondELBO[iter * 2,2] <- Sys.time()
        print(ELBO[iter * 2])
        
        Cl[iter] = length(unique(model$labels)) #Counts number of non-empty clusters
        
        if(check_convergence(ELBO, iter, maxiter, tol)) break
        
      } else if (lap %% 1 == 0){
        
        #E step
        print(paste("Iteration number ",iter))
        model = expectStep(X, model) #Expectation step
        .GlobalEnv$maxNCat <- maxNCat
        ELBO[iter * 2 - 1] = secondELBO[ELBOcounter,1] = ELBOCalc(X, model, prior) #ELBO
        secondELBO[ELBOcounter,2] <- Sys.time()
        print(secondELBO[ELBOcounter,1])
        ELBOcounter <- ELBOcounter + 1
        
        Cl[iter] = length(unique(model$labels)) #Counts number of non-empty clusters
        
        #M step
        model = maxStep(X, model, prior) #Maximisation step
        ELBO[iter * 2] = secondELBO[ELBOcounter,1] = ELBOCalcMStep(X, model, prior) #ELBO
        secondELBO[ELBOcounter,2] <- Sys.time()
        print(secondELBO[ELBOcounter,1])
        ELBOcounter <- ELBOcounter + 1
        
        Cl[iter] = length(unique(model$labels)) #Counts number of non-empty clusters
        
        if(check_convergence(ELBO, iter, maxiter, tol)) break
        
        #Merge delete: only happens every "laps" number of EM iterations
        if (iter %% lap == 0){
          model = mergeMoveDiv(X, model, prior, secondELBO[ELBOcounter - 1,1], divergence)
          secondELBO[ELBOcounter:(ELBOcounter+3),1] = rep(ELBOCalcMStep(X, model, prior), 4)
          secondELBO[ELBOcounter:(ELBOcounter+3),2] <- Sys.time()
          ELBOcounter <- ELBOcounter + 4
          model = deleteMove(X, model, prior, secondELBO[ELBOcounter - 1,1])
          secondELBO[ELBOcounter:(ELBOcounter+3),1] = rep(ELBOCalcMStep(X, model, prior), 4)
          secondELBO[ELBOcounter:(ELBOcounter+3),2] <- Sys.time()
          ELBOcounter <- ELBOcounter + 4
          model$totalmoves <- model$totalmoves + 1
        }
        
        
        
        
      } else {
        stop("Laps is not of a valid format. Please enter an integer value or 'all' if no merge/delete required.")
      }
      
    }
    end.time2 <- Sys.time()
    total_time <- difftime(end.time2, start.time2, unit = "secs") + total_time1
    all_models[[k]] <- list(ELBO = ELBO[1:(2*iter)], secondELBO = secondELBO[1:ELBOcounter,], Cl = Cl[1:iter], model = model, time = total_time, prior = prior)
    
  }
  names(all_models) <- c(laps)
  return(all_models)
  
}
