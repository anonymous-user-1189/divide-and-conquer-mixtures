######
#Implement MerDel with clusters selected for merging via correlation, time measured for each 'laps' model. *NO laps = 0 SIMULATION*
#'mixturemodel_merdelTIME' (make sure this isn't confused with simulation including laps = 0)
#####


library(klaR)
library(tidyverse)
library(matrixStats)
source('VariationalEStep.R')
source('VariationalMStep.R')
source('deleteMove.R')
source('mergeMove.R')
source('ELBOCalcMStep.R')

library(Rcpp)
library(RcppArmadillo)
Rcpp::sourceCpp('MixModel.cpp')
Rcpp::sourceCpp('MixModelArmaNew.cpp')

#Variational mixtures for categorical distributions
check_convergence<- function(ELBO, iter, maxiter, tol){
  if (iter > 1 && abs(ELBO[iter] - ELBO[iter - 1]) < tol && abs(ELBO[iter-1] - ELBO[iter-2] < tol) && 
      abs(ELBO[iter-2] - ELBO[iter-3] < tol)){
    cat("Stopped after iteration",iter, "\n", sep = " ") #make sure the last 3 ELBOs close to each other
    return(TRUE)
  }
  if (iter == maxiter){
    cat("Not converged after maximum number of iterations")
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
#laps = number of E-M iterations to do before performing merge/delete. Numeric, set higher than maxiter if you don't want merge/delete

mixturemodel_merdelTIME <- function(data, Kinit, alpha, maxiter, tol, laps, verbose = FALSE){
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
  K = Kinit #K is the current number of clusters (start with all of them), Kinit is the original number - keep prior consistent 
  
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
  
    model = list(alpha = rep(prior$alpha, Kinit),
               eps = firstepsCalc(Kinit, maxNCat, D, N, prior$eps, X, clusterInit),
               labels = clusterInit,
               totalmoves = 0,
               merges = 0,
               deletes = 0) 
  
  
  ########################################
    ELBOcounter <- 1
  
    lap <- laps[k]
    for (iter in 1:maxiter){
      if (lap %% 1 == 0){
        
        model = expectStep(X, model, prior, Kinit) #Expectation step
        model = maxStep(X, model, prior) #Maximisation step
        ELBO[iter] = secondELBO[ELBOcounter,1] = ELBOCalcMStep(X, model, prior, Kinit) #ELBO
        Cl[iter] = length(unique(model$labels)) #Counts number of non-empty clusters
        secondELBO[ELBOcounter,2] <- Sys.time()
        ELBOcounter <- ELBOcounter + 1
        
        if(verbose){
          cat("Iteration number ", iter, " ELBO : ", ELBO[iter], "\n", sep = "")
        }
        
        if(check_convergence(ELBO, iter, maxiter, tol)) break
        
        #Merge delete: only happens every "laps" number of EM iterations
        if (iter %% lap == 0){
          mergeResult = mergeMove(X, model, prior, secondELBO[ELBOcounter - 1,1], Kinit)
          model = mergeResult$model
          secondELBO[ELBOcounter:(ELBOcounter+3),1] = rep(mergeResult$ELBO, 4)
          secondELBO[ELBOcounter:(ELBOcounter+3),2] <- Sys.time()
          ELBOcounter <- ELBOcounter + 4
          
          deleteResult = deleteMove(X, model, prior, secondELBO[ELBOcounter - 1,1], Kinit)
          model = deleteResult$model
          secondELBO[ELBOcounter:(ELBOcounter+3),1] = rep(deleteResult$ELBO, 4)
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
    all_models[[k]] <- list(ELBO = ELBO[1:(iter)], secondELBO = secondELBO[1:ELBOcounter,], Cl = Cl[1:iter], model = model, time = total_time, prior = prior)
    
  }
  
  names(all_models) <- c(laps)
  return(all_models)
  
}
