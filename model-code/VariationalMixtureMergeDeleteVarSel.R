######
#Main function to implement MerDel with variable selection: mixturemodel_merdelvarsel
#####

library(klaR)
library(tidyverse)
library(matrixStats)
source('VariationalEStepVarSel.R')
source('VariationalMStepVarSel.R')
source('deleteMoveVarSel.R')
source('mergeMoveVarSel.R')
source('ELBOCalcVarSel.R')
source('ELBOCalcVarSelMStep.R')

library(Rcpp)
library(RcppArmadillo)
Rcpp::sourceCpp('VarSel.cpp')
Rcpp::sourceCpp('VarSelArma.cpp')

#Check convergence function
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

mixturemodel_merdelvarsel <- function(data, K, alpha, a, maxiter, tol, laps){
  ################################################
  #Process dataset
  
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
  categories <- lapply(1:ncol(X), function(j)levels(X[, j])) #list out categories
  
  cat_lengths <- sapply(categories, length)
  if(any(cat_lengths == 0)) {
    stop("Column(s) ", paste0(colnames(X)[cat_lengths == 0], collapse = ","), " is empty (i.e. all values are NA). Please fix this.")
  }
  #Currently work under assumption there are no missing values? unsure how to deal with this now lol
  
  #Create numeric matrix for data
  X <- data.matrix(X)
  #note that now eg. binary factors are now 1 and 2 instead of 0 and 1
  
  ##########################################
  
  N = dim(X)[1]
  D = dim(X)[2]
  
  nCat <- as.vector(apply(X, 2, max)) #number of categories in each variable
  maxNCat <- max(nCat)
  ELBO = Cl = secondELBO = rep(-Inf,maxiter*2) #prepare for L (ELBO) and number of clusters at each iteration to be recorded
  
  ########################################
  #Define priors
  #Assume there is an input alpha and all variables have a symmetric Dirichlet prior with parameters alpha
  
  prior = list(alpha = alpha) #default value of alpha - prior for clustering pi
  prior$eps = matrix(0, D, maxNCat) #prior for clustering phi
  for(d in 1:D){
    prior$eps [d,1:nCat[d]] <- 1/nCat[d]
  }
  prior$a = a
  
  #Initialising
  #cluster = sample(1:K, N, replace=T) random assignment!
  clusterInit <- klaR::kmodes(X, modes = K)$cluster
  
  #######################################
  
  EPSreshape = t(prior$eps) 
  dim(EPSreshape) = c(1,maxNCat,D) #sets dimension of matrix, D arrays (variable) where each one has 1 row and 
  #maxNCat columns (parameter for each category in the variable, all unused are 0)
  model = list(alpha = rep(prior$alpha, K),
               eps = EPSreshape[rep(1,K),,],
               c = rep(1, D), #initialise c_i - all variables initially included
               labels = clusterInit,
               totalmoves = 0,
               merges = 0,
               deletes = 0) #current parameters for alpha, epsilon: makes D arrays
  #ie. array for each variable, each array has a row for each cluster and a column for each category
  #each array represents phi_ki k*max(Li) number of columns and rows, zeroes to stand for categories unused
  #model alpha and model epsilon is therefore pi parameters and phi parameters
  
  #######################################
  
  #Just setting null phi to the rate of the parameter value in the dataset - max likelihood?
  nullphi <- array(0, dim = c(D, maxNCat))
  for (i in 1:D){
    for (j in 1:nCat[i]){
      nullphi[i, j] <- sum((X[,i] == j)) / N
    }
  }
  model$nullphi <- nullphi
  
  #######################################
  model$eps <- firstepsCalc(K, maxNCat, D, N, prior$eps, X, clusterInit)
  
  
  ########################################
  ELBOcounter <- 1
  for (iter in 1:maxiter){
    
    
    if (laps == "all"){
      print(paste("Iteration number ",iter))
      model = expectStep(X, model) #Expectation step
      .GlobalEnv$maxNCat <- maxNCat
      ELBO[iter * 2-1] = ELBOCalc(X, model, prior) #ELBO
      print(ELBO[iter * 2-1])
      
      Cl[iter] = length(unique(model$labels)) #Counts number of non-empty clusters
      
      model = maxStep(X, model, prior) #Maximisation step
      ELBO[iter * 2] = ELBOCalcMStep(X, model, prior) #ELBO
      print(ELBO[iter * 2])
      
      Cl[iter] = length(unique(model$labels)) #Counts number of non-empty clusters
      
      if(check_convergence(ELBO, iter, maxiter, tol)) break
      
    } else if (laps %% 1 == 0){
      
      #E step
      print(paste("Iteration number ",iter))
      model = expectStepVarSel(X, model) #Expectation step
      .GlobalEnv$maxNCat <- maxNCat
      ELBO[iter * 2 - 1] = secondELBO[ELBOcounter] = ELBOCalcVarSel(X, model, prior) #ELBO
      print(secondELBO[ELBOcounter])
      ELBOcounter <- ELBOcounter + 1
      
      Cl[iter] = length(unique(model$labels)) #Counts number of non-empty clusters
      
      #M step
      model = maxStepVarSel(X, model, prior) #Maximisation step
      ELBO[iter * 2] = secondELBO[ELBOcounter] = ELBOCalcVarSelMStep(X, model, prior) #ELBO
      print(secondELBO[ELBOcounter])
      ELBOcounter <- ELBOcounter + 1
      
      Cl[iter] = length(unique(model$labels)) #Counts number of non-empty clusters
      
      if(check_convergence(ELBO, iter, maxiter, tol)) break
      
      #Merge delete: only happens every "laps" number of EM iterations
      if (iter %% laps == 0){
        model = mergeMoveVarSel(X, model, prior, secondELBO[ELBOcounter - 1])
        secondELBO[ELBOcounter:(ELBOcounter+3)] = rep(ELBOCalcVarSelMStep(X, model, prior), 4)
        ELBOcounter <- ELBOcounter + 4
        model = deleteMoveVarSel(X, model, prior, secondELBO[ELBOcounter - 1])
        secondELBO[ELBOcounter:(ELBOcounter+3)] = rep(ELBOCalcVarSelMStep(X, model, prior), 4)
        ELBOcounter <- ELBOcounter + 4
        model$totalmoves <- model$totalmoves + 1
      }
      
      
      
      
    } else {
      stop("Laps is not of a valid format. Please enter an integer value or 'all' if no merge/delete required.")
    }
  }
  
  output <- list(ELBO = ELBO[1:(2*iter)], secondELBO = secondELBO[1:ELBOcounter], Cl = Cl[1:iter], model = model, prior = prior)
  return(output)
  
}
