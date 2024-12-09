#Brief code for the benchmarking results - dataset used as 'X' was a generated dataset with N=20000 and the same parameters used in 'rebuttal_globmerge'
#Comparing expectStep with parallelisation for calculating the rnk values vs usual expectStep

source("VariationalEStep.R") #usual full E step

expectStepPar <- function(X, model, prior, Kinit){
  modelpar <- lapply(1:4, function(i) {
    start_idx <- (i - 1) * 5000 + 1
    end_idx <- i * 5000
    
    list(
      alpha = model$alpha,
      eps = model$eps,
      labels = NULL,
      totalmoves = model$totalmoves,
      merges = model$merges,
      deletes = model$deletes,
      rnk = NULL
    )
  })
  
  modelpar2 <- foreach(k = 1:4) %dorng% {
    start_idx <- (k - 1) * 5000 + 1
    end_idx <- k * 5000
    model = expectStep(X[start_idx:end_idx,], modelpar[[k]], prior, Kinit)
    model
  }
  
  
  model = list(alpha = modelpar2[[1]]$alpha,
               eps = modelpar2[[1]]$eps,
               labels = unlist(lapply(modelpar2, '[[', 3)),
               totalmoves = modelpar2[[1]]$totalmoves,
               merges = modelpar2[[1]]$merges,
               deletes = modelpar2[[1]]$deletes,
               rnk = rbind(modelpar2[[1]]$rnk, modelpar2[[2]]$rnk, modelpar2[[3]]$rnk, modelpar2[[4]]$rnk)) 
  
}

benchmark_results <- microbenchmark(
  Method1 = expectStepPar(X, model, prior, Kinit),             # Code block 1
  Method2 = expectStep(X, model, prior, Kinit),     # Code block 2
  times = 100  # Number of iterations
)



