#MNIST test run

#library(tidyverse)
library(mclust)
load("mnist_test_binnew.RData")
load("mnist_test_labels.RData")
source("mixtureMergeDeleteTime2.R")

set.seed(344)

library(foreach)
library(doParallel)
library(doRNG)
registerDoParallel(5)

full_results <- list()
all_results_test1 <- data.frame(run_num = c(1:5), ARI = NA, clusters = NA, model = rep(20, 5))
all_results_test2 <- data.frame(run_num = c(1:5), ARI = NA, clusters = NA, model = rep(35, 5))

shuffled_results <- foreach(k = 1:5)%dorng%{
  mix1 <- mixturemodel_merdelTIME(mnist_test_binnew, 20, 0.01, 1000, 0.000005, 5)
  mix2 <- mixturemodel_merdelTIME(mnist_test_binnew, 35, 0.01, 1000, 0.000005, 5)

  all_results <- list(mix1[[1]], mix2[[1]])
  all_results
}

for(k in 1:5){
  all_results_test1[k, 2] <- adjustedRandIndex(shuffled_results[[k]][[1]]$model$labels, mnist_test_labels$V1)
  all_results_test1[k, 3] <- length(unique(shuffled_results[[k]][[1]]$model$labels))
  
  all_results_test2[k, 2] <- adjustedRandIndex(shuffled_results[[k]][[2]]$model$labels, mnist_test_labels$V1)
  all_results_test2[k, 3] <- length(unique(shuffled_results[[k]][[2]]$model$labels))
  
  full_results[[k]] <- list(twenty = shuffled_results[[k]][[1]]$model$labels, thirtyfive = shuffled_results[[k]][[2]]$model$labels)
}


shuffles <- rbind(all_results_test1, all_results_test2)
mnist_bin_test_results <- list(shuffles, full_results)

save(mnist_bin_test_results, file = "mnist_bin_test_results.RData")


