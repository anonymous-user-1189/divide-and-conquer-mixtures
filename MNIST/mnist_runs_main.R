#MNIST global merge and results - K_init = 20 (main results in paper)

#library(tidyverse)
library(mclust)
load("mnist_binnew.RData")
load("mnist_labels.RData")
source("mixtureMergeDeleteTime2.R")
source("globalMerge.R")
source("globalMergeAlt.R")

set.seed(383)

library(foreach)
library(doParallel)
library(doRNG)
registerDoParallel(6)

full_results <- list()
all_results_shuffle2 <- data.frame(run_num = c(1:10), total_time = NA, ARI = NA, clusters = NA, ELBO = NA, model = rep("shuffle2", 10))
all_results_shuffle1 <- data.frame(run_num = c(1:10), total_time = NA, ARI = NA, clusters = NA, ELBO = NA, model = rep("shuffle1", 10)) #1 is seq

counter = 1
while (counter < 6){
  
  print("counter")
  print(counter)
  
  row_shuffle <- sample(1:60000)
  label_shuffle <- mnist_labels[row_shuffle,1]
  mnist_split <-  mnist_binnew[row_shuffle, ]
  
  mnist_split <- list(mnist_split[1:10000,], mnist_split[10001:20000,], mnist_split[20001:30000,], mnist_split[30001:40000,], mnist_split[40001:50000,], mnist_split[50001:60000,])
  
  category_count <- lapply(mnist_split, function(x) apply(x, 2, function(col) length(unique(col))))
  all_categories <- unlist(category_count)
  if(sum(all_categories == 1)){
    print("broken")
    break #make sure there are no blank columns
  }
  
  samplemix <- mixturemodel_merdelTIME(mnist_split[[1]], 20, 0.01, 2, 0.000005, 5)
  prior <- samplemix[[1]]$prior #Same for each mixutre model anyway
  nCat <- samplemix[[1]]$nCat #Same for each mixutre model anyway
  
  start.time <- Sys.time()
  
  shuffled_results <- foreach(k = 1:6)%dorng%{
    mix <- mixturemodel_merdelTIME(mnist_split[[k]], 20, 0.01, 1000, 0.000005, 5)
    mix[[1]]$model
  }
  
  start.time2 <- Sys.time()
  
  global <- globalMerge(shuffled_results, prior)
  
  print("Global Merge 1 complete")
  
  end.time <- Sys.time()
  
  global2 <- globalMergeAlt(shuffled_results, prior)
  
  print("Global Merge 2 complete")
  
  end.time2 <- Sys.time()
  
  overallLabels <- unlist(global$labels)
  overallLabels2 <- unlist(global2$labels)
  
  all_results_shuffle1[counter, 3] <- adjustedRandIndex(overallLabels, label_shuffle)
  all_results_shuffle1[counter, 4] <- length(unique(overallLabels))
  all_results_shuffle1[counter, 2] <- difftime(end.time, start.time, units = 'secs')
  all_results_shuffle1[counter, 5] <- global$ELBO
  
  all_results_shuffle2[counter, 3] <- adjustedRandIndex(overallLabels2, label_shuffle)
  all_results_shuffle2[counter, 4] <- length(unique(overallLabels2))
  all_results_shuffle2[counter, 2] <- difftime(end.time2, end.time, units = 'secs') + difftime(start.time2, start.time, units = 'secs')
  all_results_shuffle2[counter, 5] <- global2$ELBO
  
  full_results[[counter]] <- list(shuffle1 = overallLabels, shuffle2 = overallLabels2, row_shuffle = row_shuffle, eps1 = global$eps, eps2 = global2$eps)
  
  counter = counter + 1
}

shuffles <- rbind(all_results_shuffle1, all_results_shuffle2)
mnist_bin_results <- list(shuffles, full_results)

save(mnist_bin_results, file = "mnist_bin_results.RData")




