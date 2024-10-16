#Global merge - categorical sim, N = 20000

library(tidyverse)
library(mclust)
library(mcclust)
library(mcclust.ext)
source("mixtureMergeDeleteTimeDiv2.R")
source("GenerateSampleDataCat.R")
source("globalMergeDiv.R")

set.seed(5425)

library(foreach)
library(doParallel)
library(doRNG)
registerDoParallel(20)


globalmergeSimCat1 <- foreach(k = 1:20) %dorng% { #20 simulated datasets
  
  sampledata <- GenerateSampleDataCat(20000, 20, rep(0.1, 20), 100, 0, yout = FALSE, cat = 3)
  data <- sampledata[[1]]
  
  all_results_full <- data.frame(run_num = c(1:10), total_time = NA, ARI = NA, clusters = NA, model = rep("full", 10), merge_time = NA, mode_time = NA, mean_ARI = NA)
  all_results_shuffle1 <- data.frame(run_num = c(1:10), total_time = NA, ARI = NA, clusters = NA, model = rep("shuffle1", 10), merge_time = NA, mode_time = NA, mean_ARI = NA) #1 is seq

  for (w in 1:10){ #10 different shuffles of the data
    
    mainmix <- mixturemodel_merdeldiv(data, 30, 0.01, 1000, 0.000005, 5, "KL") 
    all_results_full[w, 2] <- mainmix[[1]]$time
    all_results_full[w, 3] <- adjustedRandIndex(mainmix[[1]]$model$labels, sampledata[[2]]$Cluster)
    all_results_full[w, 4] <- length(unique(mainmix[[1]]$model$labels))
    
    row_shuffle <- sample(1:nrow(data))
    label_shuffle <- sampledata[[2]]$Cluster[row_shuffle]
    data_split <-  data[row_shuffle, ]
    
    data_split <- list(data_split[1:5000,], data_split[5001:10000,], data_split[10001:15000,], data_split[15001:20000,])
    
    shuffle_time <- 0
    shuffled_results <- list()
    mix_times <- c()
    indiv_ARI <- c()
    for(i in 1:4){
      mix <- mixturemodel_merdeldiv(data_split[[i]], 30, 0.01, 1000, 0.000005, 5, "KL")
      shuffled_results[[i]] <- mix[[1]]$model
      shuffle_time <- shuffle_time + mix[[1]]$time
      mix_times[i] <- mix[[1]]$time
      indiv_ARI[i] <- adjustedRandIndex(mix[[1]]$model$labels, label_shuffle[((i - 1) * 5000 + 1) :  (i * 5000)])
      prior <- mix[[1]]$prior #Same for each mixutre model anyway
    }
    
    
    all_results_shuffle1[w, 7] <- max(mix_times)
    all_results_shuffle1[w, 8] <- mean(indiv_ARI)
    
    
    start.time <- Sys.time()
    global <- globalMergeDiv(shuffled_results, prior)
    end.time <- Sys.time()
    
    overallLabels <- unlist(global$labels)
    
    all_results_shuffle1[w, 2] <- shuffle_time + difftime(end.time, start.time, units = 'secs')
    all_results_shuffle1[w, 3] <- adjustedRandIndex(overallLabels, label_shuffle)
    all_results_shuffle1[w, 4] <- length(unique(overallLabels))
    all_results_shuffle1[w, 6] <- difftime(end.time, start.time, units = 'secs')
    
  }
  
  all_results <- rbind(all_results_shuffle1, all_results_full)
  all_results
  
}

save(globalmergeSimCat1, file = 'globalmergeSimCat1.RData')
