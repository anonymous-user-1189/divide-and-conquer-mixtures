#Global merge - variable selection sim, N = 50000

library(tidyverse)
library(mclust)
library(mcclust)
library(mcclust.ext)
source("mixtureMergeDeleteTimeVarSel2.R")
source("GenerateSampleData.R")
source("globalMerge.R")
source("globalMergeAlt.R")

set.seed(5420)

library(foreach)
library(doParallel)
library(doRNG)
registerDoParallel(20)


globalmergeSimVS2 <- foreach(k = 1:20) %dorng% { 
  
  sampledata <- GenerateSampleData(50000, 10, rep(0.1, 10), 75, 25)
  data <- sampledata[[1]]
  
  all_results_full <- data.frame(run_num = c(1:5), total_time = NA, ARI = NA, clusters = NA, model = rep("full", 5), merge_time = NA, mode_time = NA, mean_ARI = NA)
  all_results_shuffle1 <- data.frame(run_num = c(1:5), total_time = NA, ARI = NA, clusters = NA, model = rep("shuffle1", 5), merge_time = NA, mode_time = NA, mean_ARI = NA) #1 is seq
  
  for (w in 1:5){ #5 different shuffles of the data
    
    mainmix <- mixturemodel_merdelTIMEVarSel(data, 20, 0.01, 2, 1000, 0.000005, 5) 
    all_results_full[w, 2] <- mainmix[[1]]$time
    all_results_full[w, 3] <- adjustedRandIndex(mainmix[[1]]$model$labels, sampledata[[2]]$Cluster)
    all_results_full[w, 4] <- length(unique(mainmix[[1]]$model$labels))
    
    row_shuffle <- sample(1:nrow(data))
    label_shuffle <- sampledata[[2]]$Cluster[row_shuffle]
    data_split <-  data[row_shuffle, ]
    
    data_split <- list(data_split[1:10000,], data_split[10001:20000,], data_split[20001:30000,], data_split[30001:40000,], data_split[40001:50000,])
    
    shuffle_time <- 0
    shuffled_results <- list()
    mix_times <- c()
    indiv_ARI <- c()
    for(i in 1:5){
      mix <- mixturemodel_merdelTIMEVarSel(data_split[[i]], 20, 0.01, 2, 1000, 0.000005, 5)
      shuffled_results[[i]] <- mix[[1]]$model
      shuffle_time <- shuffle_time + mix[[1]]$time
      mix_times[i] <- mix[[1]]$time
      indiv_ARI[i] <- adjustedRandIndex(mix[[1]]$model$labels, label_shuffle[((i - 1) * 10000 + 1) :  (i * 10000)])
      prior <- mix[[1]]$prior #Same for each mixutre model anyway
    }
    
    
    all_results_shuffle1[w, 7] <- max(mix_times)
    all_results_shuffle1[w, 8] <- mean(indiv_ARI)
    
    all_results_shuffle2[w, 7] <- max(mix_times)
    all_results_shuffle2[w, 8] <- mean(indiv_ARI)
    
    
    start.time <- Sys.time()
    global <- globalMerge(shuffled_results, prior)
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

save(globalmergeSimVS2, file = 'globalmergeSimVS2.RData')
