#Simulations for global merges - N = 10000, split into 5 equal batches of 2000 - first simulation

library(tidyverse)
library(mclust)
library(mcclust)
library(mcclust.ext)
source("mixtureMergeDeleteTime2.R")
source("GenerateSampleData.R")
source("globalMerge.R")
source("globalMergeAlt.R")

set.seed(34643)

library(foreach)
library(doParallel)
library(doRNG)
registerDoParallel(20)


globalmergeSim4 <- foreach(k = 1:20) %dorng% { #20 simulated datasets
  
  sampledata <- GenerateSampleData(10000, 8, c(0.2, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1), 100, 0)
  data <- sampledata[[1]]
  
  all_results_full <- data.frame(run_num = c(1:10), total_time = NA, ARI = NA, clusters = NA, model = rep("full", 10), merge_time = NA, mode_time = NA, mean_ARI = NA)
  all_results_shuffle1 <- data.frame(run_num = c(1:10), total_time = NA, ARI = NA, clusters = NA, model = rep("shuffle1", 10), merge_time = NA, mode_time = NA, mean_ARI = NA) #1 is seq
  all_results_shuffle2 <- data.frame(run_num = c(1:10), total_time = NA, ARI = NA, clusters = NA, model = rep("shuffle2", 10), merge_time = NA, mode_time = NA, mean_ARI = NA) #2 is random
  
  for (w in 1:10){ #10 different shuffles of the data
    
    mainmix <- mixturemodel_merdelTIME(data, 15, 0.01, 1000, 0.000005, 5) #Run merge delete model 10 times for good measure
    all_results_full[w, 2] <- mainmix[[1]]$time
    all_results_full[w, 3] <- adjustedRandIndex(mainmix[[1]]$model$labels, sampledata[[2]]$Cluster)
    all_results_full[w, 4] <- length(unique(mainmix[[1]]$model$labels))
    
    row_shuffle <- sample(1:nrow(data))
    label_shuffle <- sampledata[[2]]$Cluster[row_shuffle]
    data_split <-  data[row_shuffle, ]
    
    data_split <- list(data_split[1:2000,], data_split[2001:4000,], data_split[4001:6000,], data_split[6001:8000,], data_split[8001:10000,])
    
    shuffle_time <- 0
    shuffled_results <- list()
    mix_times <- c()
    indiv_ARI <- c()
    for(i in 1:5){
      mix <- mixturemodel_merdelTIME(data_split[[i]], 15, 0.01, 1000, 0.000005, 5)
      shuffled_results[[i]] <- mix[[1]]$model
      shuffle_time <- shuffle_time + mix[[1]]$time
      mix_times[i] <- mix[[1]]$time
      indiv_ARI[i] <- adjustedRandIndex(mix[[1]]$model$labels, label_shuffle[((i - 1) * 2000 + 1) :  (i * 2000)])
      prior <- mix[[1]]$prior #Same for each mixutre model anyway
      nCat <- mix[[1]]$nCat #Same for each mixutre model anyway
    }
    
    all_results_shuffle1[w, 7] <- max(mix_times)
    all_results_shuffle2[w, 7] <- max(mix_times)
    all_results_shuffle1[w, 8] <- mean(indiv_ARI)
    all_results_shuffle2[w, 8] <- mean(indiv_ARI)
    
    start.time <- Sys.time()
    global <- globalMerge(shuffled_results, prior, nCat)
    end.time <- Sys.time()
    
    overallLabels <- unlist(global$labels)
    
    all_results_shuffle1[w, 2] <- shuffle_time + difftime(end.time, start.time, units = 'secs')
    all_results_shuffle1[w, 3] <- adjustedRandIndex(overallLabels, label_shuffle)
    all_results_shuffle1[w, 4] <- length(unique(overallLabels))
    all_results_shuffle1[w, 6] <- difftime(end.time, start.time, units = 'secs')
    
    start.time <- Sys.time()
    global2 <- globalMergeAlt(shuffled_results, prior, nCat)
    end.time <- Sys.time()
    
    overallLabels <- unlist(global2$labels)
    
    all_results_shuffle2[w, 2] <- shuffle_time + difftime(end.time, start.time, units = 'secs')
    all_results_shuffle2[w, 3] <- adjustedRandIndex(overallLabels, label_shuffle)
    all_results_shuffle2[w, 4] <- length(unique(overallLabels))
    all_results_shuffle2[w, 6] <- difftime(end.time, start.time, units = 'secs')
    
  }
  
  all_results <- rbind(all_results_shuffle1, all_results_shuffle2, all_results_full)
  all_results
  
}

save(globalmergeSim4, file = 'globalmergeSim4.RData')
