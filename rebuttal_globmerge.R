#Simulations for global merges - N = 20000, split into 4 equal batches of 5000 - second simulation

library(tidyverse)
library(mclust)
library(mcclust)
library(mcclust.ext)
source("mixtureMergeDeleteTime2.R")
source("GenerateSampleData.R")
source("globalMerge.R")
source("globalMergeAlt.R")
source('mixtureParallelAlt.R')

set.seed(34643)

library(foreach)
library(doParallel)
library(doRNG)
registerDoParallel(4)

globalmergeSim5again <- list()
for(k in 1:5) { #10 simulated datasets
  
  sampledata <- GenerateSampleData(20000, 12, c(0.05, 0.05, 0.05, 0.05, 0.05, 0.1, 0.15, 0.1, 0.1, 0.2, 0.04, 0.06), 100, 0)
  data <- sampledata[[1]]
  
  all_results_full <- data.frame(run_num = c(1:5), total_time = NA, ARI = NA, clusters = NA, model = rep("full", 5), merge_time = NA, mode_time = NA, mean_ARI = NA)
  all_results_shuffle1 <- data.frame(run_num = c(1:5), total_time = NA, ARI = NA, clusters = NA, model = rep("shuffle1", 5), merge_time = NA, mode_time = NA, mean_ARI = NA) #1 is seq
  all_results_shuffle2 <- data.frame(run_num = c(1:5), total_time = NA, ARI = NA, clusters = NA, model = rep("shuffle2", 5), merge_time = NA, mode_time = NA, mean_ARI = NA) #2 is random
  all_results_par <- data.frame(run_num = c(1:5), total_time = NA, ARI = NA, clusters = NA, model = rep("par", 5), merge_time = NA, mode_time = NA, mean_ARI = NA) #2 is rando
  
  for (w in 1:5){ #5 different shuffles of the data
    
    #Normal
    mainmix <- mixturemodel_merdelTIME(data, 15, 0.01, 1000, 0.000005, 5) #Run merge delete model 5 times for good measure
    all_results_full[w, 2] <- mainmix[[1]]$time
    all_results_full[w, 3] <- adjustedRandIndex(mainmix[[1]]$model$labels, sampledata[[2]]$Cluster)
    all_results_full[w, 4] <- length(unique(mainmix[[1]]$model$labels))
    
    #Parallel
    mainmix2 <- mixturemodel_merdelPar(data, 15, 0.01, 1000, 0.000005, 5) #Run merge delete model 5 times for good measure
    all_results_par[w, 2] <- mainmix2[[1]]$time
    all_results_par[w, 3] <- adjustedRandIndex(mainmix2[[1]]$model$labels, sampledata[[2]]$Cluster)
    all_results_par[w, 4] <- length(unique(mainmix2[[1]]$model$labels))
    
    
    #Shuffle (proposed model)
    
    row_shuffle <- sample(1:nrow(data))
    label_shuffle <- sampledata[[2]]$Cluster[row_shuffle]
    data_split <-  data[row_shuffle, ]
    
    data_split <- list(data_split[1:5000,], data_split[5001:10000,], data_split[10001:15000,], data_split[15001:20000,])
    
    shuffle_time <- 0
    mix_times <- c()
    indiv_ARI <- c()
    
    start.timeshuffle <- Sys.time()
    shuffled_results <- foreach(i = 1:4) %dorng% {
      mix <- mixturemodel_merdelTIME(data_split[[i]], 15, 0.01, 1000, 0.000005, 5)
      mix[[1]]
    }
    end.timeshuffle <- Sys.time()
    
    shuffled_results2 <- lapply(shuffled_results, '[[', 4)
    mix_times <- unlist(lapply(shuffled_results, '[[', 5))
    prior <- shuffled_results[[1]]$prior
    #indiv_ARI <- lapply(shuffled_results, '[[', 3)
    #for (i in 1:4){
      #indiv_ARI[[i]] <- adjustedRandIndex(indiv_ARI[[i]]$labels, label_shuffle[((i - 1) * 5000 + 1) :  (i * 5000)])
    #}
    #indiv_ARI <- unlist(indiv_ARI)
    
    shuffle_time <- shuffle_time + difftime(end.timeshuffle, start.timeshuffle, units = 'secs')
    

    all_results_shuffle1[w, 7] <- max(mix_times)
    all_results_shuffle2[w, 7] <- max(mix_times)
    #all_results_shuffle1[w, 8] <- mean(indiv_ARI)
    #all_results_shuffle2[w, 8] <- mean(indiv_ARI)
    
    
    start.time <- Sys.time()
    global <- globalMerge(shuffled_results2, prior, 60)
    end.time <- Sys.time()
    
    overallLabels <- unlist(global$labels)
    
    all_results_shuffle1[w, 2] <- shuffle_time + difftime(end.time, start.time, units = 'secs')
    all_results_shuffle1[w, 3] <- adjustedRandIndex(overallLabels, label_shuffle)
    all_results_shuffle1[w, 4] <- length(unique(overallLabels))
    all_results_shuffle1[w, 6] <- difftime(end.time, start.time, units = 'secs')
    
    start.time <- Sys.time()
    global2 <- globalMergeAlt(shuffled_results2, prior, 60)
    end.time <- Sys.time()
    
    overallLabels <- unlist(global2$labels)
    
    all_results_shuffle2[w, 2] <- shuffle_time + difftime(end.time, start.time, units = 'secs')
    all_results_shuffle2[w, 3] <- adjustedRandIndex(overallLabels, label_shuffle)
    all_results_shuffle2[w, 4] <- length(unique(overallLabels))
    all_results_shuffle2[w, 6] <- difftime(end.time, start.time, units = 'secs')
    
  }
  
  all_results <- rbind(all_results_shuffle1, all_results_shuffle2, all_results_full, all_results_par)
  
  globalmergeSim5again[[k]] <- all_results
  
}

save(globalmergeSim5again, file = 'globalmergeSim5again.RData')
