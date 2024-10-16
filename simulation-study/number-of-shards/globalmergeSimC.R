#Simulations for global merges - N = 200000

library(tidyverse)
library(mclust)
library(mcclust)
library(mcclust.ext)
source("mixtureMergeDeleteTime2.R")
source("GenerateSampleData.R")
source("globalMerge.R")
source("globalMergeAlt.R")

set.seed(12)

library(foreach)
library(doParallel)
library(doRNG)
registerDoParallel(20)


globalmergeSimC <- foreach(k = 1:20) %dorng% { #20 simulated datasets
  
  sampledata <- GenerateSampleData(200000, 10, c(0.1, 0.1, 0.05, 0.05, 0.05, 0.1, 0.15, 0.1, 0.1, 0.2), 100, 0)
  data <- sampledata[[1]]
  
  all_results_shuffle1 <- data.frame(run_num = c(1:5), model = rep("shuffle1", 10), ARI_4 = NA, clusters_4 = NA, ARI_8 = NA, clusters_8 = NA) #1 is seq
  all_results_shuffle2 <- data.frame(run_num = c(1:5), model = rep("shuffle2", 10), ARI_4 = NA, clusters_4 = NA, ARI_8 = NA, clusters_8 = NA) #1 is seq
  
  for (w in 1:5){ #10 different shuffles of the data
    
    row_shuffle1 <- sample(1:nrow(data))
    label_shuffle1 <- sampledata[[2]]$Cluster[row_shuffle1]
    data_split1 <-  data[row_shuffle1, ]
    
    data_split1 <- list(data_split1[1:50000,], data_split1[50001:100000,], data_split1[100001:150000,], data_split1[150001:200000,])
    
    shuffled_results <- list()
    mix_times <- c()
    indiv_ARI <- c()
    for(i in 1:4){
      mix <- mixturemodel_merdelTIME(data_split1[[i]], 20, 0.01, 1000, 0.000005, 5)
      shuffled_results[[i]] <- mix[[1]]$model
      prior <- mix[[1]]$prior #Same for each mixutre model anyway
    }
    
    start.time <- Sys.time()
    global <- globalMerge(shuffled_results, prior)
    end.time <- Sys.time()
    
    overallLabels <- unlist(global$labels)
    
    all_results_shuffle1[w, 3] <- adjustedRandIndex(overallLabels, label_shuffle1)
    all_results_shuffle1[w, 4] <- length(unique(overallLabels))
    
    start.time <- Sys.time()
    global2 <- globalMergeAlt(shuffled_results, prior)
    end.time <- Sys.time()
    
    overallLabels <- unlist(global2$labels)
    
    all_results_shuffle2[w, 3] <- adjustedRandIndex(overallLabels, label_shuffle1)
    all_results_shuffle2[w, 4] <- length(unique(overallLabels))
    
    
    row_shuffle2 <- sample(1:nrow(data))
    label_shuffle2 <- sampledata[[2]]$Cluster[row_shuffle2]
    data_split2 <-  data[row_shuffle2, ]
    
    data_split2 <- list(data_split2[1:25000,], data_split2[25001:50000,], data_split2[50001:75000,], data_split2[75001:100000,], data_split2[100001:125000,], data_split2[125001:150000,], data_split2[150001:175000,], data_split2[175001:200000,])
    
    shuffled_results <- list()
    mix_times <- c()
    indiv_ARI <- c()
    for(i in 1:8){
      mix <- mixturemodel_merdelTIME(data_split2[[i]], 20, 0.01, 1000, 0.000005, 5)
      shuffled_results[[i]] <- mix[[1]]$model
      prior <- mix[[1]]$prior #Same for each mixutre model anyway
    }
    
    start.time <- Sys.time()
    global <- globalMerge(shuffled_results, prior)
    end.time <- Sys.time()
    
    overallLabels <- unlist(global$labels)
    
    all_results_shuffle1[w, 5] <- adjustedRandIndex(overallLabels, label_shuffle2)
    all_results_shuffle1[w, 6] <- length(unique(overallLabels))
    
    start.time <- Sys.time()
    global2 <- globalMergeAlt(shuffled_results, prior)
    end.time <- Sys.time()
    
    overallLabels <- unlist(global2$labels)
    
    all_results_shuffle2[w, 5] <- adjustedRandIndex(overallLabels, label_shuffle2)
    all_results_shuffle2[w, 6] <- length(unique(overallLabels))
    
  }
  
  all_results <- rbind(all_results_shuffle1, all_results_shuffle2)
  all_results
  
  
}

save(globalmergeSimC, file = 'globalmergeSimC.RData')
