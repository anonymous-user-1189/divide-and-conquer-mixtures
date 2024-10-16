#Compare correlation/divergence/random simulation

library(tidyverse)
library(mclust)
library(mcclust)
library(mcclust.ext)
source("mixtureMergeDeleteTime.R")
source("mixtureMergeDeleteTimeDiv.R")
source("mixtureMergeDeleteTimeRand.R")
source("GenerateSampleData.R")

set.seed(139)

library(foreach)
library(doParallel)
library(doRNG)
registerDoParallel(10)

all_sims2b <- list()
for (a in 1:15){
  generation <- GenerateSampleData(2000, 20, rep(c(0.1, 0.2, 0.3, 0.2, 0.2), 4), 100, 0)
  data <- generation[[1]]
  truelabels <- generation[[2]]
  resultforpsm <- list()
  timeresults <- list()
  mergesdeletes <- list()
  clusters <- list()
  clustlabels <- list()
  
  all_results <- foreach(k = 1:10) %dorng%{
    mix1 <- mixturemodel_merdelTIME(data, 40, 0.01, 1500, 0.00000005, c(1, 5))
    mix2 <- mixturemodel_merdeldiv(data, 40, 0.01, 1500, 0.00000005, c(1, 5), divergence = "KL")
    mix3 <- mixturemodel_merdeldiv(data, 40, 0.01, 1500, 0.00000005, c(1, 5), divergence = "bha")
    mix4 <- mixturemodel_merdelrand(data, 40, 0.01, 1500, 0.00000005, c(1, 5))
    
    #ARI
    merdelLabels1 <- lapply(mix1, "[[", 4)
    merdelLabels1 <- lapply(merdelLabels1, "[[", 3)
    merdelLabels1 <- lapply(merdelLabels1, function(x) adjustedRandIndex(x, truelabels$Cluster))
    merdelLabels1 <- data.frame(
      laps = names(merdelLabels1),
      ARI = unlist(merdelLabels1),
      model = "corr"
    )
    merdelLabels2 <- lapply(mix2, "[[", 4)
    merdelLabels2 <- lapply(merdelLabels2, "[[", 3)
    merdelLabels2 <- lapply(merdelLabels2, function(x) adjustedRandIndex(x, truelabels$Cluster))
    merdelLabels2 <- data.frame(
      laps = names(merdelLabels2),
      ARI = unlist(merdelLabels2),
      model = "KL"
    )
    merdelLabels3 <- lapply(mix3, "[[", 4)
    merdelLabels3 <- lapply(merdelLabels3, "[[", 3)
    merdelLabels3 <- lapply(merdelLabels3, function(x) adjustedRandIndex(x, truelabels$Cluster))
    merdelLabels3 <- data.frame(
      laps = names(merdelLabels3),
      ARI = unlist(merdelLabels3),
      model = "bha"
    )
    merdelLabels4 <- lapply(mix4, "[[", 4)
    merdelLabels4 <- lapply(merdelLabels4, "[[", 3)
    merdelLabels4 <- lapply(merdelLabels4, function(x) adjustedRandIndex(x, truelabels$Cluster))
    merdelLabels4 <- data.frame(
      laps = names(merdelLabels4),
      ARI = unlist(merdelLabels4),
      model = "rand"
    )
    merdelLabels <- rbind(merdelLabels1, merdelLabels2, merdelLabels3, merdelLabels4)
    
    
    #merges
    merdelmdrate1 <- lapply(mix1, "[[", 4)
    merdelmdrate1 <- lapply(merdelmdrate1, function(x) x[4:6])
    merdelmergerate1 <- lapply(merdelmdrate1, function(x) x[[2]] / x[[1]])
    merdelmergerate1 <- data.frame(
      laps = names(merdelmergerate1),
      merge_rate = unlist(merdelmergerate1),
      model = "corr"
    )
    
    merdelmdrate2 <- lapply(mix2, "[[", 4)
    merdelmdrate2 <- lapply(merdelmdrate2, function(x) x[4:6])
    merdelmergerate2 <- lapply(merdelmdrate2, function(x) x[[2]] / x[[1]])
    merdelmergerate2 <- data.frame(
      laps = names(merdelmergerate2),
      merge_rate = unlist(merdelmergerate2),
      model = "KL"
    )
    
    merdelmdrate3 <- lapply(mix3, "[[", 4)
    merdelmdrate3 <- lapply(merdelmdrate3, function(x) x[4:6])
    merdelmergerate3 <- lapply(merdelmdrate3, function(x) x[[2]] / x[[1]])
    merdelmergerate3 <- data.frame(
      laps = names(merdelmergerate3),
      merge_rate = unlist(merdelmergerate3),
      model = "bha"
    )
    
    merdelmdrate4 <- lapply(mix4, "[[", 4)
    merdelmdrate4 <- lapply(merdelmdrate4, function(x) x[4:6])
    merdelmergerate4 <- lapply(merdelmdrate4, function(x) x[[2]] / x[[1]])
    merdelmergerate4 <- data.frame(
      laps = names(merdelmergerate4),
      merge_rate = unlist(merdelmergerate4),
      model = "rand"
    )
    merdelmerges <- rbind(merdelmergerate1, merdelmergerate2, merdelmergerate3, merdelmergerate4)
    
    
    #times
    
    merdeltimes1 <-  lapply(mix1, "[[", 5)
    merdeltimes2 <-  lapply(mix2, "[[", 5)
    merdeltimes3 <-  lapply(mix3, "[[", 5)
    merdeltimes4 <-  lapply(mix4, "[[", 5)
    merdeltimes1 <- data.frame(
      laps = names(merdeltimes1),
      total_time = unlist(merdeltimes1),
      model = "corr"
    )
    merdeltimes2 <- data.frame(
      laps = names(merdeltimes2),
      total_time = unlist(merdeltimes2),
      model = "KL"
    )
    merdeltimes3 <- data.frame(
      laps = names(merdeltimes3),
      total_time = unlist(merdeltimes3),
      model = "bha"
    )
    merdeltimes4 <- data.frame(
      laps = names(merdeltimes4),
      total_time = unlist(merdeltimes4),
      model = "rand"
    )
    merdeltimes <- rbind(merdeltimes1, merdeltimes2, merdeltimes3, merdeltimes4)
    
    #clusters
    merdelClusters1 <- lapply(mix1, "[[", 3)
    merdelClusters1 <- lapply(merdelClusters1, function(x) x[length(x)])
    merdelClusters1 <- data.frame(
      laps = names(merdelClusters1),
      clusters = unlist(merdelClusters1),
      model = "corr"
    )
    
    merdelClusters2 <- lapply(mix2, "[[", 3)
    merdelClusters2 <- lapply(merdelClusters2, function(x) x[length(x)])
    merdelClusters2 <- data.frame(
      laps = names(merdelClusters2),
      clusters = unlist(merdelClusters2),
      model = "KL"
    )
    
    merdelClusters3 <- lapply(mix3, "[[", 3)
    merdelClusters3 <- lapply(merdelClusters3, function(x) x[length(x)])
    merdelClusters3 <- data.frame(
      laps = names(merdelClusters3),
      clusters = unlist(merdelClusters3),
      model = "bha"
    )
    
    merdelClusters4 <- lapply(mix4, "[[", 3)
    merdelClusters4 <- lapply(merdelClusters4, function(x) x[length(x)])
    merdelClusters4 <- data.frame(
      laps = names(merdelClusters4),
      clusters = unlist(merdelClusters4),
      model = "KL"
    )
    merdelClusters <- rbind(merdelClusters1, merdelClusters2, merdelClusters3, merdelClusters4)
    
    init_results <- list(merdelLabels, merdelmerges, merdeltimes, merdelClusters)
    init_results
  }  
  
  labels <- lapply(all_results, "[[", 1)
  labels <- bind_rows(labels, .id = "init_num")
  
  merges <- lapply(all_results, "[[", 2)
  merges <- bind_rows(merges, .id = "init_num")
  
  times <- lapply(all_results, "[[", 3)
  times <- bind_rows(times, .id = "init_num")
  
  clusters <- lapply(all_results, "[[", 4)
  clusters <- bind_rows(clusters, .id = "init_num")
  
  all_sims2b[[a]] <- list(labels, merges, times, clusters)
  
}


save(all_sims2b, file="all_sims2b.RData")


