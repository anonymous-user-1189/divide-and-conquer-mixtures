#Frequency of merge/delete with categorical data simulations - N=2000, P=100, 8 clusters, 4 categories

library(tidyverse)
library(mclust)
library(mcclust)
library(mcclust.ext)
source("mixtureMergeDeleteTimeDiv.R")
source("GenerateSampleDataCat.R")

set.seed(35001)

library(foreach)
library(doParallel)
library(doRNG)
registerDoParallel(20)

merdelsim1Cat <- foreach(k = 1:20) %dorng% {
  generation <- GenerateSampleDataCat(2000, 8, c(0.025, 0.025, 0.05, 0.1, 0.1, 0.4, 0.1, 0.2), 100, 0, yout = FALSE, cat = 4)
  data <- generation[[1]]
  truelabels <- generation[[2]]
  resultforpsm <- list()
  ELBOresults <- list()
  timeresults <- list()
  mergesdeletes <- list()
  clusters <- list()
  clustlabels <- list()
  for (j in 1:10){
    mix <- mixturemodel_merdeldiv(data, 20, 0.01, 1500, 0.00000005, c(1, 2,5, 10, 10000), divergence = "KL")
    
    merdelELBO <- lapply(mix, "[[", 2)
    merdelELBO <- lapply(merdelELBO, function(df) {
      df[, 2] <- difftime(df$Time, df$Time[1], unit = 'secs')
      return(df)
    })
    merdelELBO <- bind_rows(merdelELBO, .id = "laps")
    
    merdelLabels <- lapply(mix, "[[", 4)
    merdelLabels <- lapply(merdelLabels, "[[", 3)
    merdelLabels <- lapply(merdelLabels, function(x) adjustedRandIndex(x, truelabels$Cluster))
    merdelLabels <- data.frame(
         laps = names(merdelLabels),
         ARI = unlist(merdelLabels)
    )
    
    merdelclustLabels <- lapply(mix, "[[", 4)
    merdelclustLabels <- lapply(merdelclustLabels, "[[", 3)
    
    merdelmdrate <- lapply(mix, "[[", 4)
    merdelmdrate <- lapply(merdelmdrate, function(x) x[4:6])
    
    merdeltimes <-  lapply(mix, "[[", 5)
    merdeltimes <- data.frame(
      laps = names(merdeltimes),
      total_time = unlist(merdeltimes)
    )
    
    merdelClusters <- lapply(mix, "[[", 3)
    df_list <- list()
    
    # Iterate over each named vector in the list
    for (a in 1:6) {
      vec <- merdelClusters[[a]]
      df <- data.frame(
        laps = names(merdelClusters)[a],
        iter = seq_along(vec),
        clusters = vec
      )
      df_list[[a]] <- df
    }
    merdelClusters <- do.call(rbind, df_list)
    
    resultforpsm[[j]] <- merdelLabels
    ELBOresults[[j]] <- merdelELBO
    timeresults[[j]] <- merdeltimes
    mergesdeletes[[j]] <- merdelmdrate
    clusters[[j]] <- merdelClusters
    clustlabels[[j]] <- merdelclustLabels
  }
  
  resultforpsm <- bind_rows(resultforpsm, .id = "init_num")
  ELBOresults <- bind_rows(ELBOresults, .id = "init_num")
  timeresults <- bind_rows(timeresults, .id = "init_num")
  clusters <- bind_rows(clusters, .id = "init_num")
  
  all_result <- list(resultforpsm, ELBOresults, timeresults, clusters, mergesdeletes, truelabels, clustlabels)
  all_result
}

save(merdelsim1Cat, file="/home/jr951/VariationalMixtures/merdelsim1Cat.RData")

