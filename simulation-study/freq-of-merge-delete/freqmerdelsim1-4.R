#Frequency of merge/delete Simulation 1.4 - N=1000, P=100, 10 clusters - varsel, 75%

library(tidyverse)
library(mclust)
library(mcclust)
library(mcclust.ext)
source("mixtureMergeDeleteTimeVarSel.R")
source("GenerateSampleData.R")

set.seed(1233)

library(foreach)
library(doParallel)
library(doRNG)
registerDoParallel(20)

merdelsim1vs <- foreach(k = 1:20) %dorng% {
  generation <- GenerateSampleData(1000, 10, c(0.05, 0.05, 0.1, 0.1, 0.25, 0.05, 0.1, 0.1, 0.1, 0.1), 75, 25)
  data <- generation[[1]]
  truelabels <- generation[[2]]
  resultforpsm <- list()
  ELBOresults <- list()
  timeresults <- list()
  mergesdeletes <- list()
  clusters <- list()
  clustlabels <- list()
  selected_vars <- list()
  for (j in 1:10){
    mix <- mixturemodel_merdelTIMEVarSel(data, 20, 0.01, 2, 1500, 0.00000005, c(1, 2, 5, 10, 10000))
    
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
    
    merdelselect <- lapply(mix, "[[", 4)
    merdelselect <- lapply(merdelselect, "[[", 7)
    
    resultforpsm[[j]] <- merdelLabels
    ELBOresults[[j]] <- merdelELBO
    timeresults[[j]] <- merdeltimes
    mergesdeletes[[j]] <- merdelmdrate
    clusters[[j]] <- merdelClusters
    clustlabels[[j]] <- merdelclustLabels
    selected_vars[[j]] <- merdelselect
  }
  
  resultforpsm <- bind_rows(resultforpsm, .id = "init_num")
  ELBOresults <- bind_rows(ELBOresults, .id = "init_num")
  timeresults <- bind_rows(timeresults, .id = "init_num")
  clusters <- bind_rows(clusters, .id = "init_num")
  
  all_result <- list(resultforpsm, ELBOresults, timeresults, clusters, mergesdeletes, truelabels, clustlabels, selected_vars)
  all_result
}

save(merdelsim1vs, file="merdelsim1vs.RData")

