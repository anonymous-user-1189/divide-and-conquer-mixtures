#Example of running MerDel
library(tidyverse)
library(mclust)
library(mcclust)
library(mcclust.ext)
source("VariationalMixtureMergeDelete.R")
source("GenerateSampleData.R")
source("VariationalMixtureMergeDeleteVarSel.R")
source("globalMerge.R")
source("globalMergeAlt.R")

#1) Generate sample data:

generation <- GenerateSampleData(1000, 10, c(0.1, 0.2, 0.1, 0.2, 0.05, 0.05, 0.05, 0.05, 0.1, 0.1), 80, 20)
data <- generation[[1]] #This is the data matrix
truelabels <- generation[[2]] #The 'true' simulated labels are saved here

#2) Run MerDel with no variable selection

#20 is the initialised number of clusters
#0.01 is the initialised alpha value (set as some value below 1)
#1500 is the maximum number of VICatMix iterations - usually does not go above 200/300 regardless
#0.00000005 is the convergence parameter - if values are within this difference for 3 iterations, alg ends
#5 represents the frequency of merge/delete moves. (5 VICatMix moves before 1 merge/delete move)
#If you don't want any merge/delete, you can write 'all' or just replace 5 with 100000

mix <- mixturemodel_merdel(data, 20, 0.01, 1500, 0.00000005, 5)

#Labels are saved in mix$model$labels
#All model parameters are saved in mix$model
#ELBO and clusters (throughout the algorithm) are also saved, 'secondELBO' is attempting to take into account extra EM moves in the merge/delete moves
#mix$prior and mix$nCat are used for the global merge

#eg. finding ARI
library(mclust)
adjustedRandIndex(mix$model$labels, truelabels$Cluster)


#3) Splitting data into shards, running global merge - sMerDel

sampledata <- GenerateSampleData(5000, 8, c(0.2, 0.2, 0.15, 0.15, 0.05, 0.05, 0.1, 0.1), 100, 0)
data <- sampledata[[1]] #Data matrix

row_shuffle <- sample(1:nrow(data)) #Shuffling the rows
label_shuffle <- sampledata[[2]]$Cluster[row_shuffle] #Shuffling 'true labels' in line with data matrix
data_split <-  data[row_shuffle, ] #Shuffling rows of data

#SHARDS
data_split <- list(data_split[1:1000,], data_split[1001:2000,], data_split[2001:3000,], data_split[3001:4000,], data_split[4001:5000,])

#For global merge, need all results of models in a list

library(foreach)
library(doParallel)
library(doRNG)
registerDoParallel(5)

shuffled_results <- foreach(i = 1:5) %dorng%{
  mix <- mixturemodel_merdel(data_split[[i]], 15, 0.01, 1000, 0.000005, 5)
  prior <- mix$prior #*Should* be the same for each mixture model
  mix$model #each item in shuffled_results is model for one of the shards
}

#Results saved in 'shuffled_results' should be a list of B models
#prior must be the same as used in each mixture model 

#Global Merge: Greedy search (see paper)

global1 <- globalMerge(shuffled_results, prior)
overallLabels1 <- unlist(global1$labels) #Final cluster labels

#Global Merge: Random search (see paper)

global2 <- globalMergeAlt(shuffled_results, prior)
overallLabels2 <- unlist(global2$labels) #Final cluster labels

#Variable selection:
#4) Run one model with variable selection
#Extra parameter '2' is a variable selection hyperparameter; also shouldn't need to adjust this much

mix2 <- mixturemodel_merdelvarsel(data, 20, 0.01, 2, 1500, 0.00000005, 5)

#Output is similar to before
adjustedRandIndex(mix2$model$labels, truelabels$Cluster)
#selected variables can be found via
mix2$model$c



