#####
#Generate categorical clustered data
#####
library(gtools)

#Example settings for parameters in function
n    <- 1000 # The overall population size
K    <- 4    # The number of clusters/subpopulations
w    <- c(0.1, 0.2, 0.3, 0.4) # The mixture weights (proportion of population in each cluster)
w    <- w/sum(w)   # This ensures that the elements of w sum to 1 
p    <- 20 # The number of variables 
Irrp <- 5 #The number of irrelevant variables
yout <- TRUE #Indicate whether a binary outcome is required
cat  <- 3
#Have that the proportion of each category in each cluster is random - use Dirichlet dist
#Total data dimension is D = p + Irrp

#Returns a list of the data matrix as well as the 'true' cluster labels for each observation
GenerateSampleDataCat <- function(n, K, w, p, Irrp, yout = FALSE, cat = 2, ycat = 2){
  # Variable called "clusterLabels" to store the cluster labels, 
  
  clusterLabels <- sample(1:K, n, replace = T, prob = w)
  
  clusterParameters <- array(0, dim = c(p, cat, K))
  for(i in 1:K){
    for(j in 1:p)
    {
      for(c in 1:cat){
        clusterParameters[, ,i] <- rdirichlet(p, c(1:cat))
      }
    }
  }
  
  if(Irrp != 0){
    unclusteredParameters <- rdirichlet(Irrp, c(1:cat))
  }
  
  # Use cluster labels and cluster parameters to generate sample data
  
  dataMatrix    <- matrix(nrow = n, ncol = p + Irrp)
  for(i in 1:n)
  {
    currentClusterLabel <- clusterLabels[i]  # Get the cluster label for the i-th person
    for(j in 1:p)
    {
      # Simulate data according to the parameters for the current cluster
      dataMatrix[i,j] <- sample(1:cat, 1, replace = TRUE, prob = clusterParameters[j, , currentClusterLabel])
    }
    if (Irrp != 0){
      for (j in 1:Irrp){
        #Simulate data according to random probability, no clustering structure
        dataMatrix[i,p + j] <- sample(1:cat, 1, replace = TRUE, prob = unclusteredParameters[j,])
      }
    }
  }
  
  if (yout == TRUE){
    y <- vector(mode = "numeric", length = n)
    yParams <- matrix(0, nrow = K, ncol = ycat) 
    
    prob1 <- runif(1, 0, 10)
    prob2 <- runif(1, 0, 10)
    yParams <- rdirichlet(K, c(1:ycat)) #generate random Dirichlet probability (extension of Beta)
    
    for (i in 1:n){
      currentClusterLabel <- clusterLabels[i]
      y[i] <- sample(1:ycat, 1, replace = TRUE, prob = yParams[currentClusterLabel,])
    }
    
  }
  
  #Create an "annotation row" to show the cluster label for each person
  annotationRow <- data.frame(
    Cluster = factor(clusterLabels)
  )
  #Create data frame for outcome for each person if needed
  if (yout == TRUE){
    y <- data.frame(
      outcome = factor(y)
    )
  }
  
  #We require the data, outcome and annotation row to have the same rownames:
  rownames(annotationRow) <- rownames(dataMatrix) <- paste0("Person", seq(1,n))
  if(yout == TRUE){
    rownames(y) <- rownames(dataMatrix)
  }
  
  if(yout == TRUE){
    return(list(dataMatrix, annotationRow, y))
  } else{
    return(list(dataMatrix, annotationRow))
  }
}
