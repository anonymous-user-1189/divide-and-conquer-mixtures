#####
#MCA transformation and kmeans on the MNIST data
#####
load("mnist_binnew.RData")
library(FactoMineR)
dataMatrix <- data.matrix(mnist_binnew)
tmp <- apply(dataMatrix, 2, as.factor)
mcaResult <- MCA(tmp, graph=FALSE)
plot(mcaResult$ind$coord[,1], mcaResult$ind$coord[,2], type = "p")

mcaToPlot        <- as.data.frame(mcaResult$ind$coord)
mcaToPlot$labels <- as.factor(mnist_labels$V1)

#More dimensions?

mcaResult2 <- MCA(tmp, ncp = 20, graph=FALSE)
mcaToPlot2 <- as.data.frame(mcaResult2$ind$coord)

library(factoextra)
get_eigenvalue(mcaResult2)
fviz_screeplot(mcaResult2, addlabels = TRUE, ylim = c(0, 45)) #visualise how much variance is explained

#20% variance explained by 5 dims, 30% by 10 dims, 40% by 17 dims

set.seed(153)
km.out1 <- kmeans(mcaToPlot[,1:5], iter.max = 20, centers = 20, nstart = 20)
km.out2 <- kmeans(mcaToPlot[,1:5], iter.max = 20, centers = 35, nstart = 20)
km.out3 <- kmeans(mcaToPlot2[,1:10], iter.max = 20, centers = 20, nstart = 20)
km.out4 <- kmeans(mcaToPlot2[,1:10], iter.max = 20, centers = 35, nstart = 20)
km.out5 <- kmeans(mcaToPlot2[,1:17], iter.max = 20, centers = 20, nstart = 20)
km.out6 <- kmeans(mcaToPlot2[,1:17], iter.max = 20, centers = 35, nstart = 20)





