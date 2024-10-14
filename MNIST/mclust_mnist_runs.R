#MCA mnist with mclust (EM clustering)
load("mnist_binnew.RData")
set.seed(30033)
library(mclust)

dataMatrix <- data.matrix(mnist_binnew)
tmp <- apply(dataMatrix, 2, as.factor)
mcaResult2 <- MCA(tmp, ncp = 20, graph=FALSE)
mcaToPlot2 <- as.data.frame(mcaResult2$ind$coord)

start.time <- Sys.time()
mod1 <- Mclust(mcaToPlot2[,1:5], G = 20)
end.time <- Sys.time()

start.time2 <- Sys.time()
mod2 <- Mclust(mcaToPlot2[,1:5], G = 35)
end.time2 <- Sys.time()

time_tot <- difftime(end.time, start.time)
time_tot2 <- difftime(end.time2, start.time2)

mnistmclust <- list(mod1, mod2, time_tot, time_tot2)

#Use this to then find misclassification rates