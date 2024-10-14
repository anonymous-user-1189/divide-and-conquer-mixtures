#MNIST with flexmix: EM algorithm on binary data
set.seed(373)
library(flexmix)
load("mnist_binnew.RData")

dataMatrix <- data.matrix(mnist_binnew)
rownames(dataMatrix) <- paste0("Number", seq(1,60000))

start.time <- Sys.time()
flex <- initFlexmix(dataMatrix ~ 1, k = 20, model = FLXMCmvbinary(), control = list(minprior = 0), nrep = 10)
end.time <- Sys.time()
time_tot <- difftime(end.time, start.time)

start.time <- Sys.time()
flex2 <- initFlexmix(dataMatrix ~ 1, k = 35, model = FLXMCmvbinary(), control = list(minprior = 0), nrep = 10)
end.time <- Sys.time()
time_tot3 <- difftime(end.time, start.time)

mnistflex <- list(flex, flex2, time_tot, time_tot3)

#Use this to then find misclassification rates