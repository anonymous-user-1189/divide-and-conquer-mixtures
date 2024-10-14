#Preprocessing MNIST data (AFTER Jupyter notebook)

mnist_bin <- read.csv("bin_mnist.csv", header = FALSE) #saved from Jupyter notebook
mnist_labels <- read.csv("labels_mnist.csv", header = FALSE) #saved from Jupyter notebook

category_count2 <- apply(mnist_bin, 2, function(col) sum(col != 0))
mnist_binnew <- mnist_bin[,category_count2 > 100]
#176 columns now, removing all pixels with less than 100 non-0 values

save(mnist_binnew, file = "mnist_binnew.RData") 
#This is used in all main MNIST simulations, apart from mnist_runs_testset.R

#Visualising some of the numbers
mnist_labels_2000 <- mnist_labels[1:2000,]
mnist_bin_2000 <-mnist_binnew[1:2000,]

library(pheatmap)
pheatmap(mnist_bin_2000[order(mnist_labels_2000),], show_rownames = FALSE, show_colnames = FALSE, cluster_cols = FALSE, cluster_rows = FALSE, color = colorRampPalette(colors = c("white", "black"))(4))

#####
#Preprocess MNIST test set (used in 'mnist_runs_testset.R')

mnist_test_bin <- read.csv("bin_mnist_test.csv", header = FALSE)
mnist_test_labels <- read.csv("labels_mnist_test.csv", header = FALSE)

category_count3 <- apply(mnist_test_bin, 2, function(col) sum(col != 0))
mnist_test_binnew <- mnist_test_bin[,category_count2 > 50]
#181 columns now. Reduce to 50 non-zero values becuse of smaller dataset

save(mnist_test_binnew, file = "mnist_test_binnew.RData")