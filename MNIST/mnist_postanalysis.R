#####
#Analysis of results from MNIST runs
#####

#####
#Heatmaps
#####
#Have all rows now in the same order as original dataset + labels
mnist_bin_results_analyse <- lapply(mnist_bin_results[[2]], function(x){
  x$shuffle1 <- x$shuffle1[order(x$row_shuffle)]
  x$shuffle2 <- x$shuffle2[order(x$row_shuffle)]
  x
})


annotationColor <- list(
  trueDigits = c('0' = "cyan", "1" = "red", "2" = "#487391", "3" = "violet", "4" = "#c4e2b7", "5" = "blue", "6" = "#FF7F00", "7" = "green","8" = '#cbbfec', "9" = "yellow")
)

currentAnnotationRow <- data.frame(
  #outcome = dataMatrix[[3]],
  #Cluster = factor(mnist_bin_results_analyse[[4]]$shuffle2[1:2000]),
  trueDigits = factor(mnist_labels[1:2000,])
)
rownames(currentAnnotationRow) <- rownames(mnist_binnew[1:2000,])

as.ggplot(pheatmap(mnist_bin_2000[order(mnist_labels[1:2000,]),], show_rownames = FALSE, show_colnames = FALSE, cluster_cols = FALSE, cluster_rows = FALSE, color = colorRampPalette(colors = c("white", "black"))(2), 
         annotation_row = currentAnnotationRow, annotation_colors = annotationColor))

ggsave("mnist_sorted.png", width = 13.5, height = 8) #sorted by true digit

currentAnnotationRow <- data.frame(
  #outcome = dataMatrix[[3]],
  Cluster = factor(mnist_bin_results_analyse[[4]]$shuffle2[1:2000]),
  trueDigits = factor(mnist_labels[1:2000,])
)
rownames(currentAnnotationRow) <- rownames(mnist_binnew[1:2000,])

as.ggplot(pheatmap(mnist_bin_2000[order(mnist_bin_results_analyse[[4]]$shuffle2[1:2000], mnist_labels[1:2000,]),], show_rownames = FALSE, show_colnames = FALSE, cluster_cols = FALSE, cluster_rows = FALSE, color = colorRampPalette(colors = c("white", "black"))(2), 
                   annotation_row = currentAnnotationRow, annotation_colors = annotationColor))

ggsave("mnist_sorted.png", width = 13.5, height = 8) #Sorted by cluster

#####
#Looking at some of the 0's and 1's, visualising in Supplement as heatmaps
#####
mnist_labels_10000 <- mnist_labels[1:10000,]

mnist_0and1 <- which(mnist_labels_10000 == 1 | mnist_labels_10000 == 0)
mnist_bin_0and1 <- mnist_binnew[mnist_0and1,]

annotationColor <- list(
  trueDigits = c('0' = "cyan", "1" = "red")
)

currentAnnotationRow <- data.frame(
  #outcome = dataMatrix[[3]],
  Cluster = factor(mnist_bin_results_analyse[[4]]$shuffle2[mnist_0and1]),
  trueDigits = factor(mnist_labels_10000[mnist_0and1])
)
rownames(currentAnnotationRow) <- rownames(mnist_bin_0and1)

as.ggplot(pheatmap(mnist_bin_0and1[order(mnist_labels_10000[mnist_0and1], mnist_bin_results_analyse[[4]]$shuffle2[mnist_0and1]), ], show_rownames = FALSE, show_colnames = FALSE, cluster_cols = FALSE, cluster_rows = FALSE, color = colorRampPalette(colors = c("white", "black"))(2), 
                   annotation_row = currentAnnotationRow, annotation_colors = annotationColor))
ggsave("mnist_0and1.png", width = 10, height = 6.5)

as.ggplot(pheatmap(mnist_bin_0and1[order(mnist_labels_10000[mnist_0and1], mnist_bin_results_analyse[[4]]$shuffle2[mnist_0and1]), sample(1:176)], show_rownames = FALSE, show_colnames = FALSE, cluster_cols = FALSE, cluster_rows = FALSE, color = colorRampPalette(colors = c("white", "black"))(2), 
                   annotation_row = currentAnnotationRow, annotation_colors = annotationColor))

ggsave("mnist_0and1shuffle.png", width = 10, height = 6.5)

######
#heatmap showing the correspondence between clusters and true numbers in the clustering model with the
#best ELBO for the MNIST data (in the Supplement)
#####
MNISTPercentHeat <- matrix(0, nrow = 10, ncol = 27)
numbersandclusters <- as.data.frame(cbind(mnist_labels$V1, mnist_bin_results_analyse[[4]]$shuffle2))
all_numbers <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)

for (i in 1:10){
  for (j in 1:27){
    MNISTPercentHeat[i, j] <- sum(numbersandclusters$V1 == all_numbers[i] & numbersandclusters$V2 == j) / sum(numbersandclusters$V2 == j)
  }
}

rownames(MNISTPercentHeat) <- all_numbers
colnames(MNISTPercentHeat) <- paste0(seq(ncol(MNISTPercentHeat)))

pheatmap(MNISTPercentHeat, display_numbers = F, 
         cluster_rows = F, cluster_cols = T, fontsize_number = 10,
         color = colorRampPalette(colors = c("white", "black"))(198),
         show_rownames = T, show_colnames = T, angle_col = 0, treeheight_col = 0)

#####
#Generating images of numbers
#####

num_gen <- mnist_bin_results[[2]][[4]]$eps2[[1]][1,,] #4 or a 9 with a bit missing
num_gen <- mnist_bin_results[[2]][[4]]$eps2[[1]][2,,] #8 but weird
num_gen <- mnist_bin_results[[2]][[4]]$eps2[[1]][3,,] #3
num_gen <- mnist_bin_results[[2]][[4]]$eps2[[1]][4,,] #1 
num_gen <- mnist_bin_results[[2]][[4]]$eps2[[1]][5,,] #6 - this one?
num_gen <- mnist_bin_results[[2]][[4]]$eps2[[1]][6,,] #7 or a 9 missing a bit
num_gen <- mnist_bin_results[[2]][[4]]$eps2[[1]][7,,] #9 missing bits or maybe a 4
num_gen <- mnist_bin_results[[2]][[4]]$eps2[[1]][8,,] #2
num_gen <- mnist_bin_results[[2]][[4]]$eps2[[1]][9,,] #3 or 5
num_gen <- mnist_bin_results[[2]][[4]]$eps2[[1]][10,,] #1
num_gen <- mnist_bin_results[[2]][[4]]$eps2[[1]][11,,] #7
num_gen <- mnist_bin_results[[2]][[4]]$eps2[[1]][12,,] #2
num_gen <- mnist_bin_results[[2]][[4]]$eps2[[1]][13,,] #0
num_gen <- mnist_bin_results[[2]][[4]]$eps2[[1]][14,,] #0 - this one
num_gen <- mnist_bin_results[[2]][[4]]$eps2[[1]][15,,] #7
num_gen <- mnist_bin_results[[2]][[4]]$eps2[[1]][16,,] #5
num_gen <- mnist_bin_results[[2]][[4]]$eps2[[1]][17,,] #4 or 9
num_gen <- mnist_bin_results[[2]][[4]]$eps2[[1]][18,,] #6
num_gen <- mnist_bin_results[[2]][[4]]$eps2[[1]][19,,] #6
num_gen <- mnist_bin_results[[2]][[4]]$eps2[[1]][20,,] #8
num_gen <- mnist_bin_results[[2]][[4]]$eps2[[2]][1,,] #1
num_gen <- mnist_bin_results[[2]][[4]]$eps2[[2]][2,,] #3
num_gen <- mnist_bin_results[[2]][[4]]$eps2[[3]] #5
num_gen <- mnist_bin_results[[2]][[4]]$eps2[[4]][1,,] #6
num_gen <- mnist_bin_results[[2]][[4]]$eps2[[4]][2,,] #3??
num_gen <- mnist_bin_results[[2]][[4]]$eps2[[4]][3,,] #5??
num_gen <- mnist_bin_results[[2]][[4]]$eps2[[4]][4,,] #6 but sideways

num_gen2 <- round(apply(num_gen, 2, function(x) rbeta(1, x[1], x[2]))) #generate random vector from beta dists

num_gennew <- numeric(256) + 1 #make new vector for the image
num_gennew[as.vector(which(category_count2 > 100))] <- num_gen2 #insert random vector into all the cells without 0 in the pixel

num_gen2m <- matrix(num_gennew, nrow = 16, ncol = 16)
num_gen2m <- num_gen2m[, ncol(num_gen2m):1]  # reverse cols
par(mar = rep(0, 4)) 
image(num_gen2m, col = c("black", "white"))    

#Find original images from clusters
row_num <- sample(which(mnist_bin_results_analyse[[4]]$shuffle2 == 27 
                        #& mnist_labels$V1 == 9
                        ), 1) #take random obs from cluster
og_num <- matrix(as.numeric(mnist_bin[row_num,]), nrow = 16, ncol = 16)
og_num <- og_num[, ncol(og_num):1]  # reverse cols
par(mar = rep(0, 4)) 
mnist_labels$V1[row_num] #what is the true label of the image?
image(og_num, col = c("white", "black"))  



#####
#Misclassification rates
#####

#'best' model: already made heat map for this!
cluster_char <- apply(MNISTPercentHeat, 2, which.max) - 1
#Gives dominant number in each cluster

all_digits = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
misclass = c()

for (i in 1:10){
  dig = all_digits[i]
  
  mnist_dig = which(mnist_labels == dig)
  mnist_dig_clusts <- mnist_bin_results_analyse[[4]]$shuffle2[mnist_dig]
  
  which_dig_clusts <- as.vector(which(cluster_char == dig))
  
  total = 0
  for (k in 1:length(which_dig_clusts)){
    total = total + sum(mnist_dig_clusts == which_dig_clusts[k])
  }
  
  misclass[i] = total / length(mnist_dig)
  
}
misclass

#For other models (change mnist_bin_results_analyse[[4]]$shuffle2... to relevant model)
MNISTPercentHeat2 <- matrix(0, nrow = 10, ncol = length(unique(mnist_bin_results_analyse[[4]]$shuffle2)))
numbersandclusters2 <- as.data.frame(cbind(mnist_labels$V1, mnist_bin_results_analyse[[4]]$shuffle2))
all_numbers <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)

for (i in 1:10){
  for (j in 1:length(unique(mnist_bin_results_analyse[[4]]$shuffle2))){
    MNISTPercentHeat2[i, j] <- sum(numbersandclusters2$V1 == all_numbers[i] & numbersandclusters2$V2 == j) / sum(numbersandclusters2$V2 == j)
  }
}
cluster_char2 <- apply(MNISTPercentHeat2, 2, which.max) - 1
for (i in 1:10){
  dig = all_digits[i]
  
  mnist_dig = which(mnist_labels == dig)
  mnist_dig_clusts <- mnist_bin_results_analyse[[4]]$shuffle2[mnist_dig]
  
  which_dig_clusts <- as.vector(which(cluster_char2 == dig))
  
  total = 0
  for (k in 1:length(which_dig_clusts)){
    total = total + sum(mnist_dig_clusts == which_dig_clusts[k])
  }
  
  misclass[i] = total / length(mnist_dig)
  
}
misclass

