#####
#Example of post-processing and analysis of results from 'Frequency of Merge/Delete Moves' simulation
#'merdelsim1' and 'merdelsim1vs' used as examples (Simulations 1.1 and 1.4, without and with variable selection)
#These can be replaced with the relevant simulation
#####
library(tidyverse)
###
#Line graphs of time vs ELBO (supplementary figure 2)
###

merdelELBO <- lapply(merdelsim1, "[[", 2)
merdelELBO <- bind_rows(merdelELBO, .id = "data_no")
merdelELBO$laps <- factor(merdelELBO$laps, levels = c(0, 1, 2, 5, 10, 10000))
colnames(merdelELBO)[4] <- "log-ELBO"

condensed <- merdelELBO[merdelELBO$data_no == 1,] #Change for a different simulated dataset
#Graph gets too crowded if all datasets are on one graph!

#Getting vertical line for where it end
merdeltimesMax <- condensed %>%
  group_by(laps, data_no, init_num) %>%
  slice(which.max(Time)) %>% 
  group_by(laps, data_no) %>% #comment these last two lines out for just maxes
  summarise(mean_time=mean(Time), sd_time=sd(Time))
merdeltimesMax$laps <- factor(merdeltimesMax$laps, levels = c(0, 1, 2, 5, 10, 10000))

merdeltimesMax$lower.ci <- merdeltimesMax$mean_time - 1.96 * merdeltimesMax$sd_time / sqrt(10)
merdeltimesMax$upper.ci <- merdeltimesMax$mean_time + 1.96 * merdeltimesMax$sd_time/sqrt(10)


ggplot() +
  geom_line(data=condensed[condensed$Time != 0,], aes(x=Time, y=`log-ELBO`, group=interaction(laps, init_num), color = laps)) +
  theme_bw() + 
  geom_rect(data=merdeltimesMax, aes(xmin=lower.ci, xmax=upper.ci, ymin=-Inf, ymax=Inf, fill = laps), alpha = 0.25) + 
  geom_vline(data = merdeltimesMax, aes(xintercept = mean_time, group = laps, color = laps),  lty = "dashed")

###
#Number of clusters: plot example in Supplementary Figure 1
###

merdelclusters <- lapply(merdelsim1, "[[", 4)
merdelclusters <- bind_rows(merdelclusters, .id = "data_num")
merdelclusters$laps <- factor(merdelclusters$laps, levels = c(0, 1, 2, 5, 10, 10000))

#Plots how clusters evolve through time - plot not included
ggplot() +
  geom_line(data=merdelclusters, aes(x=iter, y=clusters, group=interaction(laps, init_num, data_num), color = laps)) +
  theme_bw()

#Find final number of clusters
merdelclustersMax <- merdelclusters %>% group_by(laps, init_num, data_num) %>%
  filter(iter == max(iter))
merdelclustersMax$laps = as.factor(merdelclustersMax$laps)

#Supplementary Figure 1
plotclusters1 <- ggplot(merdelclustersMax, aes(x = laps, y = clusters)) + geom_hline(yintercept = 5, lty = 'dotted') + geom_count() +
  scale_size_area(max_size = 3) +  # Optionally use scale_size_area for better size scaling
  labs(size = "count") +
  theme_classic() + theme(legend.position="none") 

merdelclustersMax %>% group_by(laps) %>% #comment these last two lines out for just maxes
  summarise(mean_clusters=mean(clusters), sd_clusters=sd(clusters)) %>%
  mutate(lower.ci = mean_clusters - 1.96 * sd_clusters / sqrt(200), upper.ci = mean_clusters + 1.96* sd_clusters/ sqrt(200))

###
#Comparing ARI (eg. Figure 1 in main text)
###

merdelARI <- lapply(merdelsim1, "[[", 1)
merdelARI<- bind_rows(merdelARI, .id = "data_no")
rownames(merdelARI) <- 1:nrow(merdelARI)
merdelARI$laps <- factor(merdelARI$laps, levels = c(0, 1, 2, 5, 10, 10000))

#Plot example (Figure 1 in main text)
plotARI1 <- ggplot(merdelARI, aes(x = laps, y = ARI, shape = laps, color = laps)) + geom_jitter(width = 0.1, height = 0, size = 0.4) + theme_classic() + theme(legend.position="none") 

#Quantiles, confidence intervals etc.
ARItable <- merdelARI %>%
  group_by(laps, data_no) %>%
  dplyr::reframe(x = quantile(ARI, c(0.25, 0.5, 0.75), na.rm = TRUE), q = c(0.25, 0.5, 0.75))
merdelARI %>% group_by(laps) %>% #comment these last two lines out for just maxes
  summarise(mean_ARI=mean(ARI), sd_ARI=sd(ARI)) %>%
  mutate(lower.ci = mean_ARI - 1.96 * sd_ARI / sqrt(200), upper.ci = mean_ARI + 1.96* sd_ARI/ sqrt(200))

###
#Comparing total time taken
###

merdelTime <- lapply(merdelsim1, "[[", 3)
merdelTime<- bind_rows(merdelTime, .id = "data_no")
merdelTime %>%
  group_by(laps) %>%
  dplyr::reframe(x = quantile(total_time, c(0.25, 0.5, 0.75), na.rm = TRUE), q = c(0.25, 0.5, 0.75))

merdelTime %>% group_by(laps) %>% #comment these last two lines out for just maxes
  summarise(mean_time=mean(total_time), sd_time=sd(total_time)) %>%
  mutate(lower.ci = mean_time - 1.96 * sd_time / sqrt(200), upper.ci = mean_time + 1.96* sd_time/ sqrt(200))

###
#Comparing final ELBO value (in tables in Supplement)
###

merdelELBO %>%
  group_by(laps, data_no, init_num) %>%
  slice(which.max(secondELBO)) %>% 
  group_by(laps) %>% #comment these last two lines out for just maxes
  summarise(mean_ELBO=mean(secondELBO), sd_ELBO=sd(secondELBO))

###
#VARIABLE SELECTION: Comparing F1 scores for Simulation 1.4 
###

#Find selected variables
merdelsel1vs <- lapply(merdelsim1vs, "[[", 8)

for (i in 1:20){
  currentlist <- merdelsel1vs[[i]]
  for (j in 1:10){
    currentlist2 <- currentlist[[j]]
    for (a in 1:6){
      currentlist2[[a]] <- data.frame(
        laps = names(currentlist2)[a],
        Rel = sum(currentlist2[[a]][1:75] > 0.5),
        Irrel = sum(currentlist2[[a]][76:100] <= 0.5)
      )
    }
    currentlist2 <- bind_rows(currentlist2)
    currentlist[[j]] <- currentlist2
  }
  currentlist <- bind_rows(currentlist, .id = "init_num")
  merdelsel1vs[[i]] <- currentlist
}

merdelsel1vs <- bind_rows(merdelsel1vs, .id = "data_num")
merdelsel1vs$F1 <- (2 * merdelsel1vs$Rel) / ((2 * merdelsel1vs$Rel) + (25 - merdelsel1vs$Irrel) + (75 - merdelsel1vs$Rel))
merdelsel1vs %>% group_by(laps) %>% dplyr::summarize(Mean_F1 = mean(F1, na.rm=TRUE), Mean_Rel = mean(Rel, na.rm=TRUE), Mean_Irrel = mean(Irrel, na.rm=TRUE))

