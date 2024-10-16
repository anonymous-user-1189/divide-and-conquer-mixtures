###
#Analysis of results for 'number of shards' simulation
#Only second simulation is shown: replace C with B for first simulation
###
globalmergeSimC <-  bind_rows(globalmergeSimC, .id = "data_number")

#Means/confidence intervals
test <- globalmergeSimC %>%
  group_by(model) %>%
  dplyr::summarize(meanARI = mean(ARI_8, na.rm=TRUE), sd_ARI = sd(ARI_8, na.rm = TRUE),     lower_CI = meanARI - 1.96 * sd_ARI * (1/sqrt(5)),
                   upper_CI = meanARI + 1.96 * sd_ARI * (1/sqrt(5))) 

#Medians/quantiles
globalmergeSimC %>%
  group_by(model) %>%
  dplyr::reframe(x = quantile(clusters_10, c(0.25, 0.5, 0.75), na.rm = TRUE), q = c(0.25, 0.5, 0.75))

#Perform Kruskal-Wallis test to compare number of shards
test <- globalmergeSimC %>%
  pivot_longer(
    cols = starts_with("ARI_"),      
    names_to = "shards",            #Shards column
    names_prefix = "ARI_",          
    values_to = "ARI"                # Name new ARI column
  ) %>%  mutate(
    shards = factor(shards, levels = c("4", "8")) 
  )
test <- test[!is.na(test$ARI),]

kruskal.test(ARI ~ shards, data = test) #Kruskal Wallis Test

#Plot in supplement (Figure 4)
ggplot(test, aes(x=shards, y=ARI, fill = shards)) + 
  geom_boxplot(notch=TRUE) + theme_bw() + ggtitle("N = 200,000")

ggsave("number_shards_C.png", width = 3, height = 3)

