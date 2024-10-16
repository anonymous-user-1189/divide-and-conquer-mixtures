#Code to plot results from simulation comparing merge criteria

library(tidyverse)
library(patchwork)

###
#Time
###
corrdivrand2 <- lapply(all_sims1a, "[[", 3)
corrdivrand2 <- bind_rows(corrdivrand2, .id = "data_num")
corrdivrand2 <- corrdivrand2[(corrdivrand2$laps == 1 | corrdivrand2$laps == 5),]

corrdivrand22 <- lapply(all_sims1b, "[[", 3)
corrdivrand22 <- bind_rows(corrdivrand22, .id = "data_num")
corrdivrand22 <- corrdivrand22[(corrdivrand22$laps == 1 | corrdivrand22$laps == 5),]

corrdivrand2 <- rbind(corrdivrand2, corrdivrand22)

ggplot(corrdivrand2, aes(x=model, y = log(total_time))) + geom_boxplot()

cdr_time2 <- ggplot(corrdivrand2, aes(x = laps, y = log(total_time), fill = model)) + stat_boxplot(geom = "errorbar", width = 0.5, position = position_dodge(width = 0.76)) + geom_boxplot(inherit.aes = TRUE, fatten = 1) + theme_bw()

#ggsave("corrdivrand2_time.png", width = 3.5, height = 3)

corrdivrand3 <- lapply(all_sims2a, "[[", 3)
corrdivrand3 <- bind_rows(corrdivrand3, .id = "data_num")
corrdivrand3 <- corrdivrand3[(corrdivrand3$laps == 1 | corrdivrand3$laps == 5),]

corrdivrand33 <- lapply(all_sims2b, "[[", 3)
corrdivrand33 <- bind_rows(corrdivrand33, .id = "data_num")
corrdivrand33 <- corrdivrand33[(corrdivrand33$laps == 1 | corrdivrand33$laps == 5),]

corrdivrand3 <- rbind(corrdivrand3, corrdivrand33)

cdr_time3 <- ggplot(corrdivrand3, aes(x = laps, y = log(total_time), fill = model)) + stat_boxplot(geom = "errorbar", width = 0.5, position = position_dodge(width = 0.76)) + geom_boxplot(inherit.aes = TRUE, fatten = 1) + theme_bw()

cdr_time2 + cdr_time3 + plot_annotation(tag_levels = 'A') #Supp figure 5(a)

###
#ARI plots
###

corrdivrand2 <- lapply(all_sims1a, "[[", 1)
corrdivrand2 <- bind_rows(corrdivrand2, .id = "data_num")
corrdivrand2 <- corrdivrand2[(corrdivrand2$laps == 1 | corrdivrand2$laps == 5),]

corrdivrand22 <- lapply(all_sims1b, "[[", 1)
corrdivrand22 <- bind_rows(corrdivrand22, .id = "data_num")
corrdivrand22 <- corrdivrand22[(corrdivrand22$laps == 1 | corrdivrand22$laps == 5),]

corrdivrand2 <- rbind(corrdivrand2, corrdivrand22)

cdr_ARI2 <- ggplot(corrdivrand2, aes(x = laps, y = ARI, fill = model)) + stat_boxplot(geom = "errorbar", width = 0.5, position = position_dodge(width = 0.76)) + geom_boxplot(inherit.aes = TRUE, fatten = 1) + theme_bw()

#ggsave("corrdivrand2_time.png", width = 3.5, height = 3)

corrdivrand3 <- lapply(all_sims2a, "[[", 1)
corrdivrand3 <- bind_rows(corrdivrand3, .id = "data_num")
corrdivrand3 <- corrdivrand3[(corrdivrand3$laps == 1 | corrdivrand3$laps == 5),]

corrdivrand33 <- lapply(all_sims2b, "[[", 1)
corrdivrand33 <- bind_rows(corrdivrand33, .id = "data_num")
corrdivrand33 <- corrdivrand33[(corrdivrand33$laps == 1 | corrdivrand33$laps == 5),]

corrdivrand3 <- rbind(corrdivrand3, corrdivrand33)

cdr_ARI3 <- ggplot(corrdivrand3, aes(x = laps, y = ARI, fill = model)) + stat_boxplot(geom = "errorbar", width = 0.5, position = position_dodge(width = 0.76)) + geom_boxplot(inherit.aes = TRUE, fatten = 1) + theme_bw()

cdr_ARI2 + cdr_ARI3 + plot_annotation(tag_levels = 'A') #Supp Figure 5(b)

###
#Compare merge rates 
###
corrdivrand2 <- lapply(all_sims1a, "[[", 2)
corrdivrand2 <- bind_rows(corrdivrand2, .id = "data_num")
corrdivrand2 <- corrdivrand2[(corrdivrand2$laps == 1 | corrdivrand2$laps == 5),]

corrdivrand22 <- lapply(all_sims1b, "[[", 2)
corrdivrand22 <- bind_rows(corrdivrand22, .id = "data_num")
corrdivrand22 <- corrdivrand22[(corrdivrand22$laps == 1 | corrdivrand22$laps == 5),]

corrdivrand2 <- rbind(corrdivrand2, corrdivrand22)

cdr_mergerate2 <- ggplot(corrdivrand2, aes(x = laps, y = merge_rate, fill = model)) + stat_boxplot(geom = "errorbar", width = 0.5, position = position_dodge(width = 0.76)) + geom_boxplot(inherit.aes = TRUE, fatten = 1) + theme_bw()

#ggsave("corrdivrand2_time.png", width = 3.5, height = 3)

corrdivrand3 <- lapply(all_sims2a, "[[",2)
corrdivrand3 <- bind_rows(corrdivrand3, .id = "data_num")
corrdivrand3 <- corrdivrand3[(corrdivrand3$laps == 1 | corrdivrand3$laps == 5),]

corrdivrand33 <- lapply(all_sims2b, "[[", 2)
corrdivrand33 <- bind_rows(corrdivrand33, .id = "data_num")
corrdivrand33 <- corrdivrand33[(corrdivrand33$laps == 1 | corrdivrand33$laps == 5),]

corrdivrand3 <- rbind(corrdivrand3, corrdivrand33)

cdr_mergerate3 <- ggplot(corrdivrand3, aes(x = laps, y = merge_rate, fill = model)) + stat_boxplot(geom = "errorbar", width = 0.5, position = position_dodge(width = 0.76)) + geom_boxplot(inherit.aes = TRUE, fatten = 1) + theme_bw()

cdr_mergerate2 + cdr_mergerate3 + plot_annotation(tag_levels = 'A') # Supp Figure 5(c)

