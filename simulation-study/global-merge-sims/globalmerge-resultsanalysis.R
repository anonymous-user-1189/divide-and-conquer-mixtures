#Global merge simulation - analysis of results - example with results from 'globalmergeSim7.R'
library(tidyverse)
globalmergeSim7 <- bind_rows(globalmergeSim7, .id = "data_number")
globalmergeSim7$merge_total = globalmergeSim7$mode_time + globalmergeSim7$merge_time

boxplot(ARI ~ model, data = globalmergeSim7)

globalmergeSim7 %>%
  group_by(model) %>%
  dplyr::reframe(x = quantile(merge_time, c(0.25, 0.5, 0.75), na.rm = TRUE), q = c(0.25, 0.5, 0.75))






