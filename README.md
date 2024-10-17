# divide-and-conquer-mixtures
Code to reproduce results from "Divide and Conquer Variational Mixture Models"

The folder model-code contains all R/cpp functions needed to implement MerDel and sMerDel, and code for generating sample data as detailed in Supplementary Section 8.1.
In general, all code implementing model on data (simulated and real) will 'source' all R scripts needed from the 'model-code' folder in order to run the simulation.

'example_code.R' contains basic examples of how to implement MerDel and sMerDel.

All hyperparameters etc. for running sMerDel for the EHR data are detailed in the paper, but note the code is not available as the data is not available publicly. 

## Dependencies

**Dependencies for Model**: Rcpp, RcppArmadillo, klaR, matrixStats, abind

**Additional Dependencies for R Simulations + Results Analysis**: tidyverse, pheatmap, doParallel, doRNG, foreach, parallel, mclust, flexmix, FactoMineR, factoextra, ggplotify

## Simulation study folder
-**freq-of-merge-delete**: 'freqmerdelsim1-x.R' refers to Simulation 1.x as described in Supplementary Section 8.3. 'freqmerdelsim-postprocess.R' demonstrates how to make plots and find confidence intervals etc. in the paper.

-**global-merge-sims**: 'globalmergeSim4.R', 'globalmergeSim5.R' and 'globalmergeSim7.R' refer to simulations with N=10000, N=20000, N=50000 respectively. 'globalmerge-resultsanalysis.R' has code for finding median/quantiles of metrics.

-**number-of-shards**: 'globalmergeSimB.R', 'globalmergeSimC.R' are simulations for N=100000, N=200000 respectively. 'noshards-analysis.R' has code for analysing results and creating plots.

-**comparing-merge-criteria**: 'comparingmergesim1a.R, comparingmergesim1b.R' are two halves of the same simulation, as are 'comparingmergesim2a.R, comparingmergesim2b.R'. Demonstrations for making plots etc. are in 'comparingmerge_plots.R'

-**categorical-sims**: No scripts for plotting are here - refer to the equivalent simulations for binary data to create the plots. 'freqmerdelsim-cat-x.R' are simulations for the frequency of merge/delete moves, 'globalmergeSim1Cat.R' and 'globalmergeSim2Cat.R' are simulations equivalent to 'global-merge-sims'.

-**global-merge-varsel**: No scripts for plotting are here - refer to the equivalent simulations for binary data to create the plots.

## MNIST folder
-Code to format MNIST data from Tensorflow into the data matrix used for analysis is in 'MNIST_preprocess1.ipynb' and 'MNIST_preprocess2.R'

-Main implementation of model: 'mnist_runs_main.R'

-Implementing model with test set: 'mnist_runs_testset.R'

-Comparator methods (for supplement): 'mnist_runs_merdelsupp.R', 'mclust_mnist_runs.R', 'kmeans_mnist_runs.R', 'flexmix_mnist_runs.R'

-Analysis of results and creating plots: 'mnist_postanalysis.R'




