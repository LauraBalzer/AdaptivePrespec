# Adaptive Prespecification in TMLE
Code for "Adaptive Selection of the Optimal Strategy to Improve Precision and Power in Randomized Trials"
By Laura B. Balzer (laura.balzer@berkeley.edu), Erica Cai, 	Lucas Godoy Garraza, and Pracheta Amaranath

ArXiv: https://arxiv.org/abs/2210.17453

For demonstration using real data from ACTG Study 175, please see the `vignette.Rmd`. The vignette also includes references to other applications using APS with TMLE and the Statistical Analysis Plan (SAP) for those applications.

To reproduce the simulation studies as seen in the paper: (1) Adjust inputs in ``Main.R`` by specifying the sample size, type of outcome (continuous or binary), data generating process for the simulated data, presence of effect, stratification, number of replications, number of folds for cross validation. (2) Run ``Rscript Main.R``. This results in a table computing metrics for all estimators stored in the ``OUTPUT`` folder. 

