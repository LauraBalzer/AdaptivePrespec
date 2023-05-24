# AdaptivePrespec
Code for "Adaptive Selection of the Optimal Strategy to Improve Precision and Power in Randomized Trials"
By Laura B. Balzer (laura.balzer@berkeley.edu), Erica Cai, 	Lucas Godoy Garraza, and Pracheta Amaranath

ArXiv: https://arxiv.org/abs/2210.17453

To reproduce the experiments as seen in the paper: 
1. Adjust inputs in ``Main.R`` by specifying the sample size, type of outcome (continuous or binary), data generating process for the simulated data, presence of effect, stratification, number of replications, number of folds for cross validation.
2. Run ``Rscript Main.R``. This results in a table computing metrics for all estimators stored in the ``OUTPUT`` folder. 
3. Run ``Rscript MakePretty.R`` to reproduce the plots and tables as seen in the paper. 


