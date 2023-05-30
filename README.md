# Adaptive Selection of the Optimal Strategy to Improve Precision and Power in Randomized Trials
By Laura B. Balzer (laura.balzer@berkeley.edu), Erica Cai, 	Lucas Godoy Garraza, and Pracheta Amaranath

ArXiv: https://arxiv.org/abs/2210.17453

In `vignette.Rmd` we provide worked examples of using Adaptive Prespecification (APS) with targeted minimum loss-based estimation (TMLE) for empirical efficiency maximization in randomized trials. From a pre-specified set, APS data-adaptively selects the optimal combination of estimators of the *outcome regression* (i.e., conditional expectation of the outcome, given the randomized intervention and candidate covariates) and of the known *propensity score* (i.e., conditional probability of the intervention, given the candidate covariates) to minimize the cross-validated variance estimate. 

Key methods references include

- Balzer et al., [Adaptive pre-specification](https://pubmed.ncbi.nlm.nih.gov/27436797/) in randomized trials with and without pair-matching, *Statistics in Medicine*, 2016
- Balzer et al., [Two-Stage TMLE](https://pubmed.ncbi.nlm.nih.gov/34939083/) to reduce bias and improve efficiency in cluster randomized trials, *Biostatistics*, 2021
- Balzer et al., [Adaptive Selection](https://arxiv.org/abs/2210.17453)  of the Optimal Strategy to Improve Precision and Power in Randomized Trials, *arXiv*, 2022

Example applications include

- Havlir et al., [HIV Testing and Treatment](https://pubmed.ncbi.nlm.nih.gov/31314966/) with the Use of a Community Health Approach in Rural Africa, *NEJM*, 2019 with corresponding Statistical Analysis Plan [(SAP)](https://arxiv.org/abs/1808.03231)
- Kakende et al., [A mid-level health manager intervention](https://pubmed.ncbi.nlm.nih.gov/35908553/) to promote uptake of isoniazid preventive therapy among people with HIV in Uganda: a cluster randomised trial, *LancetHIV*, 2022 with corresponding [SAP](https://arxiv.org/abs/2111.10467)
- Hickey et al., [Effect of a one-time financial incentive](https://pubmed.ncbi.nlm.nih.gov/36342940/) on linkage to chronic hypertension care in Kenya and Uganda: A randomized controlled trial, *PLoSOne*, 2022  (corresponding SAP included in article's supplementary materials)

Please see `vignette.Rmd` for worked in examples with corresponding output in `vignette.html`.

To reproduce the simulation studies: (1) Adjust inputs in ``Sims_Main.R`` by specifying the sample size, type of outcome (continuous or binary), data generating process for the simulated data, presence of effect, stratification, number of replications, number of folds for cross validation, etc. (2) Run ``Rscript Sims_Main.R`` whose output will be stored in the ``OUTPUT`` folder. This output can summarized usings `Sims_MakePretty.R` to generate the tables and figures.

