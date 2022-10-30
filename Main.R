#####################################################################
# Function to run Adaptive Prespecification for simulated data 
# This code reproduces all the experiments that are run in the paper 
# "Adaptive Selection of the Optimal Strategy to Improve Precision and Power in Randomized Trials"

# input: 
#   main_dir: Set working directory as the current folder to import all source files
#   outcome_flag: If true, sets the outcome as continuous outcome. If false, outcome is binary 
#   expt_type: "noisy_linear_1", "noisy_multicollinear_cand1", "noisy_polynomial" DGP used to generate synthetic data
#   effect: If False, generates outcome without an effect 
#   stratify: If True, applies stratified randomization
#   n: Sample size 
#   nReps: Number of replications 
#   V: number of folds used in cross validation 
#   incl.mars: If True, uses MARS as a candidate algorithm
#   verbose: If True, prints output in greater detail 

# output: 
#   computes metrics like MSE, relative efficiency, coverage etc for all four estimators 
#   1. Unadjusted
#   2. Fixed adjustment
#   3. Simple Adaptive Prespecification
#   4. Fancy Adaptive Prespecification 
#   and stores the resulting metrics for each estimator into a dataframe. 
#####################################################################
rm(list=ls())

#=====================================================
# Import libaries
#=====================================================
set.seed(1)
library('glmnet')
library('earth')
library(configr)

#=====================================================
# Set working directory
#=====================================================
main_dir <- "~/AdaptivePrespec"     # To be updated by user
setwd(main_dir)

#=====================================================
# Import all the necessary functions
#=====================================================
source("Sim_Functions.R") # New simulation data and functions for different synthetic data types
source('Stage2_Functions_Meta.R')
source('TMLE_Functions_Meta.R')
source('Adapt_Functions_Meta.R')

#=====================================================
# Input arguments to run experiments 
#=====================================================
outcome_flag <- TRUE

if(outcome_flag == FALSE){
  sim <- 'binY'
  goal <- 'aRR'
  # alt.smaller <- T
} else if(outcome_flag == TRUE){
  sim <- 'contY'
  goal <- 'RD'
  # alt.smaller <- F
}
expt_type <- "noisy_linear_1"
effect <- TRUE
stratify <- FALSE
n <- 500
nReps <- 500
V <- 5
incl.mars <- TRUE
verbose <- FALSE

#=====================================================
# Specify filename where we store the output
#=====================================================
file.name <- paste( "OUTPUT/", sim, paste0('Effect', effect), paste0('N', n), paste0('V',V), paste0('mars', incl.mars),
                    paste0('nReps', nReps),paste0('stratify', stratify), 
                    paste0('type', expt_type), paste('.RData'), sep = "_")
print(paste0("Experiment file name is: ", file.name))


#=====================================================
# Specify all possible candidates for adjustment 
#=====================================================

all.cand <- c('W1','W2','W3','W4','W5')

# Estimator 1: Simple Adaptive Prespecification 
# This estimator will automatically consider GLMs with a main term for one element of all.cand (and nothing=U)
AP.simple <- get.cand.adj(all.cand=all.cand, cand.Qform.fancy=NULL, cand.gform.fancy=NULL)

# Estimator 2: Fancy Adaptive Prespecification
# For this estimator, we need to specify candidate algorithms for adjusting for multiple covariates
# By default, we specify "glm", "stepwise", "step.interaction", "lasso" and optionally "mars"
if(incl.mars){
  cand.Qform.fancy <- c('glm', 'stepwise','step.interaction','lasso', 'mars') 
}else{
  cand.Qform.fancy <- c('glm', 'stepwise','step.interaction','lasso') 
}
cand.gform.fancy <- c('glm', 'stepwise','lasso')
AP.fancy <- get.cand.adj(all.cand=all.cand, cand.Qform.fancy=cand.Qform.fancy, 
                         cand.gform.fancy=cand.gform.fancy)

#=====================================================
# Defaults; for debugging purposes
#=====================================================
if(F){
  QAdj <- gAdj <- NULL
  Qform <- gform <- 'glm'
  one.sided <- F; sig.level=0.05; 
  scale_value <- 1; scale_value_min <- 0
}

#=====================================================
# Specify the metrics that we compute for all possible estimators 
#   "est": Estimate 
#   "psi": Effect
#   "CI.lo": lower bound for confidence interval (95%)
#   "CI.hi": upper bound for confidence interval (95%)
#   "se": standard error 
#   "bias": bias
#   "cover": coverage for 95% confidence interval 
#   "reject": whether the null hypothesis is rejected or not
#=====================================================
these <- c('Txt.est', 'Con.est', 'psi','est', 'CI.lo', 'CI.hi', 'se', 'bias', 'cover','reject')
# Output the results into a data.frame OUT
OUT <- data.frame(matrix(NA, nrow=nReps, ncol=2+length(these)) )
colnames(OUT) <- c('psi1','psi0', these)

SELECT <- data.frame(matrix(NA, nrow=nReps, ncol=4))
colnames(SELECT) <- c('QAdj', 'Qform', 'gAdj', 'gform')

UNADJ <- FORCE <-  OUT.AP <- OUT; SELECT.AP <- SELECT

#=====================================================
# For nReps replications, generate results for all estimators
#=====================================================
for(k in 1:nReps){
  full <- generate.data.wrapper(n=n, effect=effect, sim=sim, stratify=stratify, 
                                expt_type = expt_type, verbose=verbose)
  psi1 <- mean(full$Y1)
  psi0 <- mean(full$Y0)
  psi <- ifelse(goal=='aRR', psi1/psi0, psi1-psi0)
  data.input <- subset(full, select=-c(Y1,Y0))
  
  # unadjusted
  unadj <- Stage2(goal=goal, data.input=data.input, 
                  do.data.adapt =F, do.cv.variance=F, one.sided = F, # alt.smaller=alt.smaller,
                  verbose=verbose, psi=psi)
  # unadj
  UNADJ[k,] <- c(psi1, psi0, unadj[,these] )
  
  # force adjusment for W1 
  force <- Stage2(goal=goal, data.input=data.input, QAdj='W1',
                  do.data.adapt =F, do.cv.variance=F, one.sided = F, # alt.smaller=alt.smaller,
                  verbose=F, psi=psi)
  # 
  FORCE[k,] <- c(psi1, psi0, force[,these] )
  
  # simple adaptive prespec
  simp <- Stage2(goal=goal, data.input=data.input,
                 do.data.adapt=T, do.cv.variance=F, V=V, one.sided = F, #alt.smaller=alt.smaller,
                 cand.QAdj=AP.simple$cand.QAdj, cand.Qform=AP.simple$cand.Qform,
                 cand.gAdj=AP.simple$cand.gAdj, cand.gform=AP.simple$cand.gform,
                 verbose=verbose, psi=psi)
  # simp
  OUT.AP[k,] <- c(psi1, psi0, simp[,these] )
  SELECT.AP[k,] <- simp[,c('QAdj', 'Qform', 'gAdj', 'gform')]
  
  # fancy adaptive prespec
  fancy <- Stage2(goal=goal, data.input=data.input,
                  do.data.adapt=T, do.cv.variance=F, V=V, one.sided = F, # alt.smaller=alt.smaller,
                  cand.QAdj=AP.fancy$cand.QAdj, cand.Qform=AP.fancy$cand.Qform,
                  cand.gAdj=AP.fancy$cand.gAdj, cand.gform=AP.fancy$cand.gform,
                  verbose=verbose, psi=psi)
  
  OUT[k,] <- c(psi1, psi0, fancy[,these] )
  SELECT[k,] <- fancy[,c('QAdj', 'Qform', 'gAdj', 'gform')]
}

#=====================================================
# Print results summarizing all nReps replications
#=====================================================
print(paste0("File name is: ",file.name))
colMeans(UNADJ, na.rm=T)
colMeans(FORCE, na.rm=T)
colMeans(OUT.AP, na.rm=T)
colMeans(OUT, na.rm=T)

#=====================================================
# Compute additional metrics
#   get.MSE: Computes Mean Squared error 
#   get.RE: Computes the relative efficiency which is the ratio of the MSE of an estimator to the unadjusted estimate
#=====================================================
get.MSE <- function(output){
  mean( (output$est - output$psi)^2  ) 
}
get.RE <- function(UNADJ, FORCE, OUT.AP, OUT){
  data.frame( unadj= get.MSE(UNADJ)/get.MSE(UNADJ),
              force= get.MSE(FORCE)/get.MSE(UNADJ),
              simple= get.MSE(OUT.AP)/get.MSE(UNADJ),
              fancy= get.MSE(OUT)/get.MSE(UNADJ)
              )
}
get.RE(UNADJ, FORCE, OUT.AP, OUT)

#=====================================================
# Compute Risk Ratio for binary outcome
#=====================================================
ifelse(goal=='aRR', sqrt(var(log(OUT$est), na.rm=T)), sqrt(var(OUT$est, na.rm=T) ))
#table(SELECT.AP$QAdj)
#table(SELECT.AP$gAdj)

table(SELECT$QAdj)
table(SELECT$Qform)
table(SELECT$gAdj)
table(SELECT$gform)

#=====================================================
# Save the resulting data.frames to the output file
#=====================================================
save(UNADJ, FORCE, OUT.AP, OUT, SELECT.AP, SELECT, 
     AP.simple, AP.fancy,
     gen.data.contY, gen.data.binY,
     file=file.name)
