#####################################################################
# Produces the plots and tables as seen in the paper 
# "Adaptive Selection of the Optimal Strategy to Improve Precision and Power in Randomized Trials"
#####################################################################
rm(list = ls())
#####################################################################
# Inputs for the required plots and tables
# 1. Sample size n = 40, 100, 500
n <- 500
# 2. Number of replications 
nReps <- 500
# 3. Number of folds in cross validation, V = 10 when n = 40, 100 and V = 5 otherwise
V <- 5
# 4. Specify whether we should include MARS 
incl.mars <- T
# 5. Verbose setting describes in detail the outputs produced 
verbose <- F
# 6. Flag for continuous outcome 
sim_flag <- T
if(sim_flag == TRUE){
  sim <- "contY"
} else {
  sim <- "binY"
}
# 7. Specify whether there is an effect or whether the effect is null 
effect <- T
# 8. Specify the Data Generating Process used for the simulated data in Sim_Functions.R 
if(sim=='contY'){
  expt_type <- c('noisy_linear_1','noisy_multicollinear_cand1', 'noisy_polynomial')
  null.value=0
} else{
  expt_type <- c('noisy_linear','noisy_multicollinear', 'noisy_polynomial')
  null.value=1
}
#####################################################################

#=====================================================
# Function to return the Mean Squared error for the output 
#=====================================================
get.MSE <- function(output){
  mean( (output$est - output$psi)^2 ,na.rm=TRUE )
}

#=====================================================
# Function that returns other metrics such as 
# cover: 95% confidence interval contained the truth?
# reject: null hypo of no effect rejected
# Bias: ave deviation between pt and truth
# Variance: variance of point estimates
# MSE
#=====================================================
get.metrics <- function(estimator){
 
  yay <- c( colMeans(estimator[,c('cover', 'reject', 'bias')], na.rm=T),
        var(as.numeric(unlist(estimator["est"])),na.rm=TRUE),
        get.MSE(estimator)
        )
  yay <- data.frame(t(yay))
  colnames(yay) <- c('cover','power','bias', 'var', 'mse')
  yay
}

#=====================================================
# Produce the Tables containing all the metrics for
# 1. Unadjusted Estimator 
# 2. Fixed Estimator 
# 3. Small Adaptive Prespecification 
# 4. Large Adapative Prespecification 

# We do the same for both the simple design as well as the stratified design
#=====================================================
YAY <- NULL
ests <- c('Unadjusted', 'Static', 'Small APS', 'Large APS')
STRATIFY <- c(F,T)
dgp <- c('Linear', 'Interactive', 'Polynomial')


#=====================================================
# Load the output files that are produced by running Main.R and compute the metrics 
#=====================================================
for(j in 1:length(expt_type)){
  for(k in 1:2){
  file.name <- paste( sim, paste0('Effect', effect),
                      paste0('N', n), paste0('V',V), paste0('mars', incl.mars),
                      paste0('nReps', nReps),paste0('stratify', STRATIFY[k]),
                      paste0('type', expt_type[j]), sep = "_")
  
  file.nameD <- paste( "OUTPUT/", file.name, paste('.RData'), sep = "_")
  print(paste0("Experiment file name is: ", file.nameD))
  load(file.nameD)
  
  # After loading the outputs, read the metrics 
  SIMPLE <- OUT.AP
  FANCY <- OUT 
  yay <- data.frame(rbind( 
            get.metrics(UNADJ), get.metrics(FORCE), 
            get.metrics(SIMPLE),
            get.metrics(FANCY)))
  yay <- cbind(expt=expt_type[j], stratify=STRATIFY[k], 
               ests, yay, var.ratio=yay[1,'var']/yay[,'var'], re=yay[,'mse']/yay[1,'mse'] )
  yay <- cbind(yay, savings=(1-yay$re))
  print(paste0("Unadjusted Psi", round(mean(UNADJ$psi),2)))
  YAY <- rbind(YAY, yay)
  
  # Create the data frame that stores all the metrics
  data <- data.frame(
    x=c(1:2000), 
    value1=c(UNADJ[["CI.lo"]],FORCE[["CI.lo"]],SIMPLE[["CI.lo"]],FANCY[["CI.lo"]]), 
    value2=c(UNADJ[["CI.hi"]],FORCE[["CI.hi"]],SIMPLE[["CI.hi"]],FANCY[["CI.hi"]]),
    ests=c(rep('Unadjusted',500), rep('Static',500), rep('Small APS', 500), rep('Large APS', 500))
  )
  
  psi <- mean(UNADJ$psi)
  # Plot
  ggplot(data) +
    geom_segment( aes(x=x, xend=x, y=value1, yend=value2), color="grey") +
    #geom_point( aes(x=x, y=value1), color=rgb(0.2,0.7,0.1,0.5), size=3 ) +
    #geom_point( aes(x=x, y=value2), color=rgb(0.7,0.2,0.1,0.5), size=3 ) +
    geom_point( aes(x=x, y=value1, color=factor(ests)), size=1 ) +
    geom_point( aes(x=x, y=value2, color=factor(ests)), size=1 ) +
    geom_hline(yintercept=null.value, linetype='dashed', color='black', size=0.5) +
    coord_flip()+
    theme_ipsum() +
    labs(title = dgp[j])+
    theme(
      legend.position = "none",
      #plot.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
  ggsave(filename = paste0("PLOTS/", file.name, '.eps'))
  rm(yay, data)
}
}

#=====================================================
# Create the tables for the metrics for inputs specified 
#=====================================================
library(xtable)
this.order <- c('ests', 'cover', 'power', 'mse', 'bias', 'var', 're')
print(paste0("Table for metrics: "))
YAY

# Generate table for latex version 
xtable(YAY[,this.order], digits=c(1, 1, rep(3, 6) ))

#=====================================================
# Print savings obtained while using APS compared to unadjusted estimator
#=====================================================
round( summary(YAY[YAY$ests=='Large APS', 're']), 3)
round( summary(YAY[YAY$ests=='Large APS', 'savings'])*100, 0)

round( summary(YAY[YAY$ests=='Small APS', 're']), 3)
round( summary(YAY[YAY$ests=='Small APS', 'savings'])*100, 0)


#=====================================================
# Find propoortion of times when different candidate algorithms where chosen
# Print the tables for Outcome and PScore
#=====================================================
xtable(WINNERQ)
xtable(WINNERG)