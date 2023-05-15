#####################################################################
# Produces the plots and tables as seen in the paper 
# "Adaptive Selection of the Optimal Strategy to Improve Precision and Power in Randomized Trials"
#####################################################################
rm(list = ls())

# Import libraries 
library(ggplot2)
library(hrbrthemes)

#####################################################################
# Inputs to the file 
# 1. n: Sample size 
n <- 500
# 2. Number of replications 
nReps <- 5000
# 3. Number of folds in cross validation, V = 10 when n = 40, 100 and V = 5 otherwise
V <- 5
# 4. Specify whether we should include MARS 
incl.mars <- T
# 5. Verbose setting describes in detail the outputs produced 
verbose <- F
# 6. Flag for continuous outcome 
sim_flag <- F
if(sim_flag == TRUE){
  sim <- "contY"
} else {
  sim <- "binY"
}
# 7. Specify whether there is an effect or whether the effect is null 
effect <- T
# 8. Specify the Data Generating Process used for the simulated data in Sim_Functions.R 
if(sim=='contY'){
  expt_type <- c('noisy_only_predictor_2','noisy_linear_1_r_less','noisy_multicollinear_cand1_r_less', 'noisy_polynomial_r_less')
  null.value=0
} else{
  expt_type <- c('treatment_only','noisy_linear','noisy_multicollinear', 'noisy_polynomial')
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
# Function that computes the metrics for selected candidate algorithms 
# winner is a data.frame that computes the proportion of times each candidate algorithm was selected for adjustment
#=====================================================
get.selection <- function(cand, this.var, this.form ){
  winner <- data.frame(matrix(0, nrow=1, ncol=length(cand)+1)) 
  colnames(winner) <- c('unadj',cand)
  winner['unadj'] <- sum(this.var==1 & this.form=='glm')
  winner['glm']<- sum(this.var!=1 & this.var!=-99 & this.form=='glm')
  winner['main']<- sum(this.var==-99 & this.form=='glm')
  winner['stepwise'] <- sum(this.form=='stepwise')
  winner['step.interaction'] <- sum(this.form=='step.interaction')
  winner['lasso'] <- sum(this.form=='lasso')
  winner['mars'] <- sum(this.form=='mars')
  winner
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
dgp <- c('Sad', 'Linear', 'Interactive', 'Polynomial')
WINNERQ <- WINNERG <- NULL


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
  
  winnerq <- cbind(expt=expt_type[j], stratify=STRATIFY[k], 
                  get.selection(cand=unique(AP.fancy$cand.Qform),
                    this.var= SELECT$QAdj,
                    this.form= SELECT$Qform)
  )
  winnerq[, 3:8] <- paste0(round(winnerq[,3:8]/nReps*100, 1), '%')
  WINNERQ <- rbind(WINNERQ, winnerq)
  
  
  winnerg <- cbind(expt=expt_type[j], stratify=STRATIFY[k],
                   get.selection(cand=unique(AP.fancy$cand.gform),
                                 this.var= SELECT$gAdj,
                                 this.form= SELECT$gform)
  )
  winnerg[, 3:8] <- paste0(round(winnerg[,3:8]/nReps*100, 1), '%')
  WINNERG <- rbind(WINNERG, winnerg)
  
  # Create the data frame that stores all the 95%CI metrics
  data <- data.frame(
    x=c(1: (nReps*4)), 
    value1=c(UNADJ[["CI.lo"]],FORCE[["CI.lo"]],SIMPLE[["CI.lo"]],FANCY[["CI.lo"]]), 
    value2=c(UNADJ[["CI.hi"]],FORCE[["CI.hi"]],SIMPLE[["CI.hi"]],FANCY[["CI.hi"]]),
    ests=c(rep('Unadjusted',nReps), rep('Static',nReps), rep('Small APS', nReps), rep('Large APS', nReps))
  )

  if(effect){
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
    # ggsave(filename = paste0("PLOTS/", file.name, '.eps'))
    ggsave(filename = paste0("PLOTS/", file.name, '.png'))
  }

  
  rm(yay, winnerq, winnerg, data)
}
}

#=====================================================
# Create the tables for the metrics for inputs specified 
#=====================================================
library(xtable)

# Generate table for latex version 
YAY$DGP <- c('Sad', rep('', 7),  'Linear', rep('', 7), 'Interactive', rep('', 7), 
             'Polynomial', rep('', 7) )
YAY$Design <- rep( c('Simple', '','','', 'Stratified', '','',''), 4)
this.order <- c('DGP','Design','ests', 'cover', 'power', 'mse', 'bias', 'var', 're')
YAY[,this.order]

print(xtable(YAY[,this.order], digits=c(1, 1,1,1, rep(3, 6) )), include.rownames=FALSE)

#=====================================================
# Print savings obtained while using APS compared to unadjusted estimator
#=====================================================
# drop sad
if(sim=='contY'){
  x <- YAY[YAY$expt!='noisy_only_predictor_2',]
}else{
  x <- YAY[YAY$expt!='treatment_only',]
  
}
round( summary(x[x$ests=='Large APS', 're']), 3)
round( summary(x[x$ests=='Large APS',  'savings'])*100, 0)
round( summary(x[x$ests=='Small APS', 're']), 3)
round( summary(x[x$ests=='Small APS',  'savings'])*100, 0)


if(effect){
  #rm(dd)
  dd <- x[,c('expt','stratify', 'ests', 'savings')]

  dd <- dd[dd$ests!='Unadjusted',]
  dd$savings <- dd$savings*100
  
  dd$DGP <- c( rep( 'Linear',6), rep('Interactive', 6), rep('Polynomial',6) )
  dd$Value <- round(dd$savings)
  colnames(dd) <- c('expt', 'Stratify', 'Estimator',  'savings', 'DGP','Value')
  
  dd$DGP <- factor(dd$DGP, levels=c('Linear', 'Interactive', 'Polynomial'))

  dd$Estimator <- factor(dd$Estimator, levels=c('Static', 'Small APS', 'Large APS') )
  dd$savings <- as.numeric(dd$savings)
  dd$Stratify2 <- 'Simple'
  dd[dd$Stratify,'Stratify2'] <- 'Stratified'
  text.size <- 20
  these.colors <- c('#FFC107', '#1E88E5', '#004D40')
  text.color <-'black'
  adder <-  5
  
  this.legend.position <- 'bottom'
  
  g <- ggplot(dd, aes(fill=Estimator, y=savings, x=DGP)) + 
    geom_bar(position="dodge", stat="identity") + 
    labs(
      y = 'Estimated Sample Size Savings (%)',
      x = element_blank() 
    ) +
    scale_fill_manual(values=these.colors) +
    theme_classic() +
    theme(  plot.title = element_blank(),
            legend.title=element_blank(), 
            text = element_text(size = text.size, face="bold"),
            axis.text=element_text(size=text.size, face="bold"),
            legend.text=element_text(size=text.size, face="bold"),
            #  legend.position='',
            legend.position = this.legend.position) +
    facet_wrap(~ Stratify2)

  g <-  g+ geom_text(aes(y =savings+1, label =Value),
                     col = text.color, size = 6,  # fontface = "bold",
                     position=position_dodge(.9))
  
  
  g
  file.name <- paste0('PLOTS/',sim,'savings', '.png')
  ggsave(file.name)
  
  
}


#=====================================================
# Find propoortion of times when different candidate algorithms where chosen
# Print the tables for Outcome and PScore
#=====================================================
xtable(WINNERQ)
xtable(WINNERG)