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
# 3. Number of folds in cross validation, 
if(n==500){
  V<- 5
  folder <- 'OUTPUT/'
} else if(n==100){
  V <- 10
  folder <- 'OUTPUT_other/'
} else if(n==40){
  V <- 20
  folder <- 'OUTPUT_other/'
}
# 4. Specify whether MARS was included as a candidate
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
  if(!effect) {
    # include setting with no prognostic variables
    expt_type <- c('noisy_only_predictor_2', 'noisy_linear_1_r_less','noisy_multicollinear_cand1_r_less', 'noisy_polynomial_r_less')
  }else{
    expt_type <- c('noisy_linear_1_r_less','noisy_multicollinear_cand1_r_less', 'noisy_polynomial_r_less')
  }
  null.value=0
} else{
  if(!effect){
    # include setting with no prognostic variables
    expt_type <- c('treatment_only', 'noisy_linear','noisy_multicollinear', 'noisy_polynomial')
  } else{
    expt_type <- c('noisy_linear','noisy_multicollinear', 'noisy_polynomial')
  }
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
get.selection <- function(this.var, this.form, outcome.reg ){
  all.cand <- c('U','W1','W2','W3','W4','W5')
  cand <- c('Unadj','GLM', 'Main', 'Step', 'StepInt',
            'LASSO', 'MARS')
  winner <- data.frame(matrix(0, nrow=1, ncol=length(cand))) 
  colnames(winner) <- cand
  winner['Unadj'] <- sum(this.var==1 & this.form=='glm')
  winner['GLM']<- sum(this.var!=1 & this.var<length(all.cand) & this.form=='glm')
  winner['Main']<- sum(this.var>length(all.cand) & this.form=='glm')
  winner['Step'] <- sum(this.form=='stepwise')
  winner['StepInt'] <- sum(this.form=='step.interaction')
  winner['LASSO'] <- sum(this.form=='lasso')
  winner['MARS'] <- sum(this.form=='mars')
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
if(!effect){
  dgp <- c('Txt only','Linear', 'Interactive', 'Polynomial')
} else{
  dgp <- c('Linear', 'Interactive', 'Polynomial')
}
WINNERQ <- WINNERG <- psi.ave <- power.gain <-  NULL


#=====================================================
# Load the output files that are produced by running Main.R and compute the metrics 
#=====================================================
for(j in 1:length(expt_type)){
  for(k in 1:2){
  file.name <- paste( sim, paste0('Effect', effect),
                      paste0('N', n), paste0('V',V), paste0('mars', incl.mars),
                      paste0('nReps', nReps),paste0('stratify', STRATIFY[k]),
                      paste0('type', expt_type[j]), sep = "_")
  
  file.nameD <- paste( folder, file.name, paste('.RData'), sep = "_")
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
  psi.ave <- c(psi.ave, round(mean(UNADJ$psi),3))
  if(effect){ #skip under the null
    power.gain <- c(power.gain, 
                    round( (yay[yay$ests=='Large APS', 'power'] - 
                              yay[yay$ests=='Unadjusted', 'power'])*100))
  }

  YAY <- rbind(YAY, yay)
  
  winnerq <- cbind(expt=expt_type[j], stratify=STRATIFY[k], 
                  get.selection(this.var= SELECT$QAdj, this.form= SELECT$Qform)
  )
  winnerq[, 3:ncol(winnerq)] <- paste0(round(winnerq[,3:ncol(winnerq)]/nReps*100, 1), '%')
  WINNERQ <- rbind(WINNERQ, winnerq)
  
  
  winnerg <- cbind(expt=expt_type[j], stratify=STRATIFY[k],
                   get.selection(this.var= SELECT$gAdj,
                                 this.form= SELECT$gform)
  )
  winnerg[, 3:ncol(winnerg)] <- paste0(round(winnerg[,3:ncol(winnerg)]/nReps*100, 1), '%')
  WINNERG <- rbind(WINNERG, winnerg)
  

  # if(effect & n==500){
  #   # Create the data frame that stores all the 95%CI metrics
  #   data <- data.frame(
  #     x=c(1: (nReps*4)), 
  #     value1=c(UNADJ[["CI.lo"]],FORCE[["CI.lo"]],SIMPLE[["CI.lo"]],FANCY[["CI.lo"]]), 
  #     value2=c(UNADJ[["CI.hi"]],FORCE[["CI.hi"]],SIMPLE[["CI.hi"]],FANCY[["CI.hi"]]),
  #     ests=c(rep('Unadjusted',nReps), rep('Static',nReps), rep('Small APS', nReps), rep('Large APS', nReps))
  #   )
  #   
  #   psi <- mean(UNADJ$psi)
  #   # Plot
  #   ggplot(data) +
  #     geom_segment( aes(x=x, xend=x, y=value1, yend=value2), color="grey") +
  #     #geom_point( aes(x=x, y=value1), color=rgb(0.2,0.7,0.1,0.5), size=3 ) +
  #     #geom_point( aes(x=x, y=value2), color=rgb(0.7,0.2,0.1,0.5), size=3 ) +
  #     geom_point( aes(x=x, y=value1, color=factor(ests)), size=1 ) +
  #     geom_point( aes(x=x, y=value2, color=factor(ests)), size=1 ) +
  #     geom_hline(yintercept=null.value, linetype='dashed', color='black', size=0.5) +
  #     coord_flip()+
  #     theme_ipsum() +
  #     labs(title = dgp[j])+
  #     theme(
  #       legend.position = "none",
  #       #plot.title = element_blank(),
  #       axis.title.x = element_blank(),
  #       axis.title.y = element_blank())
  #   # ggsave(filename = paste0("PLOTS/", file.name, '.eps'))
  #   ggsave(filename = paste0("PLOTS/", file.name, '.png'))
  # }

  
  rm(yay, winnerq, winnerg)
}
}
psi.ave
power.gain
summary(power.gain)
round ( summary( YAY[YAY$ests=='Large APS', 're'] ), 2) 
 1- round ( summary( YAY[YAY$ests=='Large APS', 're'] ), 2) 
round ( summary( YAY[ YAY$ests=='Small APS', 're'] ), 2) 

#=====================================================
# Create the tables for the metrics for inputs specified 
#=====================================================
library(xtable)

# Generate table for latex version 
DGP <- NULL
for(j in dgp){
  DGP <- c(DGP, c(j, rep('',7)))
}
YAY$DGP <- DGP
YAY$Design <- rep( c('Simple', '','','', 'Stratified', '','',''), length(dgp))
this.order <- c('DGP','Design','ests', 'cover', 'power', 'mse', 'bias', 'var', 're')
YAY[,this.order]

print(xtable(YAY[,this.order], digits=c(1, 1,1,1, rep(3, 6) )), include.rownames=FALSE)


#=====================================================
# Find proportion of times when different candidate algorithms where chosen
# Print the tables for Outcome and PScore
#=====================================================
factorize_me <- function(dd,  strata.col='stratify'){
  this.rep <- nrow(dd)/3
  
  DGP <- c( rep( 'Linear',this.rep), rep('Interactive', this.rep), rep('Polynomial',this.rep) )
  dd<- cbind(DGP, Design='Simple',dd)
  dd[dd[,strata.col],'Design'] <- 'Stratified'
  dd <- subset(dd, select=-c(expt, stratify))
  dd
}


# ESTIMATED SAMPLE SIZE SAVINGS GRAPH
if(effect & n==500 ){
  
  WINNERQ <- cbind( Target=c('OutReg', rep('', (nrow(WINNERQ)-1) )), 
                    factorize_me(WINNERQ))
  WINNERG <- cbind( Target=c('PScore', rep('', (nrow(WINNERG)-1) )), 
                    factorize_me(WINNERG))
  WINNERG[,'StepInt'] <- WINNERG[,'MARS'] <- '-'
  print(xtable(rbind(WINNERQ,WINNERG)), include.rownames=F)
  
  #rm(dd)
  dd <- YAY[,c('expt','stratify', 'ests', 'savings')]
  dd <- dd[dd$ests!='Unadjusted',]
  
  dd <- factorize_me(dd)
  dd$savings <- dd$savings*100
  dd$Value <- round(dd$savings)

  dd$ests <- factor(dd$ests, levels=c('Static', 'Small APS', 'Large APS') )
  dd$savings <- as.numeric(dd$savings)
  text.size <- 16
  these.colors <- c('#bdd7e7','#3182bd', '#08519c')
  text.color <-'black'
  adder <-  5
  
  this.legend.position <- 'bottom'
  if(sim_flag){
    ylab <- 'Estimated Sample Size Savings (%) - Continous Outcome'
  }else{
    ylab <- 'Estimated Sample Size Savings (%) - Binary Outcome'
  }

  g <- ggplot(dd, aes(fill=ests, y=savings, x=DGP)) + 
    geom_bar(position="dodge", stat="identity") + 
    labs(
      y = ylab,
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
    facet_wrap(~ Design)

  g <-  g+ geom_text(aes(y =savings+1, label =Value),
                     col = text.color, size = 6,  # fontface = "bold",
                     position=position_dodge(.9))
  
  
  g
  file.name <- paste0('PLOTS/',sim,'savings', '.eps')
  ggsave(file.name, w=10, h=8)
  
  
}


# save(YAY, file=paste0('Summary_', sim, paste0('Effect', effect),
#                             paste0('N', n),'.Rdata') )


# if(effect ){
#   MAIN <- YAY
#   # load in the null
#   load(paste0('Summary_', sim, paste0('Effect', F),
#                             paste0('N', n),'.Rdata') )
#   
#   print(xtable(cbind(MAIN[,this.order], YAY$power), 
#                      digits=c(1, 1,1,1, rep(2, 7) )), include.rownames=FALSE)
#   
# }